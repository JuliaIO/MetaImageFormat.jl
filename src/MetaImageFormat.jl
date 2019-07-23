module MetaImageFormat

# Packages needed to return the possible range of element types
using FixedPointNumbers, Colors, ColorVectorSpace
# Other packages
using AxisArrays, MappedArrays
using FileIO
# import Libz   # what compression does MetaImageFormat use? seems undocumented

using Colors: AbstractGray
using AxisArrays: HasAxes

using LinearAlgebra

string2type = Dict(
    "MET_UCHAR" => UInt8,
    "MET_SHORT" => Int16,
    "MET_USHORT" => UInt16,
    "MET_INT" => Int32,
    "MET_UINT" => UInt32,
    "MET_FLOAT" => Float32,
    "MET_UCHAR_ARRAY" => RGB{N0f8},
    "MET_USHORT_ARRAY" => RGB{N0f16},
)

# The opposite of string2type
type2string(::Type{UInt8})  = "MET_UCHAR"
type2string(::Type{Int16})  = "MET_SHORT"
type2string(::Type{UInt16}) = "MET_USHORT"
type2string(::Type{Int32})  = "MET_INT"
type2string(::Type{UInt32}) = "MET_UINT"
type2string(::Type{Float32}) = "MET_FLOAT"
type2string(::Type{RGB{T}}) where T = type2string(T)*"_ARRAY"

function space2axes(str)
    syms = Symbol.(tuple(str...))
    if length(syms) != 3 || !allunique(syms)
        @warn "AnatomicalOrientation string unrecognized, got $str. Using xyz instead."
        (:x, :y, :z)
    else
        syms
    end
end

axes2space(syms::NTuple{3,Symbol}) = join(syms)

# Parsing
VTuple{T} = Tuple{Vararg{T}}          # space-delimited tuple: 80 150

# This should list anything that DOESN'T parse to a string
# We replace True/False with the corresponding Bool
const parse_type = Dict(
    "NDims"=>Int,
    "ID"=>Int,
    "ParentID"=>Int,
    "BinaryData"=>Bool,
    "CompressedData"=>Bool,
    "BinaryDataByteOrderMSB"=>Bool,
    "ElementByteOrderMSB"=>Bool,
    "Color"=>RGBA{Float32},
    "Origin"=>VTuple{Float32},
    "Position"=>VTuple{Float32},
    "Offset"=>VTuple{Float32},
    "Orientation"=>Matrix{Float32},
    "Rotation"=>Matrix{Float32},
    "TransformMatrix"=>Matrix{Float32},
    "ElementSpacing"=>VTuple{Float32},
    "DimSize"=>VTuple{Int},
    "HeaderSize"=>Int,
    "SequenceID"=>VTuple{Int},
    "ElementMin"=>Float32,
    "ElementMax"=>Float32,
    "ElementNumberOfChannels"=>Int,
    "ElementSize"=>VTuple{Float32},
)

# const fieldorder = ["ObjectType", "NDims", "BinaryData",
#                     "BinaryDataByteOrderMSB", "ElementByteOrderMSB",
#                     "CompressedData",
#                     "TransformMatrix", "Offset", "CenterOfRotation",
#                     "sizes", "spacings", "space directions", "kinds",
#                     "centers", "centerings", "thickness",
#                     "axis mins", "axismins", "axis maxs", "axismaxs",
#                     "labels", "units",
#                     "min", "max", "old min", "oldmin", "old max", "oldmax",
#                     "block size", "blocksize", "endian", "encoding",
#                     "space units", "space origin", "measurement frame",
#                     "line skip", "lineskip", "byte skip", "byteskip",
#                     "sample units", "sampleunits",
#                     "data file", "datafile"]

const mymsb = ENDIAN_BOM == 0x01020304

# Don't extend FileIO.load
# Set mode to "r+" if you want to be able to modify values in the
# image and have them update in the disk file
function load(f::File{format"MetaImage"}, args...; mode="r", mmap=:auto)
    open(f, mode) do io
        load(io, args...; mode=mode, mmap=mmap)
    end
end

function load(io::Stream{format"MetaImage"}, Tuser::Type=Any; mode="r", mmap=:auto)
    # Assemble all the information about the array we're about to
    # read: element type, size, and the "meaning" of axes
    header = parse_header(io)
    Traw, need_bswap = raw_eltype(header)
    szraw = (header["DimSize"]...,)  # "sizes" may change in outer_eltype!, grab it now
    T, nd, perm = outer_eltype!(header, Traw, Tuser)
    axs = get_axes(header, nd)
    sz = get_size(axs)

    # Read the data
    iodata = find_datafile(io, header; mode=mode)
    compressed = get(header, "CompressedData", false)
    if compressed
        error("The MetaImage format doesn't appear to specifiy the compression algorithm, so not supported")
        # iodata = Libz.ZlibInflateInputStream(iodata)
    end

    can_mmap = !compressed

    if mmap == true && (!can_mmap)
        error("Cannot use memory-mapped for reading a compressed file")
    end

    # Use memory-mapping for large files
    do_mmap = can_mmap && (prod(szraw) > 10^8) && (mmap == :auto)
    do_mmap |= can_mmap && (mmap == true)

    if !compressed
        szraw = checked_size(Traw, szraw, sz, iodata)
    end

    if do_mmap
        A = Mmap.mmap(iodata, Array{Traw,length(szraw)}, szraw, position(iodata);
                      grow=false)
        if need_bswap
            f = mode == "r+" ? (bswap, bswap) : bswap
            A = mappedarray(f, A)
        end
    else
        A = Array{Traw}(undef, szraw...)
        read!(iodata, A)

        if need_bswap
            A = [bswap(a) for a in A]
        end
    end

    if perm == ()
        if T != eltype(A)
            A = need_bswap ? A = mappedarray(x->T(x), A) : reshape(reinterpret(T, vec(A)), sz)
        end
    else
        A = permuteddimsview(A, perm)
        if T<:Color
            A = colorview(T, A)
        end
    end

    isa(axs, Dims) ? A : AxisArray(A, axs)
end

# function save(f::File{format"MetaImage"}, img::AbstractArray; props::Dict = Dict{String,Any}(), keyvals=nothing, comments=nothing)
#     open(f, "w") do io
#         write(io, magic(format"NRRD"))
#         save(io, img; props=props, keyvals=keyvals, comments=comments)
#     end
# end

# function save{T}(io::Stream{format"NRRD"}, img::AbstractArray{T}; props::Dict = Dict{String,Any}(), keyvals=nothing, comments=nothing)
#     axs = axisinfo(img)
#     header = headerinfo(T, axs)
#     header_eltype!(header, T)
#     # copy fields from props to override those in header
#     for (k, v) in props
#         header[k] = v
#     end
#     v = version(header, !(keyvals==nothing || isempty(keyvals)))
#     write_header(io, v, header, keyvals, comments)
#     datafilename = get(props, "datafile", "")
#     if isempty(datafilename)
#         datafilename = get(props, "data file", "")
#     end
#     if isempty(datafilename)
#         nrrd_write(io, img)
#     else
#         println(io, "data file: ", datafilename)
#         if !get(props, "headeronly", false)
#             open(datafilename, "w") do file
#                 nrrd_write(file, img)
#             end
#         end
#     end
# end

# axisinfo(img) = axisinfo(HasAxes(img), img)
# axisinfo(::HasAxes{true}, img) = axes(img)
# axisinfo(::HasAxes, img) = size(img)

### Interpreting header settings

"""
    arraytype!(header) -> T, axs, perm, need_bswap

Analyze the `header` dictionary to extract the element-type `T`, size
or axes information `axs`, the permutation `perm` (if any) that julia
should use for "wrapping" the read data, and a boolean `need_bswap`
indicating whether the data need to be byte-swapped (to account for
differences in endianness). `T` includes any color information (in
which case a dimension of the array will be "consumed"). `axs` will be
a Dims-tuple in simple cases, or an `Axes` tuple (from AxisArrays.jl)
if dimensions are labeled or have their spatial information
(pixelspacing, spacedirections, etc) specified. `perm` is the
permutation needed to move the color data to the first dimension, or
an empty tuple if no permutation is required.

This function may modify the `header` dictionary (the reason for the !
in the name), so make a copy first if necessary.
"""
function arraytype!(header)
    Traw, need_bswap = raw_eltype(header)
    T, nd, perm = outer_eltype!(header, Traw)
    axs = get_axes(header, nd)
    T, axs, perm, need_bswap
end

"""
    arraytype(filename)

Parse MetaImage header and call `arraytype!(header)`. See
`arraytype!` for information about the return values.
"""
function arraytype(filename)
    header = parse_header(filename)
    arraytype!(header)
end

# function headerinfo(T, axs)
#     header = Dict{String,Any}()
#     Traw = raw_eltype(T)
#     header["type"] = type2string(Traw)
#     header["endian"] = myendian()
#     header["encoding"] = "raw"
#     if T <: Gray
#         val = gray(one(T))
#         val = isa(val, FixedPoint) ? reinterpret(val) : val
#         header["sample units"] = string("gray ", val)
#     elseif T <: Union{RGB,RGBA}
#         val = red(one(T))
#         val = isa(val, FixedPoint) ? reinterpret(val) : val
#         valfmt = isa(val, Integer) ? 'd' : 'f'
#         if T <: RGB
#             colstr = "rgb"
#             valfmtstr = "(%$valfmt,%$valfmt,%$valfmt)"
#             vals = (val,val,val)
#         else
#             colstr = "rgba"
#             valfmtstr = "(%$valfmt,%$valfmt,%$valfmt,%$valfmt)"
#             vals = (val,val,val,val)
#         end
#         fmtstr = "%s $valfmtstr"
#         x = (colstr, vals...)
#         header["sample units"] = @eval @sprintf($fmtstr, $(x...))
#     end
#     # Do the axes information
#     header["dimension"] = length(axs)
#     if isa(axs, Base.Indices)
#         axs = map(length, axs)
#     end
#     if isa(axs, Dims)
#         header["sizes"] = [axs...]
#     else
#         # axs is an Axis-tuple
#         header["sizes"] = [map(length, axs)...]
#         axnames = map(ax->axisnames(ax)[1], axs)
#         isspace = map(s->!startswith(string(s), "time"), axnames)
#         if haskey(axes2space, axnames)
#             header["space"] = axes2space[axnames]
#         else
#             header["space dimension"] = sum(isspace)
#         end
#         header["kinds"] = [isspc ? "domain" : "time" for isspc in isspace]
#         if !all(isdefaultname, axnames)
#             header["labels"] = [string(s) for s in axnames]
#         end
#         rng = map(ax->axisvalues(ax)[1], axs)
#         stepval = map(step, rng)
#         unitstr = map(x->isa(x, Quantity) ? string(unit(x)) : "", stepval)
#         spacing = map(x->isa(x, Quantity) ? ustrip(x) : x,        stepval)
#         if !all(x->x=="", unitstr)
#             header["units"] = [unitstr...]
#         end
#         if !all(x->x==1, spacing)
#             header["spacings"] = [spacing...]
#         end
#         origin = map(x->isa(x, Quantity) ? ustrip(x) : x, map(first, rng))
#         if any(x->x!=0, origin)
#             header["space origin"] = [origin[[isspace...]]...]
#         end
#     end
#     # Adjust the axes for color
#     if T <: Colorant && !(T <: AbstractGray)
#         header["dimension"] = length(axs)+1
#         unshift!(header["sizes"], length(T))
#         if haskey(header, "spacings")
#             unshift!(header["spacings"], NaN)
#         end
#         if haskey(header, "labels")
#             unshift!(header["labels"], lowercase(string(T.name.name)))
#         end
#         if haskey(header, "units")
#             unshift!(header["units"], "")
#         end
#         if !haskey(header, "kinds")
#             header["kinds"] = ["domain" for d = 1:length(axs)]
#         end
#         if T <: RGB
#             colkind = "RGB-color"
#         elseif T <: HSV
#             colkind = "HSV-color"
#         elseif T <: XYZ
#             colkind = "XYZ-color"
#         elseif T <: RGBA
#             colkind = "RGBA-color"
#         else
#             colkind = string(length(T), "-color")
#         end
#         unshift!(header["kinds"], colkind)
#     end
#     header
# end

# function version(header, has_keyvalue::Bool=false)
#     for n = length(version_reqs):-1:1
#         vr = version_reqs[n]
#         for f in vr
#             if haskey(header, f)
#                 return n
#             end
#         end
#     end
#     has_keyvalue ? 2 : 1
# end

"""
    raw_eltype(header) -> Traw, need_bswap
    raw_eltype(::Type{T}) -> Traw

Get the "basic" element type of the data, e.g., `UInt16` or
`Float32`.

This function does not try to determine whether the image is color
(`Traw` does not contain any color information), nor does it try to
interpret `Traw` as a `Normed` type.

See also: outer_eltype!, fixedtype.
"""
function raw_eltype(header)
    Traw = string2type[header["ElementType"]]
    msb = haskey(header, "BinaryDataByteOrderMSB") ? header["BinaryDataByteOrderMSB"] :
          haskey(header, "ElementByteOrderMSB") ? header["ElementByteOrderMSB"] : false
    need_bswap = msb != mymsb && sizeof(Traw) > 1
    Traw, need_bswap
end

#raw_eltype{C<:Colorant}(::Type{C}) = raw_eltype(eltype(C))
raw_eltype(::Type{C}) where C <: Colorant = raw_eltype(eltype(C))
raw_eltype(::Type{T}) where T <:FixedPoint = FixedPointNumbers.rawtype(T)
raw_eltype(::Type{T}) where T = T

"""
    fixedtype(Traw, header) -> Tu

Attempt to interpret type `Traw` in terms of FixedPoint numbers.
If `Traw` cannot be interpreted as `Normed`, `Tu = Traw`.
"""
function fixedtype(::Type{Traw}, header) where Traw <: Unsigned
    Traw == UInt8 && return N0f8
    if Traw == UInt16
        maxval = UInt16(get(header, "ElementMax", typemax(UInt16)))
        maxval == 0xffff && return N0f16
        maxval == 0x3fff && return N2f14
        maxval == 0x3fff && return N2f14
        maxval == 0x0fff && return N4f12
        maxval == 0x03ff && return N6f10
        maxval == 0x00ff && return N8f8
        error("unknown maxval $maxval")
    end
    Traw
end
fixedtype(::Type{Traw}, header) where Traw = Traw

function fixedtype_max(::Type{Traw}, mx) where Traw <: Unsigned
    fmx = log2(mx+1)
    if round(fmx) == fmx
        return Normed{Traw,round(Int,fmx)}
    end
    Traw
end

# header_eltype!(header, ::Type) = header
# function header_eltype!{T<:FixedPoint}(header, ::Type{T})
#     header["sample units"] = string("gray ", reinterpret(one(T)))
#     header
# end
# function header_eltype!{C<:Colorant}(header, ::Type{C})
#     _header_eltype!(header, C, eltype(C))
#     header
# end
# function _header_eltype!{C<:AbstractGray,T<:FixedPoint}(header, ::Type{C}, ::Type{T})
#     header["sample units"] = string("gray ", reinterpret(one(T)))
# end
# function _header_eltype!{C<:AbstractGray,T}(header, ::Type{C}, ::Type{T})
#     header["sample units"] = string("gray ", one(T))
# end
# function _header_eltype!{C<:AbstractRGB,T<:FixedPoint}(header, ::Type{C}, ::Type{T})
#     o = reinterpret(one(T))
#     header["sample units"] = "rgb ($o,$o,$o)"
# end
# function _header_eltype!{C<:AbstractRGB,T}(header, ::Type{C}, ::Type{T})
#     o = one(T)
#     header["sample units"] = "rgb ($o,$o,$o)"
# end
# function _header_eltype!{C<:XYZ}(header, ::Type{C}, ::Type)
#     header["sample units"] = "xyz (95.047,100.000,108.883)"
# end
# function _header_eltype!{C<:HSV}(header, ::Type{C}, ::Type)
#     header["sample units"] = "hsv (360, 0, 1)"
# end

"""
    colorant_eltype(C, T) -> Tc

Return a valid "inner" element type `Tc` for colorant type `C`. When
`T` != `Tc`, values must be "converted" before they can be interpreted
as type `C`.
"""
colorant_eltype(::Type{C}, ::Type{T}) where {C <: Colorant, T <: AbstractFloat} = C{T}
colorant_eltype(::Type{C}, ::Type{T}) where {C <: Colorant, T} = C{Float32}

"""
    UnknownColor{T,N}

An unknown Color. This type gets returned when one of the "kind"
settings is "3-color" or "4-color".
"""
struct UnknownColor{T,N} <: Color{T,N}
    col::NTuple{N,T}
end

"""
    outer_eltype!(header, Traw) -> T, nd, perm

Extract the julia array `eltype` `T`, the number of dimensions `nd`
**excluding** color/complex/vector/matrix element data, and any
permutation needed to put the eltype dimension first. Any dimensions
in the header corresponding to color (or if "kind" is set to one of
the vector types) will be "consumed" upon exit. `Traw` is the
element type as determined by `raw_eltype`.

See also: raw_eltype.
"""
function outer_eltype!(header, Traw, Tuser=Any)
    nd = header["NDims"]
    sz = header["DimSize"]
    length(sz) == nd || error("parsing of sizes: $sz is inconsistent with $nd dimensions")
    perm = ()
    T = fixedtype(Traw, header)
    T, nd, perm
end

function get_axes(header, nd)
    if haskey(header, "AnatomicalOrientation")
        axnames = space2axes(header["AnatomicalOrientation"])
    else
        axnames = nd == 2 ? (:x, :y) :
                  nd == 3 ? (:x, :y, :z) : error("$nd dimension not handled")
    end
    tform = haskey(header, "TransformMatrix") ? header["TransformMatrix"] :
            haskey(header, "Rotation") ? header["Rotation"] : Matrix(I, nd, nd)
    if tform != Matrix(I, nd, nd)
        error("rotations aren't yet handled, got $tform")
    end
    offset = get(header, "Offset", zeros(nd))
    spacing = get(header, "ElementSpacing", ones(nd))
    sz = header["DimSize"]
    rng = map((o, s, l) -> range(o, step=s, length=l), offset, spacing, sz)
    (map((name, r) -> Axis{name}(r), axnames, rng)...,)
end

get_size(sz::Dims) = sz
get_size(axs) = map(length, axs)

### Parsing

"""
    header = parse_header(io)

Parse the MetaImage header (an .mhd file). `io` should be positioned
at the beginning of the header file.

See also: write_header.
"""
function parse_header(io)
    header = Dict{String, Any}()
    # Read until we encounter a blank line, which is the separator
    # between the header and data
    line = strip(readline(io))
    while !isempty(line)
        idx = findfirst(isequal('='), line)
        idx == 0 && error("no colon found in $line")
        key, value = strip(line[1:idx-1]), strip(line[idx+1:end])
        T = get(parse_type, key, String)
        header[key] = meta_parse(T, value)
        line = strip(readline(io))
    end
    header
end

parse_header(s::Stream{format"MetaImage"}) = parse_header(stream(s))

function parse_header(filename::AbstractString)
    f = File{format"MetaImage"}(filename)
    open(f) do io
        parse_header(io)
    end
end

# """
#     write_header(io, version, header, [keyvals, [comments]])

# Write an NRRD header, as the top of the .nrrd or the separate .nhdr
# file. `io` should be positioned just after the initial "NRRD" in the
# file. This writes the header and a blank line, so that at the end `io`
# is positioned at the first byte of the data (if present).

# Note that if you're writing a header for a "detached" data file
# (separate .nhdr and .raw files), `header` should contain a "data file"
# or "datafile" field storing the name of the .raw file.

# Inputs:
# - `version` is a 4-character string, e.g., "0002", giving the NRRD version of the header, or a integer corresponding to a recognized NRRD header version number
# - `header` is a `Dict{String,Any}` of `field=>setting` pairs (as returned by `parse_header`)
# - `keyvals` is a `Dict{String,String}` containing `key=>value` pairs (NRRD0002 or higher, lines like key:=value; many NRRD files do not contain any of these)
# - `comments` is an array containing lines of the header that began with `#` (but with the `#` and leading whitespace stripped out)

# See also: parse_header.
# """
# function write_header(io::IO, version, header, keyvals=nothing, comments=nothing)
#     writeversionstr(io, version)

#     for fn in fieldorder
#         if haskey(header, fn)
#             print(io, fn, ": ")
#             T = get(parse_type, fn, String)
#             nrrd_format(io, T, header[fn])
#             println(io, "")
#         end
#     end
#     if keyvals != nothing
#         for (k,v) in keyvals
#             println(io, k, ":=", v)
#         end
#     end
#     if comments != nothing
#         for c in comments
#             println(io, "# ", c)
#         end
#     end
#     println(io)
#     nothing
# end
# write_header(s::Stream{format"NRRD"}, args...) = write_header(stream(s), args...)

meta_parse(::Type{T}, str) where T = parse(T, str)  # fallback
meta_parse(::Type{String}, str) = str

function meta_parse(::Type{Bool}, str)
    str == "True" && return true
    str == "False" && return false
    error("unrecognized boolean string $str")
end

function meta_parse(::Type{VTuple{T}}, s::AbstractString) where T
    ss = split(s)
    v = Vector{T}(undef, length(ss))
    for i = 1:length(ss)
        v[i] = meta_parse(T, ss[i])
    end
    return v
end

function meta_parse(::Type{Matrix{T}}, s::AbstractString) where T
    v = meta_parse(VTuple{T}, s)
    l = length(v)
    sqrtl = round(Int, sqrt(l))
    sqrtl^2 == l || error("cannot interpret transformation matrix $v")
    reshape(v, sqrtl, sqrtl)
end

meta_format(io, ::Type{T}, x) where T = print(io, x)

function meta_format(io, ::Type{VTuple{T}}, container) where T
    isfirst = true
    for x in container
        if isfirst
            isfirst = false
        else
            print(io, ' ')
        end
        nrrd_format(io, T, x)
    end
    nothing
end

### Utilities

function chksize(sz, sztarget)
    sz == sztarget || error("dimension size should be $sztarget, got $sz")
    nothing
end

function spacings(spacedirections)
    spacings = zeros(length(spacedirections))
    can_convert = true
    d = 1
    for (i, v) in enumerate(spacedirections)
        if v == "none"
            spacings[i] = NaN
        else
            n = 0
            for (j,x) in enumerate(v)
                if x != 0
                    if j == d   # make certain they are in order
                        spacings[i] = x
                        d += 1
                    else
                        can_convert = false
                    end
                    n += 1
                end
            end
            can_convert &= n < 2
        end
    end
    can_convert, spacings
end

### File-related functions

# Adjust the array size if the file is not big enough for reading to succeed
function checked_size(Traw, szraw, sz, iodata)
    cpos = position(iodata)
    datalen = div(filesize(stat(iodata)) - cpos, sizeof(Traw))
    if datalen < prod(szraw)
        if szraw != sz
            # If we've dropped a dimension due to "absorbing"
            # color, etc, don't try to figure out how to correct
            # it, just punt to the user
            error("data length $datalen too small for size $sz array of $T")
        end
        # If the data are smaller than the header suggests, read in as
        # much "complete" data as are available
        sznew = [szraw...]
        strds = [1;cumprod(sznew)]
        k = length(sznew)
        sznew[k] = div(datalen, strds[k])
        while sznew[k] == 0 && k > 1
            pop!(sznew)
            k -= 1
            sznew[k] = div(datalen, strds[k])
        end
        tsznew = (sznew...,)
        warn("header indicates an array size $szraw, but the file size is consistent with at most $tsznew")
        szraw = tsznew
    end
    szraw
end

function stream2name(s::IO)
    name = s.name
    if !startswith(name, "<file ")
        error("io name ", name, " doesn't fit expected pattern")
    end
    name[7:end-1]
end

function find_datafile(iodata::IOStream, header; mode="r", path=nothing)
    fdata = header["ElementDataFile"]
    if path == nothing
        # Check current directory first
        if isfile(fdata)
            iodata = open(fdata)
        else
            path = dirname(stream2name(iodata))
            iodata = open(joinpath(path, fdata), mode)
        end
    else
        iodata = open(joinpath(path, fdata), mode)
    end
    iodata
end

function find_datafile(s::Stream{format"MetaImage"}, header; mode="r")
    find_datafile(stream(s), header; mode=mode)
end

# nrrd_write(io, A::AxisArray) = nrrd_write(io, A.data)
# nrrd_write(io, A::AbstractArray) = nrrd_write_elty(io, A)
# nrrd_write_elty{C<:Colorant}(io, A::AbstractArray{C}) = nrrd_write_elty(io, channelview(A))
# nrrd_write_elty{T<:Fixed}(io, A::AbstractArray{T}) = nrrd_write_elty(io, rawview(A))
# nrrd_write_elty(io, A::AbstractArray) = write(io, A)

end
