using MetaImageFormat, FileIO, Images, AxisArrays
using Test

# TODO: change the next to lines to `load("example.mhd")` after registry
f = File(format"MetaImage", "example.mhd")
img = MetaImageFormat.load(f)

@test pixelspacing(img) == (25, 20, 40)
@test axisnames(img) == (:R, :A, :I)
@test eltype(img) == N0f16
@test size(img) == (8, 9, 5)

@info "Warning \"AnatomicalOrientation string unrecognized, got ???. Using xyz instead.\" is expected"
f2 = File(format"MetaImage", "example_dup_axes.mhd")
img2 = MetaImageFormat.load(f2)
@test parent(img) == parent(img2)
@test axisnames(img2) == (:x, :y, :z)
