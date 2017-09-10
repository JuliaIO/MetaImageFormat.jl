# MetaImageFormat

[![Build Status](https://travis-ci.org/JuliaIO/MetaImageFormat.jl.svg?branch=master)](https://travis-ci.org/JuliaIO/MetaImageFormat.jl)

[![codecov.io](http://codecov.io/github/JuliaIO/MetaImageFormat.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaIO/MetaImageFormat.jl?branch=master)

This package supports the
[MetaImage](https://itk.org/Wiki/ITK/MetaIO/Documentation#Reference:_Tags_of_MetaImage)
file format for the Julia language.  You can read "*.mhd" files using

```julia
using FileIO
img = load("myimage.mhd")
```

This package is a work-in-progress, and writing is not yet supported.

This implementation was started by copying the
[NRRD](https://github.com/JuliaIO/NRRD.jl) package and was based
purely on existing documentation on the MetaImage format.
Documentation on the format's definition seems somewhat sparse (e.g.,
which fields are required and which are optional? In what order should
they appear in the header? Are certain redundant combinations allowed
and other disallowed? ...).  In practice, the format appears to be
largely defined by the IO capabilities of ITK and Fiji/ImageJ. In
cases where there might be disagreements, one should check the source
code of these other projects.
