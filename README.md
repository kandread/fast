# Flood Assessment Spatial Triage

The code is contained in the `flooding.jl` file but you have to create the [Julia](https://julialang.org/) environment beforehand

```julia
import Pkg
Pkg.activate(".")
Pkg.instantiate()
```

After that you can generate the flood inundation files using

```julia
using Glob
include("flooding.jl")
files = readdir(glob"FloodFiles/*txt")
for f in files
  flooded(f, "merit.tif")
  println(f)
end
  ```
