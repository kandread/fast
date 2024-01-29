using CSV, DataFrames, DataStructures
using Rasters, ArchGDAL

"""
Read GNSS-R Rongowai observations from CSV `filename`.
"""
function readRongowai(filename::String)
    obs = CSV.File(filename;
                   select=[1, 2, 5, 6, 7, 8, 10],
                   header=["Lat", "Lon", "", "", "Water", "Fmajor", "Fminor", "Forient", "", "Angle", ""])
    obs |> DataFrame
end

"""
Estimated flooded area for Rongowai flight.

# Arguments
- `filename`: CSV file with Rongowai GNSS-R data
- `dem`: Digital Elevation Model (DEM) raster

"""
function flooded(filename::String, demfile::String)
    obs = readRongowai(filename)
    dem = Raster(demfile, missingval=-9999)
    bufdist = 0.1 # distance to buffer bounding box for GNSS-R data
    minlat = minimum(obs.Lat) - bufdist
    maxlat = maximum(obs.Lat) + bufdist
    minlon = minimum(obs.Lon) - bufdist
    maxlon = maximum(obs.Lon) + bufdist
    # extract DEM region with Rongowai flight bounding box
    ndem = setmappedcrs(view(dem, X(minlon .. maxlon), Y(minlat .. maxlat)), EPSG(4326))
    flood = bathtub(ndem, obs)
    write(replace(filename, ".txt" => ".tif"), flood, options=Dict("compress" => "LZW"))
end

"""
Initialize flood pixels from Rongowai GNSS-R observations.

# Arguments
- `dem`: DEM raster
- `obs`: DataFrame that contains GNSS-R observations

Uses an assumption of 1° equal to 110 km distance, and a planar projection of the ellipsoid GNSS-R footprint to intersect the Rongowai observation with the DEM.
See https://math.stackexchange.com/a/573174 for derivation.

"""
function initFlood(dem::Raster, obs::DataFrame)
    q = Deque{Tuple{Int, Int}}()
    nr, nc = size(dem)
    bnds = mappedbounds(dem)
    xmin, xmax = bnds[1]
    ymin, ymax = bnds[2]
    yres = (ymax - ymin) / nc
    xres = (xmax - xmin) / nr
    lats = mappedindex(dem)[2]
    lons = mappedindex(dem)[1]
    # initialize flood pixel queue
    for c=1:nrow(obs)
        if obs.Water[c] == 1
            j = floor(Int, (ymax - obs.Lat[c]) / yres)
            i = floor(Int, (obs.Lon[c] - xmin) / xres)
            push!(q, (i, j))
            # FIXME: can we use 110 km to convert from 1°?
            dj = floor(Int, (obs.Fminor[c] / 110000.) / yres)
            di = floor(Int, (obs.Fmajor[c] / 110000.) / xres)
            for k=-di:di
                for l=-dj:dj
                    ni = min(max(i+k, 1), nr)
                    nj = min(max(j+l, 1), nc)
                    # REVIEW: ellipsoid projection
                    if !((ni, nj) in q) && (3/4) * (lons[i] - lons[ni])^2 + (lats[j] - lats[nj])^2 < 1
                        push!(q, (ni, nj))
                    end
                end
            end
        end
    end
    q
end

"""
Run a bathtub model to generate flooded area.

# Arguments
- `dem`: DEM raster
- `obs`: DataFrame that contains GNSS-R observations

"""
function bathtub(dem::Raster, obs::DataFrame)
    out = zeros(Int, size(dem))
    nr, nc = size(dem)
    q = initFlood(dem, obs)
    # iterate through queue and flood downstream pixels
    while !isempty(q)
        i, j = pop!(q)
        out[i, j] = 1
        dh = 1.0
        ci, cj = i, j
        while dh > 0
            for k=-1:1
                for l=-1:1
                    ni = min(max(ci+k, 1), nr)
                    nj = min(max(cj+l, 1), nc)
                    if dem[i, j] > dem[ni, nj] && out[ni, nj] == 0 && dem[ni, nj] != missingval(dem)
                        push!(q, (ni, nj))
                    end
                end
            end
            i1, i2 = max(ci-1, 1), min(ci+1, nr)
            j1, j2 = max(cj-1, 1), min(cj+1, nc)
            dh = maximum(dem[ci, cj] .- dem[i1:i2, j1:j2])
            n = argmax(dem[ci, cj] .- dem[i1:i2, j1:j2])
            ci = collect(i1:i2)[n[1]]
            cj = collect(j1:j2)[n[2]]
        end
    end
    Raster(out, dims(dem), missingval=0)
end
