module VAWTools

include("other-pks-fixes.jl")

# General tools, maybe applicable more widely.
export read_agr, write_agr, read_xyn, inpoly, AGridded, Gridded, Gridded1d, Traj, map_onto_bands, smooth_vector,
    downsample, split_traj!, boxcar, boxcar_matrix, bin_grid, piecewiselinear, split_poly,
    transform_proj


import Base: ==, size, length, step, +, -, *

# const spatialorder = "yx"  # order of imagestorage, https://github.com/timholy/Images.jl#storage-order-and-changing-the-representation-of-images

abstract AGridded{T}
"Check whether an error-field is defined."
haserror(g::AGridded) = isdefined(g,:err)
size(g::AGridded) = size(g.v)
length(g::AGridded) = length(g.v)
step(g::AGridded) = step(g.x)

function Base.show(io::IO, g::AGridded)
    show(io, typeof(g))
    print(io, " ")
    # if haserror(g)
    #     println(io, " holding values and errors on grid:")
    # else
    #     println(io, " holding values on grid:")
    # end
    show(io, g.x);
    if isdefined(g, :y)
        print(io, " "); show(io, g.y)
    end
end


"""
Gridded holds 2D data (and its error) on a regular grid.

- rows of v correspond to x-dir, columns to ydir
- be sure to be clear whether x corresponds to cell centers or
  boundaries.


TODO:

- hold projection information?
- hold NA value?
"""
immutable Gridded{T} <: AGridded{T}
    x::FloatRange{Float64}
    y::FloatRange{Float64}
    midpoint::Bool  # if true, the (x,y) is cell midpoint, otherwise lower-left corner
    v::Matrix{T}  # values
    err::Matrix{T}  # error of values
    function Gridded(x,y,midpoint,v,err)
        nx, ny = length(x), length(y)
        @assert size(v)==size(err)
        @assert size(v)==(length(x), length(y))
        new(x,y,midpoint,v,err)
    end
    function Gridded(x,y,midpoint,v)
        @assert size(v)==(length(x), length(y))
        new(x,y,midpoint,v)
    end
    Gridded() = new()
end
Gridded{T}(x, y, mp, v::Matrix{T}) = Gridded{T}(x, y, mp, v)
Gridded{T}(x, y, mp, v::Matrix{T}, err::Matrix) = Gridded{T}(x,y,mp,v,err)
Gridded{T}(x, y, v::Matrix{T}) = Gridded{T}(x, y, true, v)
Gridded{T}(x, y, v::Matrix{T}, err::Matrix) = Gridded{T}(x,y,true,v,err)
function Base.convert{T2,T1}(::Type{Gridded{T2}}, g::Gridded{T1})
    if haserror(g)
        Gridded(g.x, g.y, g.midpoint, convert(Matrix{T2}, g.v), convert(Matrix{T2}, g.err))
    else
        Gridded(g.x, g.y, g.midpoint, convert(Matrix{T2}, g.v))
    end
end
function downsample(g::Gridded, step, start=1)
    if haserror(g)
        Gridded(g.x[start:step:end],
                g.y[start:step:end],
                g.midpoint,
                g.v[start:step:end, start:step:end],
                g.err[start:step:end, start:step:end])
    else
        Gridded(g.x[start:step:end],
                g.y[start:step:end],
                g.midpoint,
                g.v[start:step:end, start:step:end])
    end
end
# (+)(g1::Gridded, g2::Gridded) = (@assert g1.midpoint==g2.midpoint; Gridded(g1.x, g1.y, g1.midpoint, g1.v+g2.v))
# (*)(n::Number, g::Gridded) = Gridded(g.x, g.y, g.midpoint, g.v*n)
# (*)(g::Gridded, n::Number) = n*g


"""
Holds 1D fields, e.g. elevation band data.

Be sure to be clear whether x corresponds to cell centers or boundaries.
"""
immutable Gridded1d{T} <:  AGridded{T}
    x::FloatRange{Float64}
    midpoint::Bool  # if true, the x is cell midpoint, otherwise the lower bound
    v::Vector{T} # values
    err::Vector{T}  # error of values
    function Gridded1d(x,mp,v,err)
        @assert size(v)==size(err)
        @assert size(v)==(length(x), )
        new(x,mp,v,err)
    end
    function Gridded1d(x,mp,v)
        @assert size(v)==(length(x), )
        new(x,mp,v)
    end
    Gridded1d() = new()
end
Gridded1d{T}(x, mp, v::Vector{T}) = Gridded1d{T}(x,mp,v)
Gridded1d{T}(x, mp, v::Vector{T}, err::Vector) = Gridded1d{T}(x,mp,v,err)
Gridded1d{T}(x, v::Vector{T}) = Gridded1d{T}(x,true,v)
Gridded1d{T}(x, v::Vector{T}, err::Vector) = Gridded1d{T}(x,true,v,err)
(+)(g1::Gridded1d, g2::Gridded1d) = (@assert g1.midpoint==g2.midpoint; Gridded1d(g1.x, g1.midpoint, g1.v+g2.v))
(*)(n::Number, g::Gridded1d) = Gridded1d(g.x, g.midpoint, g.v*n)
(*)(g::Gridded1d, n::Number) = n*g


"""
Holds values on a trajectory or generally unstructured data.

Fields:
- x,y -- coordinates
- splits::Vector{UnitRange{Int}} -- If the trajectory consists of several sub-trajectories
- v::Vector{T} -- values
- err::Vector{T} -- errors (or second set of values)

There are constructors to leave off splits, v and err.
"""
type Traj{T}
    x::Vector{Float64}
    y::Vector{Float64}
    splits::Vector{UnitRange{Int}}
    v::Vector{T}
    err::Vector{T}
    function Traj(x,y,v,err,splits::Vector{UnitRange{Int}})
        @assert length(x)==length(y)==length(v)==length(err)
        new(x,y,splits,v,err)
    end
    function Traj(x,y,v,err)
        @assert length(x)==length(y)==length(v)==length(err)
        new(x,y,[1:length(x)],v,err)
    end
    function Traj(x,y,v,splits::Vector{UnitRange{Int}})
        @assert length(x)==length(y)==length(v)
        new(x,y,splits,v)
    end
    function Traj(x,y,v)
        @assert length(x)==length(y)==length(v)
        new(x,y,[1:length(x)],v)
    end
    function Traj(x,y,splits::Vector{UnitRange{Int}})
        @assert length(x)==length(y)
        new(x,y,splits)
    end
    function Traj(x,y)
        @assert length(x)==length(y)
        new(x,y,[1:length(x)])
    end
    Traj() = new()
end
Traj{T}(x, y, v::AbstractVector{T}, err::AbstractVector{T}) = Traj{T}(x,y,v,err)
Traj{T}(x, y, v::AbstractVector{T}) = Traj{T}(x,y,v)
Traj{T}(x::AbstractVector{T}, y::AbstractVector{T}) = Traj{T}(x,y)

Traj{T}(x, y, v::AbstractVector{T}, err::AbstractVector{T}, splits::Vector{UnitRange{Int}}) =
    Traj{T}(x,y,v,err,splits)
Traj{T}(x, y, v::AbstractVector{T}, splits::Vector{UnitRange{Int}}) = Traj{T}(x,y,v,splits)
Traj{T}(x::AbstractVector{T}, y::AbstractVector{T}, splits::Vector{UnitRange{Int}}) = Traj{T}(x,y,splits)


haserror(t::Traj) = isdefined(t, :err)
hasvalues(t::Traj) = isdefined(t, :v)
Base.length(t::Traj) = length(t.x)
# "Splits trajectory by inserting NaNs"
# function split_traj_nan!(t::Traj, dist)
#     # TODO this drops datapoints!
#     for i=1:length(t.x)-1
#         if (t.x[i]-t.x[i+1])^2+(t.y[i]-t.y[i+1])^2 > dist^2
#             t.v[i]=NaN
#             if haserror(t)
#                 t.err[i]=NaN
#             end
#         end
#     end
#     nothing
# end
"Splits trajectory into several"
function split_traj!{T<:Traj}(t::T, dist)
    ii = 1
    empty!(t.splits)
    for i=1:length(t.x)-1
        if (t.x[i]-t.x[i+1])^2+(t.y[i]-t.y[i+1])^2 > dist^2
            push!(t.splits,ii:i)
            ii = i+1
        end
    end
    push!(t.splits, ii:length(t.x))
    nothing
end


function Base.convert{T2,T1}(::Type{Traj{T2}}, g::Traj{T1})
    if typeof(g)==Traj{T2}
        return g
    end
    x = convert(Vector{T2}, g.x)
    y = convert(Vector{T2}, g.y)
    if hasvalues(g)
        if haserror(g)
            Traj(x, y, convert(Vector{T2}, g.v), convert(Vector{T2}, g.err), g.splits)
        else
            Traj(x, y, convert(Vector{T2}, g.v), g.splits)
        end
    else
        Traj(x, y, g.splits)
    end

end

"""
Holds a rectangular grid, mirrors the ASCII-grid .agr, .grid, .asc, .bin files.

Note that the indexing into the value matrix ins awkward, see Ref
below.  In general use Gridded.

Ref: https://en.wikipedia.org/wiki/Esri_grid
"""
immutable AGR{T} # Ascii GRid
    v::Matrix{T}    # values: orientation is awkward, as in the AGR file, see https://en.wikipedia.org/wiki/Esri_grid
    nc::Int64       # NCOLS TODO: remove those
    nr::Int64       # NROWS
    xll::T          # XLLCORNER
    yll::T          # YLLCORNER
    dx::T           # CELLSIZE
    NA::T           # NODATA
    extra_header::Vector{Float32} # the .bin files have space for
                                  # extra information in the header
    function AGR(va, nc, nr, xll, yll, dx, NA, extra_header)
        @assert size(va)==(nr,nc)
        @assert dx>=0
        new(va, nc, nr, xll, yll, dx, NA, extra_header)
    end
end
AGR{T}(va::Matrix{T}, nc, nr, xll, yll, dx, NA, extra_header) = AGR{T}(va, nc, nr, xll, yll, dx, NA, extra_header)
AGR{T}(va::Matrix{T}, nc, nr, xll, yll, dx, NA) = AGR{T}(va, nc, nr, xll, yll, dx, NA, zeros(Float32, 6))

"""
Transform Gridded to AGR

Optional:
- NA_g: no-data value of g:Gridded, defaults to NaN
- NA_agr: no-data value of AGR output, defaults to NaN
- write_err: write g.err instead of g.v
"""
function AGR{T}(g::Gridded{T}; NA_g=convert(T,NaN), NA_agr=convert(T,NaN), write_err=false)
    v = write_err ? g.err : g.v # whether to write the values or errors into the ascii grid
    if !isequal(NA_g,NA_agr)
        v = deepcopy(v)
        for i=eachindex(v)
            if isequal(v[i],NA_g)
                v[i] =  NA_agr
            end
        end
    end

    v = rotl90(v)
    nc,nr = size(g)
    dx = convert(T, step(g.x))
    if g.midpoint
        xll = g.x[1]-dx/2
        yll = g.y[1]-dx/2
    else
        xll = g.x[1]
        yll = g.y[1]
    end
    AGR{T}(v, nc, nr, xll, yll, dx, NA_agr, zeros(Float32, 6))
end

# xs(g::AGR) = linspace(g.xll+g.dx/2, nr*g.dx, nr)
# ys(g::AGR) = linspace(g.yll+g.dx/2, nc*g.dx, nc)

size(g::AGR) = size(g.v)
length(g::AGR) = length(g.v)
function ==(g1::AGR,g2::AGR)
    g1.v==g2.v &&
    g1.nc==g2.nc &&
    g1.nr==g2.nr &&
    g1.xll==g2.xll &&
    g1.yll==g2.yll &&
    g1.dx==g2.dx &&
    g1.NA==g2.NA &&
    g1.extra_header==g2.extra_header
end

# function my_rotr90(A::AbstractMatrix)
#     ind1, ind2 = indices(A)
#     B = similar(A, (ind2,ind1))
#     m = first(ind1)+last(ind1)
#     @inbounds for i=ind1, j=indices(A,2) #added inbounds
#         B[j,m-i] = A[i,j]
#     end
#     return B
# end


"""
Construct Gridded from AGR.  Replace AGR-no-value NA with NaN.  By
default assumes AGR-no-value==NaN.  Moves (x,y)-coords to midpoints
"""
function Gridded{T}(agr::AGR{T}; NA=convert(T,NaN))
    #    v = my_rotr90(agr.v)
    v = rotr90(agr.v)
    if !isequal(NA, agr.NA)
        for i=eachindex(v)
            if isequal(v[i], agr.NA); v[i] = NA end
        end
    end
    Gridded(range(agr.xll+agr.dx/2, agr.dx, agr.nc),
            range(agr.yll+agr.dx/2, agr.dx, agr.nr),
            true,
            v)
end

# File readers
##############

"""Check if .bin file nor not"""
isbin_file(fn::AbstractString) = splitext(fn)[2]==".bin"

"""Read a Ascii grid.
https://en.wikipedia.org/wiki/Esri_grid

In:
fn -- file name
Optional keywords:
T -- type of output array.  Defaults to Float32.
NA -- replace the fill value with this.  Defaults to use NaN.

Out:
Gridded
"""
function read_agr(fn::AbstractString, T=Float32; NA=convert(T,NaN))
    Gridded(_read_agr(fn::AbstractString, T, NA), NA=NA)::Gridded{T}
end

function _read_agr(fn::AbstractString, T=Float32, NA=nothing)
    if !isfile(fn)
        error("File $fn does not exist!")
    end
    if isbin_file(fn)
        toT = (io,T) -> convert(T, read(io, Float32))
    else
        toT = (io,T) -> parse(T, split(readline(io))[2])
    end
    open(fn, "r") do io
        local va::Matrix{T}
        # read header
        # nc = toT(io)::Int
        # nr = toT(io)::Int
        # xll = toT(io,T)::T
        # yll = toT(io,T)::T
        # dx = toT(io,T)::T
        # fill = toT(io,T)::T
        nc = toT(io, Int)
        nr = toT(io, Int)
        xll = toT(io,T)
        yll = toT(io,T)
        dx = toT(io,T)
        fill = toT(io,T)
        if isbin_file(fn)
            # read extra header
            extra_header = read(io, Float32, 6)
            # read values
            va = read(io, Float32, nc, nr).'
            if eltype(va)!=T
                error("Not implemented yet")
            end
            if !eof(io)
                warn("End-of-file was not reached!")
            end
        else
            va = Array(T, nr, nc)
            # no extra header for ascii .agr
            extra_header = zeros(Float32, 6)
            # read values
            tmp = split(readstring(io))
            if length(tmp)!=nc*nr
                error("Something's wrong with the file.  It should contain $(nc*nr) values but only has $(length(tmp))")
            end
            for i=1:nr, j=1:nc
                va[i,j] = parse(T, tmp[(i-1)*nc + j])
            end
        end
        # replace missing value with something else
        if NA!=nothing && fill!=NA
            fill = _refill!(va, fill, NA)
        end
        # make a AGR data structure
        AGR(va, nc, nr, xll, yll, dx, fill, extra_header)
    end
end
# A function barrier is needed, thus this helper function:
function _refill!(a::Array, oldfill, newfill)
    @inbounds for i in eachindex(a)
        if a[i]==oldfill
            a[i] = newfill
        end
    end
    newfill
end


"""Write Ascii grid.

Output format is determined by the extension:
- .bin -- binary
- .agr -- ASCII
- .grid or .asc are renamed to .agr

Always writes .bin in Float32

Optional:
- NA_g: no-data value of input g::Gridded, defaults to NaN
- NA_agr: no-data value of AGR-write output, defaults to NaN
- write_err: write g.err instead of g.v
- T: choose output type: Int or Float (only for ASCII files)
"""
function write_agr{T_}(g::Gridded{T_}, fn::AbstractString; T=T_, NA_g=convert(T_,NaN), NA_agr=convert(T,NaN), write_err=false)
    write_agr(AGR(g, NA_g=NA_g, NA_agr=NA_agr, write_err=write_err), fn, T=T)
end

# setting NA will transform the NA value to that
function write_agr{T_}(g::AGR{T_}, fn::AbstractString; NA=nothing, T=T_)
    if !isbin_file(fn)
        ext = splitext(fn)[2]
        if ext==".grid" || ext==".asc"
            println("Changing extension to .agr")
            ext = ".agr"
        elseif ext!=".agr"
            error("unsupported extension")
        end
    end
    if isbin_file(fn)
        convInt(x) = convert(Float32, x)
        convT(x) = convert(Float32, x)
        convMT = convT
    else # conversion to ASCII
        convInt = (x) -> @sprintf("%d\n", x)
        if T<:AbstractFloat
            convT = (x) ->  @sprintf("%f\n", x)
            convMT = (x) ->  @sprintf("%f  ", x)
        elseif T<:Integer
            convT = (x) ->  @sprintf("%d\n", x)
            convMT = (x) ->  @sprintf("%d  ", x)
        end
    end
    open(fn, "w") do io
        # write header
        isbin_file(fn) || write(io, "ncols         ")
        write(io, convInt(g.nc))
        isbin_file(fn) || write(io, "nrows         ")
        write(io, convInt(g.nr))
        isbin_file(fn) || write(io, "xllcorner     ")
        write(io, convT(g.xll))
        isbin_file(fn) || write(io, "yllcorner     ")
        write(io, convT(g.yll))
        isbin_file(fn) || write(io, "cellsize      ")
        write(io, convT(g.dx))
        isbin_file(fn) || write(io, "NODATA_value  ")
        # replace NA value if necessary
        if NA==nothing
            fill = g.NA
            va = g.v
        else
            fill = NA
            va = deepcopy(g.v)
            for i=1:length(va)
                if va[i]==g.NA
                    va[i] = fill
                end
            end
        end
        write(io, convT(fill))
        # write extra header if .bin:
        isbin_file(fn) && write(io, g.extra_header)
        # write values
        for i=1:g.nr, j=1:g.nc
            v = va[i,j]
            write(io, convMT(v))
        end
        return nothing
    end
end

"""
Read .xyn files which can contain several, joined polygons.

x              y           label
-2031744.122   833011.310  21
...

where the label marks the beginning or end of line (21 --start, 22 --
inside, 23 -- end).

There can be embedded polygons, with label sequence 23,20,21.  The
line with 20 is then dropped.

The polygon must be closed, i.e. first element equals last.

Set `hasz` if it also contains a z-coordinate, i.e. a .xyzn file.

Return:

- a list of [x,y(,z)] arrays, one for each polygon.
"""
function read_xyn(fn; hasz=false, fix=false)
    if !isfile(fn)
        error("File $fn cannot be found.")
    end
    x,y,z,l = open(fn, "r") do io
        x = Float64[]
        y = Float64[]
        z = Float64[]
        l = Int[]
        for ls in readlines(io)
            if hasz
                x_,y_,z_,l_ = split(ls)
            else
                x_,y_,l_ = split(ls)
            end
            push!(x, parse(Float64, x_))
            push!(y, parse(Float64, y_))
            if hasz
                push!(z, parse(Float64, z_))
            end
            push!(l, parse(Int, l_))
        end
        return x,y,z,l
    end
    # Split up into individual polygons
    n = length(x)
    out = Matrix{Float64}[]
    is = 1 # start index of one poly
    iend = -1  # end index of one poly
    i = 1
    drop20 = false
    while i<=n
        if l[i]!=21
            error("Malformed file: Expected 21 on line $i.")
        end
        # start a new polygon
        i += 1
        while true  # go through one poly:
            if l[i]==23 # end of a poly
                if i!=n && l[i+1]==20 # drop the 20-line
                    drop20 = true
                end
                break
            elseif l[i]==22
                i+=1
            else
                error("shouldn't get here")
            end
        end
        iend = i
        if length(is:iend)<3
            error("Polygon needs at least three vertices")
        end
        # fill:
        if hasz
            push!(out, (hcat(x[is:iend],y[is:iend],z[is:iend]))' )
        else
            push!(out, (hcat(x[is:iend],y[is:iend]))' )
        end
        if x[is]!=x[iend] || y[is]!=y[iend]
            if fix
                out[end] = hcat(out[end], out[end][:,1])
            else
                error("Malformed file: polygon $is:$iend not closed")
            end
        end

        is = iend+1
        i+=1
        if drop20
            is+=1
            i+=1
            drop20 = false
        end
    end
    return out
end

"""
    concat_poly(mpoly::Vector)

Concatenates a split poly, as returned from read_xyn, into one. Also
returns indices where to split apart again.  The concatenated polygon
is fully connected, with an edge going back to the first point.  This
allows to use `inpoly`, at least if the inner polygons have different
orientation to the outer. Also note that the input and output polygons
are closed, i.e. last point == first point.

Return:
- bigpoly -- with size==(2,n)
- splits -- ith poly has indices splits[i]:splits[i+1]-1
"""
function concat_poly(mpoly::Vector)
    T = eltype(mpoly[1])
    # total size is sum of sizes plus one extra point for all but the
    # first poly.
    totsize = mapreduce(x->size(x,2), +, mpoly) + length(mpoly) -1
    bigpoly = Array(T,size(mpoly[1],1),totsize)
    splits = Int[]
    is = 1
    for i=1:length(mpoly)
        push!(splits, is)
        @assert mpoly[i][:,1]==mpoly[i][:,end] "All input polys need to be closed."
        bigpoly[:,is:is+size(mpoly[i],2)-1] = mpoly[i]
        if i==1
            is = is+size(mpoly[i],2)
        else
            # add extra point to connect back to mpoly[1][:,1]
            is = is+size(mpoly[i],2) + 1
            bigpoly[:,is-1] = mpoly[1][:,1]
        end
    end
    push!(splits, totsize+1)
    return bigpoly, splits
end

"Split up concatenated polygon."
function split_poly{T}(bigpoly::Matrix{T}, splits)
    out = Matrix{T}[]
    if length(splits)==0
        return out
    else
        push!(out, bigpoly[:,splits[1]:splits[1+1]-1])
    end
    for i=2:length(splits)-1
        # remove the extra point joining to the first poly
        push!(out, bigpoly[:,splits[i]:splits[i+1]-2])
    end
    out
end


import RasterIO, Proj4
"""
    read_rasterio(fn::AbstractString, T=Float32; NA=convert(T,NaN))

Read various raster formats via the RasterIO.jl package.  Put output
into a Gridded instance and the Proj4 projection string.

"""
function read_rasterio(fn::AbstractString, T=Float32; NA=convert(T,NaN))
    ra = RasterIO.openraster(fn)
    nr = ra.height
    nc = ra.width
    #proj = RasterIO.getprojection(ra.dataset)
    proj4 = try # some computers may not have gdalsrsinfo installed
        strip(readstring(`gdalsrsinfo -o proj4 $fn`), ['\n', ''', ' '])
    catch
        ""
    end
    va = convert(Matrix{T}, RasterIO.fetch(ra,1))
    # get the NoData value
    aa = Cint[0]
    nodata, success = RasterIO.getrasternodatavalue(RasterIO.getrasterband(ra.dataset,1))
    if success
        for i in eachindex(va)
            if va[i]==nodata
                va[i] = NA
            end
        end
    end
    gt = RasterIO.geotransform(ra)
    origin = gt[[1,4]]
    xll,yll = RasterIO.applygeotransform(gt, 0.0, Float64(ra.height))
    pixelsz = gt[[2,6]]
    dx = pixelsz[1]
    dy = -pixelsz[2]
    @assert gt[[3,5]]==[0,0] "Can only handle North-up images"
    #Gridded(VAWTools.AGR(va', nc, nr, xll, yll, dx, nodata)), proj
    Gridded(range(xll+dx/2, dx, nc),
            range(yll+dy/2, dy, nr),
            true, flipdim(va,2)), proj4
end


## Misc helpers
###############
# """
# Cuts a smaller array out of a bigger one.
# """
#function cut_array(ar, )

# Does not work.  Probably better call out to gdalsrsinfo, see above.
# """
# Convert Well-Know-Text to Proj4.  Only works if there is a ESRI or
# ESPG number embedded in the WKT.
# """
# function wkt2proj4(wkt)
#     auth = search(wkt, "AUTHORITY")
#     if length(auth)==0
#         error("No AUTHORITY entry found in well-known-text")
#     end
#     endpos = search(wkt[auth[end]+1:end], "]")[1]
#   @show  li = eval(parse(wkt[auth[end]+1:auth[end]+endpos]))
#     nr = parse(Int, li[2])
#     if li[1]=="EPSG"
#         return Proj4.epsg[nr]
#     elseif li[1]=="EPSG"
#         return Proj4.esri[nr]
#     else
#         error("AUTHORITY entry mal-formed: $li")
#     end
# end

"""
    transform_proj(xy or xyz, from, to)

Transform between projections.  Uses Proj4.jl

Input:
- xy or xyz vector, matrix or trajectory of input points (x,y) or (x,y,z)
- `from` and `to` are either strings for Proj4.Projections (more preformant).
"""
function transform_proj(xyz, from, to)
    from = Proj4.Projection(from)
    to = Proj4.Projection(to)
    transform_proj(xyz, from, to)
end
function transform_proj(xyz, from::Proj4.Projection, to::Proj4.Projection)
    Proj4.transform(from, to, xyz)
end
function transform_proj(tr::Traj, from, to)
    xy = transform_proj(hcat(tr.x, tr.y), from, to) # TODO: avoid temporary
    if haserror(tr)
        Traj(xy[:,1], xy[:,2], tr.v, tr.err, tr.splits)
    else
        Traj(xy[:,1], xy[:,2], tr.v, tr.splits)
    end
end


"""
Convert int to string and pad with 0 to get to len
"""
int2str(i, len) = @sprintf "%05d" i

"""
    windnr(p, poly::Matrix)

Determines the winding number of a point and a polygon, i.e. how many
times a polygon winds around the point.

It follows Dan Sunday: http://geomalgorithms.com/a03-_inclusion.html.
"""
function windnr(p, poly::Matrix)
    @assert size(poly,1)==2
    @assert size(poly,2)>1
    @assert length(p)==2
    @assert poly[:,1]==poly[:,end]
    # Loop over edges
    px = p[1]
    py = p[2]
    wn = 0 # winding number, inside if odd
    len = length(poly)-2
    @inbounds for k=1:2:len #size(poly,2)-1 # @fastmath makes it worse
        # Coordinates of edge endpoints
        ex1,ey1,ex2,ey2 = poly[k], poly[k+1], poly[k+2], poly[k+3]

        # Check edges intersecting a horizontal ray:
        # rules http://geomalgorithms.com/a03-_inclusion.html#Edge-Crossing-Rules
        orient = leftorright(px,py,ex1,ey1,ex2,ey2)
        if up(ey1,ey2)
            # an upward edge includes its starting endpoint, and
            # excludes its final endpoint;
            !(py>=ey1 && py<ey2) && continue
            if orient==-1; wn-=1 end # p strictly left of oriented e
        elseif down(ey1,ey2)
            # a downward edge excludes its starting endpoint, and
            # includes its final endpoint;
            !(py<ey1 && py>=ey2) && continue
            if orient==1; wn+=1 end # p strictly right of oriented e
        else # Horizontal edges are excluded.
            continue
        end
        # NOTE: Rule "the edge-ray intersection point must be strictly
        # right of the point P."  implies that points on a right-side
        # boundary edge being outside, and ones on a left-side edge
        # being inside.
    end
    return wn
end
# edge goes up
up(ey1,ey2)   = ey1<ey2
# edge goes down
down(ey1,ey2) = ey1>ey2
function leftorright(px,py,ex1,ey1,ex2,ey2)
    # returns:
    # -1 if on left of line
    #  0 if on line
    #  1 if on right of line
    vx,vy = ex2-ex1, ey2-ey1
    px,py =  px-ex1,  py-ey1
    sign(px*vy-py*vx)
end

"""
Determines if a point is inside a polygon.

- p -- point (x,y) or [x,y]
- poly -- polygon vertices [x1 x2 ... xn x1
                            y1 y2 ... yn y1]
          (a closed poly)

Returns true if point has an odd winding number.  This should label
points as exterior which are inside outcrops.  See test for a test.

Note that points located on a right-side boundary edge are outside,
and ones on a left-side edge are inside.
"""
inpoly(p, poly::Matrix) = isodd(windnr(p,poly))

# Other inpoly algo here:
# https://github.com/helenchg/PolygonClipping.jl/blob/1fc74ab797c6585795283749b7c4cc9cb2000243/src/PolygonClipping.jl#L139

"""
Absolute value of slope angle (by finite differences)

- Averaged over 3x3 points.
- no slope is calculated for the outermost points

In:
- gridded elevation set
Out:
- slope angle (rad)

Note: slope and angle are very similar up to about ~0.4

TODO:
- allow only using points inside a mask
"""
function absslope(g::Gridded)
    nx, ny = size(g)
    dx,v = step(g.x), g.v
    alphas = zeros(Float64, nx, ny)
    for j=2:ny-1, i=2:nx-1
        dvx = ((v[i+1,j+1]+2*v[i+1,j]+v[i+1,j-1])-(v[i-1,j+1]+2*v[i-1,j]+v[i-1,j-1]))/(8*dx)
        dvy = ((v[i-1,j-1]+2*v[i,j-1]+v[i+1,j-1])-(v[i-1,j+1]+2*v[i,j+1]+v[i+1,j+1]))/(8*dx)
        alphas[i,j] = atan(sqrt(dvx^2+dvy^2))
    end
    return alphas
end


##############
# Bands
##############

"""
Bins a gird.  Often used to bin a DEM into elevation bands.

- g -- g.v of grid is binned
- binsize_or_bins -- bin size or the bins (a Range)
- mask -- specify if not all locations of a gird should be binned.]

KW:
- binround -- floor the bin-start using this many digits (see help of floor)

Return:
- bands -- a range of the bands, e.g. 0.0:10.0:100.0
- inds -- a Vector{Vector{Int}} of length(bands) with each element
          containing the indices of cells in the band
"""
function bin_grid(g::Gridded, binsize_or_bins, mask=BitArray([]); binround=-floor(Int, log10(binsize_or_bins)))
    if isempty(mask)
        v = g.v
        ginds = 1:length(v)
    else
        @assert size(mask)==size(g.v)
        v = g.v[mask]
        ginds = find(mask[:])
    end
    nv = length(v)
    if isa(binsize_or_bins, Number) # i.e. a binsize
        mi, ma = minimum(v), maximum(v)
        binstart = floor(mi, binround)
        binend = floor(ma, binround)
        bins = binstart:binsize_or_bins:binend # these are the start of the bins
    else
        bins = binsize_or_bins
    end
    _bin_grid_kernel(bins, nv, v, ginds)
end
@inbounds function _bin_grid_kernel(bins, nv, v, ginds)
    indices = Vector{Int}[]
    for b in bins
        ind = Int[]
        push!(indices, ind)
    end
    for j=1:nv
        if (bins[2]-bins[1])>0
            i = searchsortedfirst(bins, v[j])-1
        else
            # https://github.com/JuliaLang/julia/issues/18653
            i = searchsortedfirst(collect(bins), v[j], rev=true)-1
        end
        i = max(i,1) # TODO, this is a bit of a hack... otherwise if bins=0:10 and v[j]==0 -> i==0
        push!(indices[i], ginds[j])
    end
    return bins, indices
end


"""
Bins a trajectory using a grid

- tr -- trajectory to be binned
- g -- g.v of grid used for binning
- binsize_or_bins -- bin size or the bins (a Range)
- mask -- specify if not all locations of a gird should be binned.]

KW:
- binround -- floor the bin-start using this many digits (see help of floor)
"""
function bin_traj(tr::Traj, g::Gridded, binsize_or_bins, mask=trues(size(g.v)); binround=-floor(Int, log10(binsize_or_bins)))
    @assert size(mask)==size(g.v)

    demi  = interpolate((g.x, g.y), g.v, Interpolations.Gridded(Linear()) )
    maski  = interpolate((g.x, g.y), mask, Interpolations.Gridded(Constant()) ) # returns Int!

    v = [demi[x,y] for (x,y) in zip(tr.x,tr.y)]
    vm = Bool[maski[x,y] for (x,y) in zip(tr.x,tr.y)]
    v = v[vm]
    ginds = find(vm)
    nv = length(v)

    if isa(binsize_or_bins, Number)
        mi, ma = minimum(v), maximum(v)
        binstart = floor(mi, binround)
        binend = floor(ma, binround)
        bins = binstart:binsize_or_bins:binend # these are the start of the bins
    else
        bins = binsize_or_bins
    end
    _bin_grid_kernel(bins, nv, v, ginds)
end

"""
   map_onto_bands(bandi, field, fn=mean)

Map a field onto the elevation bands.  The field needs to have the
same size as the original binned-grid.

Input:
- bandi -- as returned by bin_grid
- field -- the field, either a Matrix or a Gridded
- fn -- the function to do the reduction with.  Default==mean

Output
- the value of the field on the bands
"""
function map_onto_bands(bandi, field::Matrix, fn=mean, fill=NaN)
    resT = typeof(fn(field[bandi[1]])) # to get result type
    out = zeros(resT, length(bandi))
    for i in 1:length(bandi)
        # drop any with fill values
        if any(isequal(field[bandi[i]],fill))
            out[i] = fill
        else
            out[i] = fn(field[bandi[i]])
        end
    end
    return out
end
map_onto_bands(bandi, field::Gridded, fn=mean) = map_onto_bands(bandi, field.v, fn=mean)


#################

"""
Return the index range of points within the distance dist.
"""
function around(g::Gridded, i, j, dist)
    inds = floor(Int,dist/step(g))
    nx, ny = size(g)
    return max(1,i-inds):min(nx,i+inds), max(1,j-inds):min(ny,j+inds)
end

"""
Return the linear index of points within the distance dist provided
they are within the mask.
"""
function around(g::Gridded, mask, i, j, dist)
    is, js = around(g::Gridded, i, j, dist)
    dims = size(mask)
    l = Int[] # linear index
    for j=js, i=is
        if mask[i,j]
            push!(l, ind2sub(mask, i, j))
        end
    end
    return l
end

"Mean ignoring NaNs"
function meannan(a)
    n = 0
    cum = zero(eltype(a))
    @inbounds for i in eachindex(a)
        if !isnan(a[i])
            n+=1
            cum+=a[i]
        end
    end
    cum/n # == NaN if n==0
end

"Std ignoring NaNs"
function stdnan(a)
    n = 0
    inds = Int[]
    @inbounds for i in eachindex(a)
        if !isnan(a[i])
            push!(inds,i)
        end
    end
    std(a[inds])
end


"Mean ignoring fill"
function meanfill(a,fill)
    n = 0
    cum = zero(eltype(a))
    @inbounds for i in eachindex(a)
        if a[i]!=fill
            n+=1
            cum+=a[i]
        end
    end
    if n==0
        println("Could not fill a gap! Returning fill.")
        return fill
    end
    cum/n # == NaN if n==0
end


"weighted mean"
function mean_weighted(a,w)
    # make an accumulator type closed under addition (needed for Bools):
    n = zero(one(eltype(w)) + one(eltype(w)))
    cum = zero(eltype(a))
    @inbounds for i in eachindex(a)
        n+=w[i]
        cum+=a[i]*w[i]
    end
    if n==0
        error("zero weights")
    end
    cum/n # == NaN if n==0
end

"""
    boxcar(A::AbstractArray, window[, weights])

Boxcar filter.  The two argument call ignores NaNs.  The three
argument call uses weights instead of NaNs, it can be a lot faster.

For the weights it may be faster to use non Bool arrays or BitArrays,
say Int8.

Also works for Vectors.

Smoothing occurs over x+/-window.

From http://julialang.org/blog/2016/02/iteration
"""
function boxcar(A::AbstractArray, window)
    out = similar(A)
    R = CartesianRange(size(A))
    I1, Iend = first(R), last(R)
    for I in R # @inbounds does not help
        if !isnan(A[I])
            out[I] = NaN
        end
        n, s = 0, zero(eltype(out))
        for J in CartesianRange(max(I1, I-I1*window), min(Iend, I+I1*window))
            if !isnan(A[J])
                s += A[J]
                n += 1
            end
        end
        out[I] = s/n
    end
    out
end

# this drops points which themselves have zero weight.
function boxcar(A::AbstractArray, window, weights::AbstractMatrix)
    @assert size(weights)==size(A)
    out = zeros(A)
    # make an accumulator type closed under addition (needed for Bools):
    T = typeof(one(eltype(weights)) + one(eltype(weights)))
    R = CartesianRange(size(A))
    I1, Iend = first(R), last(R)
    @inbounds @fastmath for I in R
        if weights[I]==0
            out[I] = 0
            continue
        end
        n, s = zero(T), zero(eltype(out))
        for J in CartesianRange(max(I1, I-I1*window), min(Iend, I+I1*window))
            s += A[J]*weights[J]
            n += weights[J]
        end
        out[I] = s/n
    end
    out
end

"""
    boxcar_matrix(T::DataType, window::Integer, weights::AbstractMatrix)

This produces a sparse matrix which can be used to apply the filter:
`bx*hs2d`. Expensive to do create but very fast to apply, thus use
when needing the same filter several times.
"""
function boxcar_matrix{T}(::Type{T}, window::Integer, weights::AbstractMatrix)
    # make an accumulator type closed under addition & division
    # (needed for Ints and Bools):
    Tacc = promote_type(T, eltype(weights))
    nr = size(weights,1)
    nc = size(weights,2)
    is = Int[]
    js = Int[]
    vs = T[]
    sizehint!(is, 2*window*nr*nc)
    sizehint!(js, 2*window*nr*nc)
    sizehint!(vs, 2*window*nr*nc)
    R = CartesianRange(size(weights))
    I1, Iend = first(R), last(R)
    @inbounds @fastmath for I in R
        if weights[I]==0
            # do nothing
            continue
        end
        i = (I.I[2]-1)*nr + I.I[1] # row of output matrix
        nrows = 0 # number of contributing cells
        acc = zero(Tacc) # sum of all weights for one cell
        for J in CartesianRange(max(I1, I-I1*window), min(Iend, I+I1*window))
            if weights[J]==0
                continue
            end
            j = (J.I[2]-1)*nr + J.I[1]
            push!(is, i)
            push!(js, j)
            push!(vs, weights[J])
            nrows +=1
            acc += weights[J]
        end
        # divide the current batch by the summed weights:
        lvs = length(vs)
        for n=lvs-(0:nrows-1)
            vs[n] /= acc
        end
    end
    return sparse(is, js, vs, length(weights), length(weights))
end

"""
Returns a piecewise linear function through points (xs,ys).
For x<xs[1]->ys[1] and for x>xs[end] -> ys[end]

When length(xs)==1, then it returns a constant function.

Currently extrapolates using the last value. (TODO)

Almost 2x faster than using Interpolations.jl.
"""
function piecewiselinear(xs::AbstractVector, ys::AbstractVector)
    if length(xs)==1
        @assert length(xs)==length(ys)
        return xq -> ys[1]
    end
    rats = diff(ys)./diff(xs)
    if issorted(xs)
        return function (xq)
            # xq<xs[1] && error("cannot extrapolate")
            # xq>xs[end] && error("cannot extrapolate")
            # xq==xs[end] && return ys[end]
            xq<=xs[1] && return ys[1]
            xq>=xs[end] && return ys[end]
            ii = searchsortedlast(xs, xq)
            out = ys[ii] - (xs[ii]-xq)*rats[ii]
            return out
        end
    elseif issorted(xs, Base.Order.Reverse)
        return function (xq)
            # xq>xs[1] && error("cannot extrapolate")
            # xq<xs[end] && error("cannot extrapolate")
            # xq==xs[end] && return ys[end]
            xq>=xs[1] && return ys[1]
            xq<=xs[end] && return ys[end]
            ii = searchsortedlast(xs, xq, Base.Order.Reverse)
            out = ys[ii] - (xs[ii]-xq)*rats[ii]
            return out
        end
    else
        error("xs must be sorted")
    end
end
# """
# Includes extrapolation:
# Returns a piecewise linear function through points (xs,ys).
# For x<xs[1]->ys[1] and for x>xs[end] -> ys[end]

# - assumes xs is sorted

# Almost 2x faster than using Interpolations.jl.
# """
# function piecewiselinear(xs_, ys_)
#     rats = diff(ys_)./diff(xs_)
#     ys = copy(ys_)
#     prepend!(ys,ys[1:1])
#     append!(ys,ys[end:end])

#     xs = copy(xs_)
#     prepend!(xs, xs[1:1])
#     append!(xs, xs[end:end])

#     prepend!(rats, [0.0])
#     append!(rats,  [0.0])

#     function (xq)
#         ii = searchsortedfirst(xs_, xq)
#         return ys[ii] - (xs[ii]-xq)*rats[ii]
#     end
# end

###########
# Smoothing splines
###########
import StatsBase, SmoothingSplines
const SSp = SmoothingSplines

"""

Smoothes a vector of datapoints with a spline with minimum wave-length
components `len`.  Returns a smoothed vector of the same length.

Notes:
- if len=0 no smoothing occurs
- better than boxcar if having no edge-effects is important
"""
function smooth_vector{T}(x, y::AbstractVector{T}, len, out=y)
    if len==0 && out==y
        return y
    else
        lambda0 = 0.001 # by trial and error
        lambda = convert(T, lambda0 * len^3) # the ^3 seems correct (trial and error)
        if out==y
            return SSp.predict(StatsBase.fit(SSp.SmoothingSpline, convert(Vector{T}, x), y, lambda))
        else
            return SSp.predict(StatsBase.fit(SSp.SmoothingSpline, convert(Vector{T}, x), y, lambda), out)
        end
    end
end


# using Dierckx
## using the number of nodes does funny stuff to the curve...
# """
# Smoothes a vector using a spline with nknots.  Returns a smoothed
# vector of the same length.

# Notes:
# - if nknots+2>=length(x) then it just returns y
# - better than boxcar if having no edge-effects is important
# """
# function smooth(x, y, nknots)
#     if nknots+2>=length(x)
#         return y
#     else
#         if !issorted(x)
#             perm = sortperm(x)
#             x_ = x[perm]
#             y = y[perm]
#         else
#             x_ = x
#         end
#         xknots = linspace(x_[1],x_[end],nknots+2)[2:end-1]
#         return Spline1D(x_, y, xknots)(x)
#     end
# end

end # module
