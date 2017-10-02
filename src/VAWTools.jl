__precompile__() # RasterIO is not pre-compiled
warn("BREAKING CHANGE: The Gridded, Gridded1d and Traj interface was overhauled.")

module VAWTools
using Compat
using Parameters

include("other-pks-fixes.jl")

# General tools, maybe applicable more widely.
export read_agr, write_agr, read_xyn, inpoly, AGridded, Gridded, Gridded1d, Traj,
    smooth_vector, absslope, gradient3by3,
    downsample, split_traj!, boxcar, boxcar_matrix, apply_boxcar_matrix, bin_grid, piecewiselinear, split_poly,
    transform_proj


import Base: ==, size, length, step, +, -, *

# const spatialorder = "yx"  # order of imagestorage, https://github.com/timholy/Images.jl#storage-order-and-changing-the-representation-of-images

@compat abstract type AGridded{T} end
"Check whether an error-field is defined."
haserror(g::AGridded) = size(g.err)!=(0,0)
hasproj(g::AGridded) = g.proj!=""
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
- hold NA value?
- default of err, see https://github.com/JuliaLang/julia/issues/23926
"""
@with_kw struct Gridded{T} <: AGridded{T}
    x::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
    y::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
    v::Matrix{T}  # values
    @assert size(v)==(length(x), length(y))
    err::Matrix{T}=Matrix{eltype(v)}(0,0)  # error of values
    midpoint::Bool=true  # if true, the (x,y) is cell midpoint, otherwise lower-left corner
    proj::String="" # Proj4 string
end
Gridded(x,y,v::Matrix{T}) where T = Gridded{T}(x,y,v,Matrix{T}(0,0),true,"")
Gridded(x,y,v::Matrix{T},err) where T = Gridded{T}(x,y,v,err,true,"")

"""
    downsample(g::Gridded, step, start=1, average=true, averagemask=all_points)

Makes a `Gridded` smaller by downsampling.  If averaging (by default),
the window is chosen such that there is no overlap (for odd `step`) or
one cell overlap (for even `step`).

"""
function downsample(g::Gridded, step, start=1, average=true, averagemask=UniformArray{Bool,2}(true))
    rx = start:step:size(g.v,1)
    ry = start:step:size(g.v,2)
    nx,ny = length(rx),length(ry)
    @assert nx>1
    @assert ny>1
    x = g.x[rx]
    y = g.y[ry]
    midpoint = g.midpoint
    # Window: no overlap on odd step sizes
    window = average ? Int(floor((step)/2)) : 1
    vnew = _downsample_boxcar(g.v,rx,ry,window,averagemask)

    errnew = haserror(g) ? _downsample_boxcar(g.err,rx,ry,window,averagemask) : g.err
    Gridded(x,y,vnew,errnew,g.midpoint,g.proj)
end
function _downsample_boxcar(v,rx,ry,window,averagemask)
    # make type which allows for averaging
    T = typeof(one(eltype(v))/2)
    vnew = zeros(T,length(rx), length(ry))
    if window==0
        for (j,rj) in enumerate(ry)
            for (i,ri) in enumerate(rx)
                vnew[i,j]=v[ri,rj]
            end
        end
        return vnew
    end
    # Follows boxcar filter below
    for (j,rj) in enumerate(ry)
        for (i,ri) in enumerate(rx)
            # keep NaNs or points outside averagemask as they are:
            if isnan(v[ri,rj]) || !averagemask[ri,rj]
                vnew[i,j] = v[ri,rj]
                continue
            end
            # other points make an average:
            v0 = zero(T)
            n = 0
            for jj=max(1,rj-window):min(size(v,2),rj+window)
                for ii=max(1,ri-window):min(size(v,1),ri+window)
                    tmp = v[ii,jj]
                    if !isnan(tmp) && averagemask[ii,jj]
                        v0 += tmp
                        n+=1
                    end
                end
            end
            vnew[i,j] = v0/n
        end
    end
    return vnew
end

function (+)(g1::Gridded, g2::Gridded)
    @assert g1.midpoint==g2.midpoint
    @assert g1.proj==g2.proj
    Gridded(g1.x, g1.y, g1.v.+g2.v, g1.err.+g2.err, g1.midpoint, g1.proj)
end
(*)(n::Number, g::Gridded) = Gridded(g.x, g.y, g.v*n, g.err*n, g.midpoint, g.proj)
(*)(g::Gridded, n::Number) = n*g


"""
Holds 1D fields, e.g. elevation band data.

Be sure to be clear whether x corresponds to cell centers or boundaries.
"""
@with_kw struct Gridded1d{T} <:  AGridded{T}
    x::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
    v::Vector{T} # values
    @assert size(v)==(length(x), )
    err::Vector{T}=Vector{eltype(v)}(0) # error of values
    @assert size(v)==size(err) || size(err)==(0,)
    midpoint::Bool=true  # if true, the (x) is cell midpoint, otherwise lower corner
    proj::String="" # Proj4 string
end
Gridded1d(x,v::Matrix{T}) where T = Gridded1d{T}(x,v,T[],true,"")
Gridded1d(x,v::Matrix{T},err) where T = Gridded1d{T}(x,v,err,true,"")


# function (+)(g1::Gridded1d, g2::Gridded1d)
#     @assert g1.midpoint==g2.midpoint
#     @assert g1.proj==g2.proj
#     Gridded1d(g1.x, g1.v+g2.v, g1.err+g2.err, g1.midpoint, g1.proj)
# end
(*)(n::Number, g::Gridded1d) = Gridded1d(g.x, g.v*n, g.err*n, g.midpoint, g.proj)
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
@with_kw mutable struct Traj{T}
    x::Vector{Float64}
    y::Vector{Float64}
    @assert length(x)==length(y)
    v::Vector{T}=Void[]
    @assert length(v)==length(y) || length(v)==0
    err::Vector{T}=Vector{eltype(v)}(0)
    @assert length(err)==length(y) || length(err)==0
    splits::Vector{UnitRange{Int}}=UnitRange{Int}[1:length(x)]
    proj::String="" # Proj4 string
end
Traj(x,y) = Traj{Void}(x,y,Void[],Void[],UnitRange{Int}[1:length(x)],"")
Traj(x,y,v::AV) where AV<:AbstractVector{T} where T = Traj{T}(x,y,v,T[],UnitRange{Int}[1:length(x)],"")
Traj(x,y,v::AV,err) where AV<:AbstractVector{T} where T = Traj{T}(x,y,v,err,UnitRange{Int}[1:length(x)],"")
Traj(x,y,v::AV,err,splits) where AV<:AbstractVector{T} where T = Traj{T}(x,y,v,err,splits,"")

hasvalues(g::Traj) = length(g.v)!=0
haserror(g::Traj) = length(g.err)!=0
hasproj(g::Traj) = g.proj!=""
length(g::Traj) = length(g.x)

"""
    split_traj!{T<:Traj}(t::T, dist)

Splits trajectory into several
"""
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
    return Traj(x, y, convert(Vector{T2}, g.v), convert(Vector{T2}, g.err), g.splits, proj)
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
    hasutm::Bool    # Matthias sometimes abuses the NODATA_value field as UTM-zone field
    extra_header::Vector{Float32} # the .bin files have space for
                                  # extra information in the header
    function (::Type{AGR{T}}){T}(va, nc, nr, xll, yll, dx, NA,
                                 hasutm=false, extra_header=zeros(Float32,6))
        @assert size(va)==(nr,nc)
        @assert dx>=0
        new{T}(va, nc, nr, xll, yll, dx, NA, hasutm, extra_header)
    end
end
AGR{T}(va::Matrix{T}, nc, nr, xll, yll, dx, NA, hasutm=false, extra_header=zeros(Float32,6)) = AGR{T}(va, nc, nr, xll, yll, dx, NA, hasutm, extra_header)

"""
Transform Gridded to AGR

Optional:
- NA_g: no-data value of g:Gridded, defaults to NaN
- NA_agr: no-data value of AGR output, defaults to NaN
- utmzone: if 0<utmzone<61, then above two are ignored and instead a file with UTM_ZONE is written.
- write_err: write g.err instead of g.v

"""
function AGR{T}(g::Gridded{T}; NA_g=convert(T,NaN), NA_agr=convert(T,NaN), write_err=false, utmzone=0)
    v = write_err ? g.err : g.v # whether to write the values or errors into the ascii grid
    if 0<utmzone<61
        hasutm = true
    else
        hasutm = false
        if !isequal(NA_g,NA_agr)
            v = deepcopy(v)
            for i=eachindex(v)
                if isequal(v[i],NA_g)
                    v[i] =  NA_agr
                end
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
    if hasutm
        AGR(v, nc, nr, xll, yll, dx, utmzone, hasutm)
    else
        AGR(v, nc, nr, xll, yll, dx, NA_agr)
    end
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
    proj = agr.hasutm ? "+proj=utm +zone=$(agr.NA)" : ""
    Gridded(range(agr.xll+agr.dx/2, agr.dx, agr.nc),
            range(agr.yll+agr.dx/2, agr.dx, agr.nr),
            v,
            zeros(T,0,0),
            true,
            proj)
end

# File readers
##############

"""
Check if .bin file nor not.  With IO types this may not work 100%.
"""
isbin_file(fn::AbstractString) = splitext(fn)[2]==".bin"
function isbin_file(fn::IO)
    if :name in fieldnames(typeof(fn))
        f = fn.name
    else # this should catch TranscodingStreams.jl
        f = fn.stream.name
    end
    f = strip(f, ['<','>']) # often looks like "<file asdf.ext>"
    f = strip(f)
    isbin_file(f)
end

"""
    read_agr(fl, T=Float32; NA=convert(T,NaN))

Read a Ascii grid https://en.wikipedia.org/wiki/Esri_grid.  Can also
read from .gz compressed files, if `import CodecZlib`.

In:
- fn -- file name

Optional keywords:
- T -- type of output array. Defaults to Float32.
- NA -- replace the fill value with this. Defaults to use NaN.

Out: Gridded
"""
function read_agr(fl, T=Float32; NA=convert(T,NaN))
    Gridded(_read_agr(fl, T, NA), NA=NA)::Gridded{T}
end

function _read_agr(fl::AbstractString, T=Float32, NA=nothing)
    if !isfile(fl)
        error("File $fl does not exist!")
    end
    if endswith(fl, ".gz")
        @eval import CodecZlib
        io = CodecZlib.GzipDecompressionStream(open(fl))
    else
        io = open(fl, "r")
    end
    out = _read_agr(io, T, NA)
    close(io)
    return out
end

function _read_agr(io::IO, T=Float32, NA=nothing)
    if isbin_file(io)
        toT = (io,T) -> convert(T, read(io, Float32))
    else
        toT = (io,T) -> parse(T, split(readline(io))[2])
    end

    local va::Matrix{T}
    # read header
    nc = toT(io, Int)
    nr = toT(io, Int)
    xll = toT(io,T)
    yll = toT(io,T)
    dx = toT(io,T)
    # Matthias sometimes abuses the NODATA_value field as UTM-zone field
    # in non-binary grids:
    if isbin_file(io)
        fill = toT(io,T)
        hasutm = false
    else
        prop, val = split(readline(io))
        if prop=="NODATA_value"
            hasutm = false
        elseif prop=="UTM_ZONE"
            hasutm = true
        else
            error("Unrecognized last header field: $prop")
        end
        fill = parse(T,val)
    end

    if isbin_file(io)
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
        va = Array{T}(nr, nc)
        # no extra header for ascii .agr
        extra_header = zeros(Float32, 6)
        # read values
        tmp = split(readstring(io))
        if length(tmp)!=nc*nr
            error("Something's wrong with the file/stream $io.  It should contain $(nc*nr) values but has $(length(tmp))")
        end
        for i=1:nr, j=1:nc
            va[i,j] = parse(T, tmp[(i-1)*nc + j])
        end
    end
    # replace missing value with something else
    if NA!=nothing && fill!=NA && hasutm
        fill = _refill!(va, fill, NA)
    end
    # make a AGR data structure
    AGR(va, nc, nr, xll, yll, dx, fill, hasutm, extra_header)
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
        if !isbin_file(fn)
            if g.hasutm
                write(io, "UTM_ZONE      ")
            else
                write(io, "NODATA_value  ")
            end
        end
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
        if g.hasutm
            write(io, convInt(fill))
        else
            write(io, convT(fill))
        end
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
    bigpoly = Array{T}(size(mpoly[1],1),totsize)
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

"Finds splitting points in a big-poly"
function find_poly_splits(bigpoly::Matrix)
    splits = Int[]
    p1 = bigpoly[:,1]
    for i=1:size(bigpoly,2)
        if p1==bigpoly[:,i]
            push!(splits,i)
        end
    end
    splits
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

import Proj4
# # Make RasterIO conditional as it is not Julia 0.6 compatible
# if haskey(Pkg.installed(), "RasterIO")
#     println("Enabling RasterIO function read_rasterio")
#     eval(quote
#         import RasterIO
#         """
#         read_rasterio(fn::AbstractString, T=Float32; NA=convert(T,NaN))

#         Read various raster formats via the RasterIO.jl package.  Put output
#         into a Gridded instance and the Proj4 projection string.
#         """
#         function read_rasterio(fn::AbstractString, T=Float32; NA=convert(T,NaN))
#             ra = RasterIO.openraster(fn)
#             nr = ra.height
#             nc = ra.width
#             #proj = RasterIO.getprojection(ra.dataset)
#             proj4 = try # some computers may not have gdalsrsinfo installed
#                 strip(readstring(`gdalsrsinfo -o proj4 $fn`), ['\n', '\'', ' '])
#             catch
#                 ""
#             end
#             va = convert(Matrix{T}, RasterIO.fetch(ra,1))
#             # get the NoData value
#             aa = Cint[0]
#             nodata, success = getrasternodatavalue(RasterIO.getrasterband(ra.dataset,1))
#             if success
#                 for i in eachindex(va)
#                     if va[i]==nodata
#                         va[i] = NA
#                     end
#                 end
#             end
#             gt = RasterIO.geotransform(ra)
#             origin = gt[[1,4]]
#             xll,yll = RasterIO.applygeotransform(gt, 0.0, Float64(ra.height))
#             pixelsz = gt[[2,6]]
#             dx = pixelsz[1]
#             dy = -pixelsz[2]
#             @assert gt[[3,5]]==[0,0] "Can only handle North-up images"
#             #Gridded(VAWTools.AGR(va', nc, nr, xll, yll, dx, nodata)), proj
#             Gridded(range(xll+dx/2, dx, nc),
#                     range(yll+dy/2, dy, nr),
#                     true, flipdim(va,2)), proj4
#         end
#     end)
# end

## Shapefiles
import Shapefile, GeoInterface

"""
    read_shapefile_poly(filename)

Reads a shapefile and returns a list of the polygons it contains.
"""
function read_shapefile_poly(fln)
    shp = open(fln) do fd
        read(fd, Shapefile.Handle)
    end
    outt = []
    for ii =1:length(shp.shapes)
        poly = shp.shapes[1]
        np = length(poly.points)
        pp = poly.parts
        push!(pp, np)

        out = zeros(2,np+length(pp)-2)
        extra_j = 0
        for i = 2:length(pp)
            for j = (pp[i-1]+1:pp[i])+extra_j
                out[:,j] = GeoInterface.coordinates(poly.points[j-extra_j])
            end
            if i!=length(pp)
                out[:,pp[i]+extra_j+1] = out[:,1]
                extra_j += 1
            end
        end
        push!(outt, out)
    end
    return outt
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

Transform between projections.  Uses Proj4.jl.  Example:

    longlat = "+proj=longlat"
    utm56 = "+proj=utm +zone=56 +south"
    transform_proj([6e3,197e3], utm56, longlat)

Input:
- xy or xyz vector, matrix or trajectory of input points (x,y) or (x,y,z)
- `from` and `to` are either strings for Proj4.Projections (more preformant).

See VAWTools.Projections for pre-defined projections
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

module Projections
using Proj4
"Swiss grid 1903"
const swiss1903 = Proj4.Projection("+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +k_0=1 +x_0=600000 +y_0=200000 +ellps=bessel +towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +no_defs")
"Longitude latitude"
const longlat = Proj4.Projection("+proj=longlat")

# UTM zones: utm1n, utm1s, etc
# https://github.com/JuliaGeo/Proj4.jl
for z=1:60
    zs = Symbol("utm$(z)s")
    zn = Symbol("utm$(z)n")
    strs = "+proj=utm +zone=$z +south +datum=WGS84 +units=m +no_defs"
    strn = "+proj=utm +zone=$z +north +datum=WGS84 +units=m +no_defs"
    @eval $zs = Projection($strs)
    @eval $zn = Projection($strn)
end
end

######
# Distances
######
"""
    dist(x1,y1,x2,y2)

Cartesian distance
"""
dist(x1,y1,x2,y2) = sqrt((x1-x2)^2 + (y1-y2)^2)

"""
    dist_longlat(lon1,lat1,lon2,lat2,R=6373)

Distance on a sphere (default earth size).
"""
function dist_longlat(lon1,lat1,lon2,lat2,R=6373) # the radius of the Earth
    # http://andrew.hedges.name/experiments/haversine/

    dlon = deg2rad(lon2 - lon1)
    dlat = deg2rad(lat2 - lat1)
    lat1 = deg2rad(lat1)
    lat2 = deg2rad(lat2)
    a = (sin(dlat/2))^2 + cos(lat1) * cos(lat2) * (sin(dlon/2))^2
    c = 2 * atan2( sqrt(a), sqrt(1-a) )
    return R * c
end



########
"""
Convert int to string and pad with 0 to get to len
"""
int2str5(i) = @sprintf "%05d" i

"""
Represents an array all filled with one value.

Pretty incomplete implementation.  Probably needs something similar as
https://github.com/JuliaLang/julia/blob/86bf95fe0a76e4750d41f569bfcb6fa1fb1805e4/base/linalg/uniformscaling.jl.

TODO:
- arithmetic
"""
immutable UniformArray{T,N} <: AbstractArray{T,N}
    val::T
end
Base.size{T,N}(::UniformArray{T,N}) = ntuple(x->typemax(Int), Val{N})
Base.getindex(A::UniformArray, i::Int) = A.val
@compat Base.IndexStyle{U<:UniformArray}(::Type{U}) = IndexLinear()
Base.start(::UniformArray) = error("Cannot iterate over UniformArray")
# Does not work https://github.com/JuliaLang/julia/issues/18004
#Base.show{T,N}(io::IO, u::UniformArray{T,N}) = print(io, "UniformArray{$T,$N} with value $(u.val)")
Base.show{T,N}(io::IO, ::MIME"text/plain", u::UniformArray{T,N}) = print(io, "UniformArray{$T,$N} with value $(u.val)")

#############
# Polygons
#############
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
    absslope(g::Gridded,weights::AbstractMatrix=UniformArray{Bool,2}(true),retnan=false)
    absslope(x::Range,y::Range,v,weights::AbstractMatrix=UniformArray{Bool,2}(true),retnan=true)

Absolute value of slope angle (by finite differences)

- uses a 3x3 points stencil.
- on outermost points uses one-sided stencil

In:
- gridded elevation set
Out:
- slope angle (rad)

Note: slope and angle are very similar up to about ~0.4 or 20deg
"""
absslope(g::Gridded,weights::AbstractMatrix=UniformArray{Bool,2}(true),retnan=true) =
    absslope(g.x,g.y,g.v,weights,retnan)
function absslope(x::Range,y::Range,v,weights::AbstractMatrix=UniformArray{Bool,2}(true),retnan=true)
    dvx,dvy = gradient3by3(x,y,v,weights,retnan)
    for i in eachindex(dvx)
        dvx[i] = atan(sqrt(dvx[i]^2+dvy[i]^2))
    end
    dvx
end

## Old, probably faster implementation:
# absslope(g::Gridded) = absslope(g.x,g.y,g.v)
# function absslope(x::Range,y::Range,v)
#     nx, ny = length(x),length(y)
#     dx = step(x)
#     alphas = zeros(Float64, nx, ny)
#     for j=2:ny-1, i=2:nx-1
#         dvx =  ((v[i+1,j+1]+2*v[i+1,j]+v[i+1,j-1])-(v[i-1,j+1]+2*v[i-1,j]+v[i-1,j-1]))/(8*dx)
#         dvy = -((v[i-1,j-1]+2*v[i,j-1]+v[i+1,j-1])-(v[i-1,j+1]+2*v[i,j+1]+v[i+1,j+1]))/(8*dx)
#         alphas[i,j] = atan(sqrt(dvx^2+dvy^2))
#     end
#     # on edge and corners use one-sided f-d
#     for j=[1,ny], i=1:nx
#         if i==1
#             dvx = (v[i+1,j]-v[i,j])/dx
#         elseif i==nx
#             dvx = (v[i,j]-v[i-1,j])/dx
#         else
#             dvx = (v[i+1,j]-v[i-1,j])/(2*dx)
#         end
#         if j==1
#             dvy = (v[i,j+1]-v[i,j])/dx
#         else
#             dvy = (v[i,j]-v[i,j-1])/dx
#         end
#         alphas[i,j] = atan(sqrt(dvx^2+dvy^2))
#     end
#     for j=1:ny, i=[1,nx]
#         if i==1
#             dvx = (v[i+1,j]-v[i,j])/dx
#         else
#             dvx = (v[i,j]-v[i-1,j])/dx
#         end
#         if j==1
#             dvy = (v[i,j+1]-v[i,j])/dx
#         elseif j==ny
#             dvy = (v[i,j]-v[i,j-1])/dx
#         else
#             dvy = (v[i,j+1]-v[i,j-1])/(2*dx)
#         end
#         alphas[i,j] = atan(sqrt(dvx^2+dvy^2))
#     end
#     return alphas
# end

"""
    gradient3by3(g::Gridded)
    gradient3by3(x::Range,y::Range,v)
    gradient3by3(x::Range,y::Range,v,weights::AbstractMatrix,retnan=true)

2D gradient (by finite differences)

- Averaged over 3x3 points.
- no slope is calculated for the outermost points
- with `weights` is given then it can be used for masked areas. If
  `retnan==false` a 0 is returned where no gradient can be made.

In:
- gridded elevation set or x,y,z
Out:
- dvx,dvy

Reference: gradient3by3_fast
"""
gradient3by3(g::Gridded,weights::AbstractMatrix=UniformArray{Bool,2}(true),retnan=true) =
    gradient3by3(g.x,g.y,g.v,weights,retnan)
# Implementation using weights: drops pairs which don't work.
function gradient3by3(x::Range,y::Range,v,weights::AbstractMatrix=UniformArray{Bool,2}(true),retnan=true)
    nx, ny = length(x),length(y)
    dx = step(x)
    dvx = zeros(Float64, nx, ny)
    dvy = zeros(Float64, nx, ny)

    # local weights
    lw = zeros(Float64,3,3)
    f = [1,2,1]
    for j=1:ny, i=1:nx
        lw[:] = 0
        for jj=max(1,j-1):min(ny,j+1), ii=max(1,i-1):min(nx,i+1)
            lw[ii-i+2,jj-j+2] = weights[ii,jj]
        end

        # xdir
        n = 0
        for k=-1:1
            if lw[3,k+2]!=0 && lw[1,k+2]!=0
                # centered difference
                dvx[i,j] += (v[i+1,j+k]-v[i-1,j+k])/(2*dx)*f[k+2]
                n+=1*f[k+2]
            elseif lw[2,k+2]!=0 && lw[1,k+2]!=0
                # right difference
                dvx[i,j] += (v[i,j+k]-v[i-1,j+k])/dx*f[k+2]
                n+=1*f[k+2]
            elseif lw[3,k+2]!=0 && lw[2,k+2]!=0
                # left difference
                dvx[i,j] += (v[i+1,j+k]-v[i,j+k])/dx*f[k+2]
                n+=1*f[k+2]
            end
        end
        if n==0;
            if retnan
                dvx[i,j] = NaN
            else
                dvx[i,j] = 0
            end
        else
            dvx[i,j] /= n
        end
        # ydir
        n = 0
        for k=-1:1
            if lw[k+2,3]!=0 && lw[k+2,1]!=0
                # centered difference
                dvy[i,j] += (v[i+k,j+1]-v[i+k,j-1])/(2*dx)*f[k+2]
                n+=1*f[k+2]
            elseif lw[k+2,3]!=0 && lw[k+2,2]!=0
                # right difference
                dvy[i,j] += (v[i+k,j+1]-v[i+k,j])/dx*f[k+2]
                n+=1*f[k+2]
            elseif lw[k+2,2]!=0 && lw[k+2,1]!=0
                # left difference
                dvy[i,j] += (v[i+k,j]-v[i+k,j-1])/dx*f[k+2]
                n+=1*f[k+2]
            end
        end
        if n==0;
            if retnan
                dvy[i,j] = NaN
            else
                dvy[i,j] = 0
            end
        else
            dvy[i,j] /= n
        end
    end
    return dvx, dvy
end

"Faster than gradient3by3 but does not deal with weights nor the edge."
gradient3by3_fast(g::Gridded) = gradient3by3_fast(g.x,g.y,g.v)
function gradient3by3_fast(x::Range,y::Range,v)
    nx, ny = length(x),length(y)
    dx = step(x)
    dvx = zeros(Float64, nx, ny)
    dvy = zeros(Float64, nx, ny)
    for j=2:ny-1, i=2:nx-1
        dvx[i,j] =  ((v[i+1,j+1]+2*v[i+1,j]+v[i+1,j-1]) - (v[i-1,j+1]+2*v[i-1,j]+v[i-1,j-1]))/(8*dx)
        dvy[i,j] = -((v[i-1,j-1]+2*v[i,j-1]+v[i+1,j-1]) - (v[i-1,j+1]+2*v[i,j+1]+v[i+1,j+1]))/(8*dx)
    end
    return dvx, dvy
end


##############
# Bands
##############
include("elevation-bands.jl")

#################
# Filtering
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
function meanfill(a,fill,verbose=false)
    n = 0
    cum = zero(eltype(a))
    @inbounds for i in eachindex(a)
        if a[i]!=fill
            n+=1
            cum+=a[i]
        end
    end
    if n==0
        verbose && println("Could not fill a gap! Returning fill.")
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
    boxcar(A::AbstractArray, window[, weights, dropmask])

Boxcar filter.  The two argument call ignores NaNs.  The three & four
argument call uses weights instead of NaNs, it can be a lot faster.

For the weights it may be faster to use non-Bool arrays nor BitArrays,
say Int8.  Note that even NaNs where weight==0 will poison the result!

No average is calculated for points where dropmask==true, instead
their original value will be used.

Also works for Vectors.

Smoothing occurs over +/-window indices.

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
function boxcar(A::AbstractArray, window::AbstractArray)
    out = similar(A)
    R = CartesianRange(size(A))
    I1, Iend = first(R), last(R)
    for I in R # @inbounds does not help
        if !isnan(A[I])
            out[I] = NaN
        end
        n, s = 0, zero(eltype(out))
        for J in CartesianRange(max(I1, I-I1*window[I]), min(Iend, I+I1*window[I]))
            if !isnan(A[J])
                s += A[J]
                n += 1
            end
        end
        out[I] = s/n
    end
    out
end

# # this drops points which themselves have zero weight.
# #function boxcar{T1,T2,T3,N}(A::AbstractArray{T1,N}, window::AbstractArray{T2,N}, weights::AbstractArray{T3,N}=UniformArray{T,N}(1))
# function boxcar(A::AbstractArray, window, weights::AbstractArray)
#     @assert size(weights)==size(A)
#     out = zeros(A)
#     # make an accumulator type closed under addition (needed for Bools):
#     T = typeof(one(eltype(weights)) + one(eltype(weights)))
#     R = CartesianRange(size(A))
#     I1, Iend = first(R), last(R)
#     @inbounds @fastmath for I in R
#         if weights[I]==0
#             out[I] = 0
#             continue
#         end
#         n, s = zero(T), zero(eltype(out))
#         for J in CartesianRange(max(I1, I-I1*window), min(Iend, I+I1*window))
#             s += A[J]*weights[J]
#             n += weights[J]
#         end
#         out[I] = s/n
#     end
#     out
# end

# this does not drop points which themselves have zero weight, only points on dropmask.
function boxcar{T,N}(A::AbstractArray{T,N}, window,
                     weights::AbstractArray,
                     dropmask::AbstractArray=(weights.==0))
    @assert size(weights)==size(A)
    out = zeros(A)
    # make an accumulator type closed under addition (needed for Bools):
    AT = typeof(one(eltype(weights)) + one(eltype(weights)))
    R = CartesianRange(size(A))
    I1, Iend = first(R), last(R)
    @inbounds @fastmath for I in R
        if dropmask[I]
            out[I] = A[I]
            continue
        end
        n, s = zero(AT), zero(T)
        for J in CartesianRange(max(I1, I-I1*window), min(Iend, I+I1*window))
            s += A[J] * convert(T, weights[J])
            n += weights[J]
        end
        out[I] = s/n
    end
    out
end


"""
    boxcar_matrix(T::DataType, window::Integer, weights::AbstractMatrix[, dropmask])

This produces a sparse matrix which can be used to apply the filter:
`bx*hs2d`. Relatively expensive to create but very fast to apply, thus
use when needing the same filter several times.

Note, that here NaNs where weight==0 do not poison the result.

Apply with

    apply_boxcar_matrix(M, orig)
"""
function boxcar_matrix{T,TT,N}(::Type{T}, window::Integer,
                               weights::AbstractArray{TT,N},
                               dropmask::AbstractArray=(weights.==0))
                               # dropmask::AbstractArray=BitArray{N}(ntuple(x->0, Val{N})...))
    # make an accumulator type closed under addition (needed for Ints and Bools):
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
        i = (I.I[2]-1)*nr + I.I[1] # row of output matrix
        if dropmask[I]
            # do not average at this location, preserve original value
            push!(is, i)
            push!(js, i)
            push!(vs, 1)
            continue
        end
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

"Apply the boxcar filter matrix and reshape result"
apply_boxcar_matrix(M, orig) = reshape(M*reshape(orig, length(orig)), size(orig))

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
