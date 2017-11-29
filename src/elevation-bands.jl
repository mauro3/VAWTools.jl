##############
# Bands
##############
#
# TODO: add extrapolation too (second half of make_bands)

export map_onto_bands, make_1Dglacier, map_back_to_2D, map_back_to_2D!

_binround(binsize::Number) = -floor(Int, log10(abs(binsize)))
_binround(binsize) = 0

"""
    make_1Dglacier(dem::Gridded, binsize_or_bins, glaciermask=BitArray([]);
                        binround=_binround(binsize_or_bins),
                        window_dem_smooth=0.0,
                        window_width_smooth=0.0,
                        alpha_min=deg2rad(0.4),
                        alpha_max=deg2rad(60.0),
                        FILL=-9999999.0)

Makes a 1D glacier from a 2D DEM by using Huss' elevation band trick.

Returns:
- bands -- elevation bands.  The i-th band is (bands[i], bands[i]+step(bands))
- bandi -- linear index into dem which assigns each cell to a elevation band
- alphas, areas, lengths, widths, x, xmid -- elevation band slope, area, etc
- dem -- the used DEM, a possibly smoothed version of the input DEM
- alpha2d -- slopes at each point of `dem`
"""
function make_1Dglacier(dem::Gridded, binsize_or_bins, glaciermask=trues(size(dem.v));
                        binround=_binround(binsize_or_bins),
                        window_dem_smooth=0.0,
                        window_width_smooth=0.0,
                        alpha_min=deg2rad(0.4),
                        alpha_max=deg2rad(60.0),
                        FILL=-9999999.0)
    dx = step(dem.x)
    # Smooth dem to get smooth alpha, smoother bands.  This is in
    # particular important when there is cross-flow bumpiness, such as
    # on Uaar.  However, it can also be bad.  YMMV, check!
    if window_dem_smooth>0
        dem = deepcopy(dem)
        fillmask = dem.v.!=FILL
        mask = fillmask .& glaciermask
        dem.v[:] = boxcar(dem.v, round(Int,window_dem_smooth/dx), mask, (!).(mask))
        dem.v[(!).(fillmask)] = FILL
    end
    # no FILL inside glaciermask
    @assert !any(dem.v[glaciermask].==FILL)
    @assert !any(isnan.(dem.v[glaciermask]))

    # 2D slopes
    ret_nans = false
    alpha2d = absslope(dem, glaciermask, ret_nans)

    bands, bandi = bin_grid(dem, binsize_or_bins, glaciermask, binround=binround)

    nb = length(bands)
    cellsz = step(dem.x)^2

    totalarea = sum(glaciermask)*cellsz

    malphas, widths, lengths, areas = (zeros(nb) for i=1:4)
    dzs = zeros(nb)
    for i=1:nb
        ind = bandi[i]
        if i!=nb
            dzs[i] = abs(bands[i+1]-bands[i])
        else
            # TODO: a hack. Fix
            dzs[i] = abs(bands[i]-bands[i-1])
        end
        # this is the critical step:
        malphas[i] = band_slope!(alpha2d[ind], i, alpha_min, alpha_max)
        areas[i] = length(ind)*cellsz
        lengths[i] = dzs[i]/tan(malphas[i])
        widths[i] = areas[i]/lengths[i]
    end
    # update missing bands
    for i=1:nb
        if malphas[i]==-9999
            ma1 = malphas[max(1,i-1)]
            ma2 = malphas[min(nb,i+1)]
            if ma1==-9999 && ma2==-9999
                error("Too many consecutive bands for which no slope could be calculated: $(max(1,i-1):min(nb,i+1))")
            elseif ma1==-9999
                malphas[i] = ma2
            elseif ma2==-9999
                malphas[i] = ma1
            else
                malphas[i] = 1/2*(ma1+ma2)
            end
            lengths[i] = dzs[i]/tan(malphas[i])
            widths[i] = areas[i]/lengths[i]
        end
    end

    # Smooth the width:
    # Note, this can make the malphas noisy!
    if window_width_smooth>0
        widths = boxcar(widths, round(Int, window_width_smooth/mean(dzs) ))
        for i=1:nb
            lengths[i] = areas[i]/widths[i]
            malphas[i] = atan(dzs[i]/lengths[i])
            @assert malphas[i]>=0
        end
    end

    x = vcat(0,cumsum(lengths))
    xmid = x[1:end-1] + diff(x)/2

    # tests
    if abs(totalarea-sum(areas))>1.0
        error("Something's amiss, sum of area of bands $(sum(areas)) not equal total area $totalarea,")
    end
    # check band length against diagonal
    tmp = sum(glaciermask,2)
    xextent = (findlast(tmp.>0)-findfirst(tmp.>0))*dx
    tmp = sum(glaciermask,1)
    yextent = (findlast(tmp.>0)-findfirst(tmp.>0))*dx
    box_diag = sqrt(xextent^2 + yextent^2)
    if abs(sum(lengths)-box_diag)/box_diag>0.4
        warn("Glacier length from might be wrong. Band-length: $(sum(lengths)/1e3)km, bounding-box diag: $(box_diag/1e3)km")
    end

    return bands, bandi, malphas, areas, lengths, widths, x, xmid, dem, alpha2d
end


"""
    band_slope!(alphas, bandnr, alpha_min=deg2rad(0.4), alpha_max=deg2rad(60.0))

Need to calculate a meaningful mean of the slopes in a elevation bin.
*This is tricky but critical!*

One check can be that all the bin-lengths should add up to the total
glacier length.

Input:
- alphas -- slope angles in one band (these are sorted in place, thus the ! in the function name)
- bandnr -- which band those alphas belong to (only used for error message)
- alpha_max, alpha_min -- maximal and minimal allowed slope
"""
function band_slope!(alphas, bandnr, alpha_min, alpha_max)
    # parameters
    ratio_fac = 2
    f_q5 = 0.05
    f_q20 = 0.2
    f_q80 = 0.8
    q_band = (0.55, 0.95)

    n = length(alphas)
    if n==0
        return -9999*one(alpha_min)
        # error("Band $bandnr has no elements!")
        # return deg2rad(45)
    end
    # magic slope calculation
    sort!(alphas)
    # calculate indices of 5th, 20th and 80th quantiles:
    iq5  = max(1,round(Int,n*f_q5))
    iq20 = max(1,round(Int,n*f_q20))
    iq80 = min(n,round(Int,n*f_q80))
    # angle of those quantiles:
    q5, q20, q80 = [max(rad2deg.(i),eps(alpha_min)) for i in (alphas[iq5], alphas[iq20], alphas[iq80])]
    # Now some of Matthias' magic:
    a = (q20/q80)*ratio_fac # 2x ratio of angles
    a = min(a, q_band[2])
    a = max(a, q_band[1]) # scaled to lie [0.55, 0.95]
    iq_magic = round(Int,n*a)  # set a "new" magic quantile to that value
    q_magic = rad2deg(alphas[iq_magic])
    # only use indices within those quantiles
    ind = q5 .<= rad2deg.(alphas) .< q_magic
    out = sum(ind)>1 ? mean(alphas[ind]) : mean(alphas)
    # limit alphas
    out = max(alpha_min, out)
    out = min(alpha_max, out)
    return out
end

"""
Bins a gird into bands.  Often used to bin a DEM into elevation bands.

- g -- to be binned ::Gridded or ::Matrix
- binsize_or_bins -- bin size or the bins (a Range)
- mask -- specify if not all locations of a gird should be binned.

KW:
- binround -- floor the bin-start using this many digits (see help of floor)

Return:
- bands -- a range of the bands, e.g. 0.0:10.0:100.0
- bandi -- a Vector{Vector{Int}} of length(bands) with each element
           containing the indices of cells in the band
"""
bin_grid(g::Gridded, binsize_or_bins, mask=BitArray([]); binround=_binround(binsize_or_bins)) =
    bin_grid(g.v, binsize_or_bins, mask; binround=binround)
function bin_grid(v::Matrix, binsize_or_bins, mask=BitArray([]); binround=_binround(binsize_or_bins))
    if isempty(mask)
        v = v
        ginds = 1:length(v)
    else
        @assert size(mask)==size(v)
        v = v[mask]
        ginds = find(mask[:])
    end
    nv = length(v)
    if isa(binsize_or_bins, Number) # i.e. a binsize
        mi, ma = minimum(v), maximum(v)
        if binsize_or_bins>=0
            binstart = floor(mi, binround)
            binend = floor(ma, binround) # better: `ceil(ma, binround) - binsize_or_bins` ?
        else
            binstart = ceil(ma, binround)
            binend = ceil(mi, binround) # better: `ceil(ma, binround) - binsize_or_bins` ?
        end
        @assert !isnan(binstart) && !isnan(binend)
        bins = binstart:binsize_or_bins:binend # these are the start of the bins
    else
        bins = binsize_or_bins
    end
    _bin_grid_kernel(bins, nv, v, ginds)
end
@inbounds function _bin_grid_kernel(bins, nv, v, ginds)
    # initialize output
    indices = Vector{Int}[]
    for b in bins
        ind = Int[]
        push!(indices, ind)
    end
    # fill it
    for j=1:nv
        if (bins[2]-bins[1])>0
            i = searchsortedlast(bins, v[j])
            i = i==0 ? 1 : i # if smaller then add to lowest bin
        else
            # https://github.com/JuliaLang/julia/issues/18653
            i = searchsortedlast(collect(bins), v[j], rev=true)
        end
        push!(indices[i], ginds[j])
    end
    return bins, indices
end


import Interpolations
"""
Bins a trajectory using a grid

- tr -- trajectory to be binned
- g -- g.v of grid used for binning
- binsize_or_bins -- bin size or the bins (a Range)
- mask -- specify if not all locations of a gird should be binned.]

KW:
- binround -- floor the bin-start using this many digits (see help of floor)
"""
function bin_traj(tr::Traj, g::Gridded, binsize_or_bins, mask=trues(size(g.v)); binround=_binround(binsize_or_bins))
    @assert size(mask)==size(g.v)

    demi  = Interpolations.interpolate((g.x, g.y), g.v, Interpolations.Gridded(Interpolations.Linear()) )
    maski  = Interpolations.interpolate((g.x, g.y), mask, Interpolations.Gridded(Interpolations.Constant()) ) # returns Int!

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

Map a field onto the (elevation) bands.  The field needs to have the
same size as the original binned-grid.

Input:
- bandi -- as returned by bin_grid
- field -- the field, either a Matrix or a Gridded
- fn -- the function to do the reduction with.  Default==mean
- fill -- fill value, if set, ignore those points

Output
- the value of the field on the bands.  If no values were found in a band,
  then return NaN.
"""
function map_onto_bands(bandi, field::Matrix, fn=mean, fill=nothing)
    resT = typeof(fn(field[bandi[1]])) # to get result type
    out = zeros(resT, length(bandi))
    for i in 1:length(bandi)
        count = 0
        for j=1:length(bandi[i])
            val = field[bandi[i][j]]
            if val!=fill
                # only use the ones which have no fill
                out[i] += val
                count+=1
            end
        end
        out[i] /=count
    end
    return out
end
map_onto_bands(bandi, field::Gridded, fn=mean, fill=NaN) = map_onto_bands(bandi, field.v, mean, fill)

"""
    map_back_to_2D(dims2d, bandi, field1d)

Maps 1D field back onto 2D.  More or less inverse of map_onto_bands.
"""
function map_back_to_2D(dims2d, bandi, field1d)
    out = zeros(eltype(field1d), dims2d)
    map_back_to_2D!(out, bandi, field1d)
    out
end
"""
    map_back_to_2D!(out2d, bandi, field1d)

Maps 1D field back onto 2D.  More or less inverse of map_onto_bands.
"""
function map_back_to_2D!(out, bandi, field1d)
    for (i,is) in enumerate(bandi)
        out[is] = field1d[i]
    end
    nothing
end


"""
    bins2matrix(g::Union{Gridded,Matrix}, bands, bandi) -> binmat

Return a matrix which gives the bin-number of each its (i,j) locations.
Locations not binned (i.e. masked) are ==0.
"""
bins2matrix(g::Gridded, bands, bandi) = bins2matrix(g.v, bands, bandi)
function bins2matrix(g::Matrix, bands, bandi)
    out = zeros(Int, size(g))
    for (n,b) in enumerate(bandi)
        for i in b
            out[i] = n
        end
    end
    return out
end


"""
    bandi_for_other_grid(bands, bandi, binmat, g::Gridded,
                         othergrid::Gridded, othermask=trues(size(othergrid.v))
    bandi_for_other_grid(bands, bandi, g::Gridded,
                         othergrid::Gridded, othermask=trues(size(othergrid.v)))

Returns vector of indices (bandi) to map a different grid onto the
bands encoded in `binmat` (or `bands, bandi`) and grid `g`.  It only
maps points onto bands which are within the mask applied to generate
the bands.  Additionally & optionally, a mask for the othergrid can
also be given.  Return:

    bandi
"""
function bandi_for_other_grid(bands, bandi::Vector{Vector{Int}}, g::Gridded,
                              othergrid::Gridded, othermask=trues(size(othergrid.v)))
    binmat=bins2matrix(g, bands, bandi)
    bandi_for_other_grid(bands, bandi, binmat, g, othergrid, othermask)
end
function bandi_for_other_grid(bands, bandi, binmat::Matrix{Int}, g::Gridded,
                              othergrid::Gridded, othermask=trues(size(othergrid.v)) )
    og = othergrid
    @assert size(othergrid)==size(othermask)
    if g.x!=og.x || g.y!=og.y
        bandi_ = [Int[] for i=1:length(bands)]
        dims = size(og.v)
        itpm = Interpolations.interpolate((g.x, g.y), binmat,
                            Interpolations.Gridded(Interpolations.Constant()) );
        itpm = Interpolations.extrapolate(itpm, 0);
        for j=1:size(og.v,2)
            for i=1:size(og.v,1)
                if othermask[i,j]
                    ind = convert(Int, itpm[og.x[i], og.y[j]])
                    if ind>0
                        push!(bandi_[ind], sub2ind(dims, i, j))
                    end
                end
            end
        end
    else
        bandi_ = deepcopy(bandi)
    end
    return bandi_
end


##################
# Extrapolation
##################

"""
To specify which edge of a cell is meant.  `_noedge` can be used if none
is used.
"""
@enum Loc _noedge=0 left=1 right=2 lower=3 upper=4

"Orientation of a line"
@enum Orientation nohand=0 lefthand=1 righthand=2

"""
One cell-edge of a regular gird of cells.

TODO:
Arguably not the best datastructure for what is done below.
"""
immutable Edge
    i::Int # cell ind
    j::Int # cell ind
    loc::Loc
end
const noedge = Edge(-1,-1,_noedge)
Base.show(io::IO,e::Edge) = println(io, "($(e.i),$(e.j),$(e.loc))")

"""
    get_nodes(e::Edge)

Returns the start and end nodes of the edge on the staggered grid.

The nodes are returned such that the cell is on the right of the edge.
"""
function get_nodes(e::Edge)
    i,j,loc = e.i, e.j, e.loc
    if loc==left
        return (i,j), (i,j+1)
    elseif loc==right
        return (i+1,j+1), (i+1,j)
    elseif loc==lower
        return (i+1,j), (i,j)
    elseif loc==upper
        return (i,j+1), (i+1,j+1)
    else
        error("!")
    end
end

"""
    orientation(e1::Edge) -> always left
    orientation(e1::Edge, e2::Edge)

Return whether the cells are on left, right or both of the two cells.
If not connected throws an error or return `_noedge`, depending on
`throwerror` flag.

"""
function orientation(e1::Edge, e2::Edge, throwerror=true)::Orientation
    n1,n2 = get_nodes(e1)
    m1,m2 = get_nodes(e2)
    if n2==m1
        return righthand
    elseif n1==m2
        return lefthand
    elseif n1==m1 || n2==m2
        return nohand
    elseif throwerror
        error("The two edges:\n $e1 $e2 are not connected")
    else
        return nohand
    end
end
orientation(e1::Edge)::Orientation = lefthand


"""
A line made up of a continuous (and sorted) collection of edges.

Also, the line has an orientation: the (i,j) cell lies on the "right" of the edge.
"""
immutable Line
    edges::Vector{Edge}
end
#Base.show(io::IO,l::Line) = show(io, l.edges)

Base.getindex(l::Line,i) = l.edges[i]
Base.length(l::Line) = length(l.edges)

"""
    orientation(l::Line)

Returns the orientation of a line: `left` means cell is on left of line.
"""
orientation(l::Line)::Orientation = orientation(l.edges)
function orientation(l::Vector)::Orientation
    if length(l)==0
        error("Line has length 0")
    end
    if length(l)==1
        # println("Line has length $(length(l)).  Too short to establish orientation.")
        return orientation(l[1]) # returns `left`
    end
    # establish orientation of first two segments
    ori = orientation(l[1],l[2])
    for i=3:length(l)
        ori2 = orientation(l[i-1],l[i])
        if ori2!=ori
            println("""
                    Line changes orientation between edge # $(i-2) and $i from $ori to $ori2.
                    e0: $(l[i-2])
                    e1: $(l[i-1])
                    e2: $(l[i])
                    """)
            return _noedge
        end
    end
    return ori
end

"""
    next_edge!(edges::Set{Edge}, edge::Edge, testori::Orientation)

Return an adjacent edge of input edge such that the orientation of the
two edges is equal `testori`.  Pops returned edge off `edges`.  If no
edge is found, return `Edge(-1,-1,nodege)`.

"""
function next_edge!(edges::Set{Edge}, edge::Edge, testori::Orientation)
    edge.loc==_noedge && error("Not a proper edge: $edge")
    locs = (left,right,lower,upper)
    for i=[0,-1,1],j=[0,-1,1]
        for loc in locs
            test_edge = Edge(edge.i+i,edge.j+j,loc)
            if orientation(edge, test_edge,false)==testori && (test_edge in edges)
                return pop!(edges, test_edge)
            end
        end
    end
    # nothing found
    return noedge
end

"""
    get_edges_on_boundary(bands, bandi, binmat, landmask=nothing) -> es

Return a list of sets of cell-edges which are at the boundary.

    typeof(es) == Dict{Tuple{Int,Int},Set{Edge}}()

which maps (bandnr,otherbandnr) => Set of edges

"""
function get_cells_on_boundary(bands, bandi, binmat, landmask=nothing)
    dims = size(binmat)

    if landmask!=nothing
        # encode sea cells in binmat
        binmat = copy(binmat)
        binmat[(landmask.==0) .& (binmat.==0)] = -1
    end

    # The boundaries of all bands:
    # (bandnr,otherbandnr) => Set of edges
    edges_at_boundaries = Dict{Tuple{Int,Int},Set{Edge}}()

    # Loop to identify all boundary cells of each band and determine
    # which cell-edges are on the boundary:
    for (ib,bup) in enumerate(bands)
        for I in bandi[ib]
            i,j = ind2sub(dims, I)
            loc = 1
            # left-right
            for ii=-1:2:1
                if 1<=i+ii<=dims[1] # skip if outside of binmat
                    iother = binmat[i+ii,j]
                    if iother!=ib
                        push!(get!(edges_at_boundaries, (ib, iother), Set{Edge}()),
                              Edge(i,j,loc))
                    end
                end
                loc+=1
            end
            # lower-upper
            for jj=-1:2:1
                if 1<=j+jj<=dims[2] # skip if outside of binmat
                    iother = binmat[i,j+jj]
                    if iother!=ib
                        push!(get!(edges_at_boundaries, (ib, iother), Set{Edge}()),
                              Edge(i,j,loc))
                    end
                end
                loc+=1
            end
        end
    end
    return edges_at_boundaries
end

"""
    calc_boundaries(bands, bandi, binmat, landmask=nothing)

Calculate interface between bands and also to outside (bandnr==0) and
sea-cells (bandnr==-1).

Returns `boundaries` which is `bandnr => otherbandnr => Vector{Line}`
i.e. the collection of all boundary lines for each
`(bandnr,otherbandnr)`.  Each line has the elevation band on its
*left*.  For cells outside the glacier `otherbandnr==0` except if
'landmask` is given, then sea-cells have `otherbandnf==-1`.

Of type `Vector{Dict{Int,Vector{Line}}}()`.
"""
function calc_boundaries(bands, bandi, binmat, landmask=nothing)
    dims = size(binmat)

    if landmask!=nothing
        # encode sea cells in binmat
        binmat = copy(binmat)
        binmat[(landmask.==0) .& (binmat.==0)] = -1
    end

    # The boundaries of all bands:
    # (bandnr,otherbandnr) => Set of edges
    edges_at_boundaries = get_cells_on_boundary(bands, bandi, binmat, landmask)

    # bandnr => otherbandnr => Vector{Line}
    boundaries = [Dict{Int,Vector{Line}}() for i=1:length(bands)]

    for ((band,otherband), edges) in edges_at_boundaries

        # # delete:
        # band, otherband = 2,0
        # edges = boundaries[(band,otherband)]

        output = Line[]

        # Take one starting point and go in both directions until no
        # more neighbors found.
        while length(edges)>0
            firstedge = first(edges)
            delete!(edges, firstedge)
            out = [Edge[], Edge[]]
            for ii=1:2 # first lefthand Orientation(1), then righthand Orientation(2)
                curedge = firstedge
                while curedge!=noedge
                    push!(out[ii], curedge)
                    curedge = next_edge!(edges, curedge, Orientation(ii))
                end
            end
            prepend!(out[2], reverse(out[1][2:end]))
            # make sure the elevation band is on the line's left:
            ori = orientation(out[2])
            if ori == righthand
                push!(output, Line(reverse(out[2])))
            elseif ori==lefthand
                push!(output, Line(out[2]))
            elseif ori==nohand
                error("""
                      Constructed line has not a consistent orientation.
                      Band: $band
                      Otherband: $otherband
                      """
                      )
            end
        end
        boundaries[band][otherband] = output
    end

    return boundaries
end

"""
    calc_fluxdir_source(l::Line, ux, uy, window, bands) -> fluxdir,ii,jj,fluxdirx,fluxdiry

Calculated the flux-direction on a Line by taking a running average
over some part of it.  Also returns the flux to use and the upstream
cell `(ii,jj)`.  The flux across the k-th edge is then given by:

    fluxdir[k] * h[ii[k],jj[k]] * u[ii[k],jj[k]]

The flux across the line is

    sum(fluxdir .* h_line .* u_line)

Notes:
- this is not 100% correct but probably close enough.
  - now the averaging is over all edges.  Arguably it could be done over cells?
"""
function calc_fluxdir(l::Line, ux, uy, window, dx, bands)
    edges = l.edges
    ne = length(edges)
    # x,y flux-direction on the cell
    fluxdirx = zeros(ne)
    fluxdiry = zeros(ne)
    fluxdir = zeros(ne) # the flux which goes with the edge
    # indices of upstream cell
    ii = zeros(Int,ne)
    jj = zeros(Int,ne)

    R = CartesianRange(size(edges))
    I1, Iend = first(R), last(R)
    for I in R
        fx, fy = 0.0, 0.0
        for J in CartesianRange(max(I1, I-I1*window), min(Iend, I+I1*window))
            i,j,loc = edges[J].i, edges[J].j, edges[J].loc
            # average ux,uy in the two cells
            if loc==left
                fx += (ux[i-1,j] + ux[i,j])
                fy += (uy[i-1,j] + uy[i,j])
            elseif loc==right
                fx += (ux[i+1,j] + ux[i,j])
                fy += (uy[i+1,j] + uy[i,j])
            elseif loc==lower
                fx += (ux[i,j-1] + ux[i,j])
                fy += (uy[i,j-1] + uy[i,j])
            elseif loc==upper
                fx += (ux[i,j+1] + ux[i,j])
                fy += (uy[i,j+1] + uy[i,j])
            else
                error()
            end
        end
        norm = sqrt(fx^2+fy^2)
        fx, fy = fx/norm, fy/norm
        fluxdirx[I] = fx
        fluxdiry[I] = fy
        # figure out upstream cell
        i,j,loc = edges[I].i, edges[I].j, edges[I].loc
        ii[I],jj[I],fluxdir[I] =
            if loc==left
                fx<0 ? (i,j,fx*dx) : (i-1,j,fx*dx)
            elseif loc==right
                fx<0 ? (i+1,j,-fx*dx) : (i,j,-fx*dx)
            elseif loc==lower
                fy<0 ? (i,j,fy*dx) : (i,j-1,fy*dx)
            else #if loc==upper
                fy<0 ? (i,j+1,-fy*dx) : (i,j,-fy*dx)
            end
    end
    # for reverse ordered bands, need a sign flip
    sig = sign(bands[2]-bands[1])

    return sig*fluxdir,ii,jj,fluxdirx,fluxdiry
end

"""
    calc_u(q1d, boundaries, thick, alpha, ux, uy, dx, mask, bands, lengths,
           flux_dir_window, window_frac=1.5,
           x=nothing, y=nothing) # these are only needed for plotting
    -> u2d, u2d_at_bands, scaling_factors_1d, mask_u2d_at_bands

Calculates a 2D field of depth averaged ice flow speed `ubar` which has the
same flux across elevation band boundaries as the supplied 1D flux
`q1d`.  Only `q1d` is a flow-band variable.  The flux across elevation
bands is calculated with `calc_fluxdir`.

Note:

- An assumption about the distribution of the u is needed, which will
  be scaled to conform to mass conservation.  At the moment a function
  `thick[i,j].^u_exp` is used.

TODO:
- pre-calculate the boxcar filter
"""
function calc_u(q1d, boundaries, u_trial, thick,
                ux, uy, dx, mask, bands, lengths,
                flux_dir_window, # in [m]
                boxcarM::AbstractMatrix;
                plotyes=false,
                x=nothing, y=nothing)
    ubar_, ubar, facs, mask_ubar = _calc_u(q1d, boundaries, u_trial, thick,
                                    ux, uy, dx, mask, bands,
                                    flux_dir_window, # in [m]
                                           plotyes,x,y)
    # below type assertion is needed for type-stability ?!
    return VAWTools.apply_boxcar_matrix(boxcarM, ubar_)::Array{eltype(q1d),2}, ubar, facs, mask_ubar
end
function calc_u(q1d, boundaries, u_trial, thick,
                ux, uy, dx, mask, bands, lengths,
                flux_dir_window, # in [m]
                window_frac;
                plotyes=false,
                x=nothing, y=nothing)
    ubar_, ubar, facs, mask_ubar_ = _calc_u(q1d, boundaries, u_trial, thick,
                                    ux, uy, dx, mask, bands,
                                    flux_dir_window, # in [m]
                                    plotyes,x,y)
    return boxcar(ubar_, Int((window_frac*maximum(lengths))÷dx)+1, mask_ubar_, (!).(mask) ), ubar, facs, mask_ubar_
end
# this helper function is needed for type stability
function _calc_u(q1d, boundaries, u_trial, thick,
                ux, uy, dx, mask, bands,
                flux_dir_window,
                plotyes,
                x, y) # these are only needed for plotting
    #plotyes && figure()
    dims = size(mask)

    # this calculates the u at all elevation band boundaries:
    ubar = zeros(dims)*NaN
    facs = Float64[]
    for (ib,bnd) in enumerate(boundaries)
        if ib<length(bands)
            # Find the receiving band, can be a number further than 1.
            ibb = 0
            for ib_=ib+1:length(bands)
                if haskey(bnd, ib_)
                    ibb = ib_
                    break
                end
            end
            ibb==0 && error("No band below band $ib, but $ib is not bottom band!  This means that the domain is likely disjoint.")
            bb=bnd[ibb]
        else # outflow at terminus
            # TODO: this probably needs updating where several elevation bands contribute (tide-water)
            # -> no, only ever the last band does outflow
            @show keys(bnd)
            bb = get(bnd, -1, bnd[0]) # if sea-terminating use those edges.
        end
        # loop over segments
        ffs = Float64[]
        is = Int[]
        js = Int[]
        for l in bb
            #@assert orientation(l)
            ff, fi, fj, fx, fy = calc_fluxdir(l, ux, uy, Int(flux_dir_window÷dx), dx, bands)
            # if plotyes
            #     quiver(x[fi],y[fj],fx,fy)
            # end
            append!(is, fi)
            append!(js, fj)
            append!(ffs, ff)
        end
        u_t = [u_trial[i,j] for (i,j) in zip(is,js)]
        h_line = [thick[i,j] for (i,j) in zip(is,js)]
        q_trial = sum(ffs.*u_t.*h_line)
        fac = max(q1d[ib],0)/q_trial
        push!(facs, fac)
        u = u_t*fac
        # m,ij = findmin(u)
        #@assert all(u.>=0) "$i, $ij, $(u_trial[ij]), $(ff[ij]), $(q_trial), $(q1d[i]), $(h_line[ij])"
        for (n,(i,j)) in enumerate(zip(is,js))
            if u[n]<0;
                println("Flow speed<0: $(u[n]) at location ($i,$j).  This should not happen!")
                u[n]=0
            end
            ubar[i,j] = u[n]
        end
    end
    mask_ubar_ = mask .& ((!).(isnan.(ubar))) # location of all cells for which `ubar` was calculated
    ubar_ = copy(ubar)
    ubar_[isnan.(ubar_)] = 0
    # if plotyes
    #     # imshow(binmat',origin="lower", extent=(x[1],x[end],y[1],y[end]), cmap="flag"); colorbar();
    #     imshow(ubar',origin="lower", extent=(x[1],x[end],y[1],y[end]),); colorbar(); clim(0,50)
    # end

    return ubar_, ubar, facs, mask_ubar_
end

"""
    get_iv_boxcar_M(F, dem, mask, bands, bandi, lengths, iv_window_frac)

Precalculate the boxcar operator for the IV calculation (this is the
most expensive part).
"""
function get_iv_boxcar_M(F, dem, mask, bands, bandi, lengths, iv_window_frac)
    ux,uy = (-).(gradient3by3(dem, mask))
    binmat = bins2matrix(dem, bands, bandi)
    boundaries = calc_boundaries(bands,bandi,binmat)
    q1d = ones(bands)
    u_trial = ones(dem.v)
    thick = u_trial
    dx = step(dem.x)
    flux_dir_window = 2
    ubar_, ubar, facs, mask_ubar_ = _calc_u(q1d, boundaries, u_trial, thick,
                                            ux, uy, dx, mask,
                                            bands,
                                            flux_dir_window,
                                            false,nothing,nothing)

    return boxcar_matrix(F, Int((iv_window_frac*maximum(lengths))÷dx)+1, mask_ubar_, (!).(mask)),
           boundaries, ux, uy
end


#################
# Plotting
################
"""
    plot_bands(dem, bands, bandi; bands2plot=1:length(bands))

Plots the bands in 2D.  For length(bands2plot)<=16 the default
colorscale will show one color per band.

Needs Plots.jl imported in Main (REPL)
"""
function plot_bands(dem, bands, bandi; bands2plot=1:length(bands))
    binmat = Float64.(VAWTools.bins2matrix(dem, bands, bandi))
    for i in eachindex(binmat)
        if !(binmat[i] in bands2plot)
            binmat[i] = NaN
        end
    end
    Main.contourf(dem.x,dem.y,binmat',aspect_ratio=:equal)
end
