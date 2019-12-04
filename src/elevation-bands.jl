##############
# Bands
##############
#
# TODO: add extrapolation too (second half of make_bands)

export map_onto_bands, make_1Dglacier, map_back_to_2D, map_back_to_2D!

# used to round (up or down) the binsize to the next decimal place
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
                        min_bin_number_ends=0,
                        min_bands=4,
                        window_dem_smooth=0.0,
                        window_width_smooth=0.0,
                        alpha_min=deg2rad(0.4),
                        alpha_max=deg2rad(60.0),
                        FILL=-9999999.0,
                        verbose=true)
    dx = step(dem.x)
    # Smooth dem to get smooth alpha, smoother bands.  This is in
    # particular important when there is cross-flow bumpiness, such as
    # on Uaar.  However, it can also be bad.  YMMV, check!
    if window_dem_smooth>0
        dem = deepcopy(dem)
        nofillmask = dem.v.!=FILL
        mask = nofillmask .& glaciermask
        dem.v[:] = boxcar(dem.v, round(Int,window_dem_smooth/dx), mask, (!).(mask))
        dem.v[(!).(nofillmask)] = FILL
    end
    # no FILL inside glaciermask
    @assert !any(dem.v[glaciermask].==FILL)
    @assert !any(isnan.(dem.v[glaciermask]))

    # 2D slopes
    ret_nans = false
    alpha2d = absslope(dem, glaciermask, ret_nans)

    bands, bandi = bin_grid(dem, binsize_or_bins, glaciermask,
                            binround=binround, min_bin_number_ends=min_bin_number_ends,
                            min_bands=min_bands)
    @assert length(bands)>=min_bands "Need at least $min_bands elevation bins, only got $(length(bands))"

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
            # make sure length is zero when area is zero
            lengths[i] = areas[i]==widths[i]==0 ? zero(lengths[i]) : areas[i]/widths[i]
            malphas[i] = atan(dzs[i]/lengths[i])
            @assert malphas[i]>=0
        end
    end

    x = vcat(0,cumsum(lengths))
    xmid = x[1:end-1] + diff(x)/2

    # tests
    # if abs(totalarea-sum(areas))>1.0
    #     error("Something's amiss, sum of area of bands $(sum(areas)) not equal total area $totalarea,")
    # end
    # check band length against diagonal
    tmp = sum(glaciermask,2)
    xextent = (findlast(tmp.>0)-findfirst(tmp.>0))*dx
    tmp = sum(glaciermask,1)
    yextent = (findlast(tmp.>0)-findfirst(tmp.>0))*dx
    box_diag = sqrt(xextent^2 + yextent^2)
    if verbose && abs(sum(lengths)-box_diag)/box_diag>0.4
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
    # calculate indices of 5th, 20th and 80th percentiles:
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
- min_bin_number_ends -- minimum number of elements in the uppermost and lowermost
                         band.  If below, merge those cells into the first viable band.
- min_bands -- minimum number of bands produced (3).  Not integrated with min_bin_number_ends!

Return:
- bands -- a range of the bands, e.g. 0.0:10.0:100.0
- bandi -- a Vector{Vector{Int}} of length(bands) with each element
           containing the indices of cells in the band
"""
function bin_grid(g::Gridded, binsize_or_bins, mask=BitArray([]);
                  binround=_binround(binsize_or_bins), min_bin_number_ends=0, min_bands=4)
    bin_grid(g.v, binsize_or_bins, mask;
             binround=binround, min_bin_number_ends=min_bin_number_ends,
             min_bands=min_bands)
end
function bin_grid(v::Matrix, binsize::Number, mask=BitArray([]);
                  binround=_binround(binsize), min_bin_number_ends=0, min_bands=4)/
    if isempty(mask)
        v = v
        ginds = 1:length(v)
    else
        @assert size(mask)==size(v)
        v = v[mask]
        ginds = find(mask[:])
    end
    nv = length(v)
    bins = 1.0:-1
    while length(bins)<=min_bands
        mi, ma = minimum(v), maximum(v)
        if binsize>=0
            binstart = floor(mi, binround)
            binend = floor(ma, binround) # better: `ceil(ma, binround) - binsize` ?
        else
            binstart = ceil(ma, binround)
            binend = ceil(mi, binround) # better: `ceil(ma, binround) - binsize` ?
        end
        @assert !isnan(binstart) && !isnan(binend)
        bins = binstart:binsize:binend # these are the start of the bins

        # decreas binsize for next round
        binsize = step(bins)/2
    end
    return _bin_grid_kernel(bins, nv, v, ginds, min_bin_number_ends)
end
function bin_grid(v::Matrix, bins::AbstractRange, mask=BitArray([]); min_bin_number_ends=0, kw...)
    if isempty(mask)
        v = v
        ginds = 1:length(v)
    else
        @assert size(mask)==size(v)
        v = v[mask]
        ginds = find(mask[:])
    end
    nv = length(v)
    _bin_grid_kernel(bins, nv, v, ginds, min_bin_number_ends)
end

@inbounds function _bin_grid_kernel(bins, nv, v, ginds, min_bin_number_ends)
    # initialize output
    bandi = Vector{Int}[]
    for b in bins
        ind = Int[]
        push!(bandi, ind)
    end
    # fill it
    for j=1:nv
        if step(bins)>0
            i = searchsortedlast(bins, v[j])
            i = i==0 ? 1 : i # if smaller then add to lowest bin
        else
            # https://github.com/JuliaLang/julia/issues/18653
            i = searchsortedlast(collect(bins), v[j], rev=true)
            i = i==0 ? 1 : i # if smaller then add to highest bin (i.e. bins[1])
        end
        push!(bandi[i], ginds[j])
    end
    # remove top and bottom bins if too few elements
    inds2pop = [1,length(bandi)] # remove up to and including these bandi
    num = [0,0]
    for i = 1:length(bandi)
        inds2pop[1] = i
        if length(bandi[i])+num[1]>=min_bin_number_ends
            break
        end
        num[1] +=length(bandi[i])
    end
    for i = length(bandi):-1:1
        inds2pop[2] = i
        if length(bandi[i])+num[2]>=min_bin_number_ends
            break
        end
        num[2] +=length(bandi[i])
    end
    # add the dropped cells to the next/previous band
    append!(bandi[inds2pop[1]], vcat(bandi[1:inds2pop[1]-1]...))
    append!(bandi[inds2pop[2]], vcat(bandi[inds2pop[2]+1:end]...))
    return bins[inds2pop[1]:inds2pop[2]], bandi[inds2pop[1]:inds2pop[2]]
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
    _bin_grid_kernel(bins, nv, v, ginds, 0)
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
    bandi_for_other_grid(bands, bandi, g::Gridded,
                         othergrid::Gridded, othermask=trues(size(othergrid.v))
    bandi_for_other_grid(bands, bandi, g::Gridded,
                         othergrid::Gridded, othermask=trues(size(othergrid.v)))

Returns vector of indices (bandi) to map a different grid (othergird) onto the
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
        s2i = LinearIndices(og.v)
        itpm = Interpolations.interpolate((g.x, g.y), binmat,
                            Interpolations.Gridded(Interpolations.Constant()) );
        itpm = Interpolations.extrapolate(itpm, 0);
        for j=1:size(og.v,2)
            for i=1:size(og.v,1)
                if othermask[i,j]
                    ind = convert(Int, itpm(og.x[i], og.y[j]))
                    if ind>0
                        push!(bandi_[ind], s2i[i, j])
                    end
                end
            end
        end
    else
        bandi_ = deepcopy(bandi)
    end
    return bandi_
end


#####################
# Extrapolation of IV
#####################

"""
To specify which edge of a cell is meant.  `_noedge` can be used if none
is used.
"""
@enum Loc _noedge=0 left=1 right=2 lower=3 upper=4

"Orientation of a line"
@enum Orientation nohand=0 lefthand=1 righthand=2

"""
One cell-edge of a regular gird of cells.  Specified as a cell
and which edge of the four edges of a cell.

TODO:
Arguably not the best datastructure for what is done below.
"""
struct Edge
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
A line made up of a continuous (and sorted) collection of edges.  Note
that the edges themselves also have a cell associated with them.  So,
a Line also a collection of (possibly repeated) cells.

The line has an orientation: when traversing the Line from 1 to n, all
cells are on the right of their edge.
"""
struct Line
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
    get_ux_uy(dem, mask)

The direction and magnitude of steepest descent on each point within the mask.
"""
function get_ux_uy(dem, mask)
    ux,uy = (-).(gradient3by3(dem, mask))
    return ux, uy
end

"""
    get_edges_on_boundary(bands, bandi, binmat, landmask=nothing) -> es

Return a list of sets of cell-edges which are at the boundary.

    typeof(es) == Dict{Tuple{Int,Int},Set{Edge}}()

which maps (bandnr,otherbandnr) => Set of edges

"""
function get_cells_on_boundary(bands, bandi, binmat, landmask=nothing)
    error("Not updated to Julia 1.0 yet")
    dims = size(binmat)
    i2s = CartesianIndices(binmat)

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
            i,j = i2s[I].I
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
    calc_fluxdir(l::Line, ux, uy, window, dx, bands) -> flux,ii,jj,fluxdirx,fluxdiry

Calculated the flux-direction on a Line by taking a running average of (ux,uy)
over some part of it of length window.  Also returns the flux to use and the upstream
cell `(ii,jj)`.  The flux across the k-th edge is then given by:

    flux[k] * h[ii[k],jj[k]] * u[ii[k],jj[k]]

where h is the thickness and u the depth averaged velocity.

The flux across the line is

    sum(flux .* h_line .* u_line)

Input:
- l::Line
- ux,uy: unsmoothed velocity field (calculated as fall-line of surface DEM with get_ux_uy)
- window : length of smoothing window (along band) in number of edges
- dx : grid spacing
- bands::AbstractRange : the elevation bands

Output:
- flux : the flux
- ii, jj : indices of upstream cell
- fluxdirx,fluxdiry : flux direction at upstream cell ii, jj

Notes:
- this is not 100% correct but probably close enough.
- now the averaging is over all edges.  Arguably it could be done over cells?
"""
function calc_fluxdir(l::Line, ux, uy, window::Int, dx, bands)
    edges = l.edges
    ne = length(edges)
    # flux in x and y-direction on the cell associated with a Line edge:
    fluxdirx = zeros(ne)
    fluxdiry = zeros(ne)
    # the flux which goes with the edge:
    flux = zeros(ne)
    # indices of upstream cell
    ii = zeros(Int,ne)
    jj = zeros(Int,ne)

    R = CartesianIndices(size(edges))
    I1, Iend = first(R), last(R)
    for I in R
        fx, fy = 0.0, 0.0
        for J in CartesianIndices(max(I1, I-I1*window), min(Iend, I+I1*window))
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
        ii[I],jj[I],flux[I] =
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

    return sig*flux,ii,jj,fluxdirx,fluxdiry
end

"""
    calc_u(q1d, boundaries, u_trial, thick,
                ux, uy, dx, mask, bands, bandi, lengths,
                flux_dir_window, # in [m]
                boxcarM::AbstractMatrix or window_frac,
                scale_u=ones(q1d);
                plotyes=false,
                x=nothing, y=nothing)

Calculates a 2D field of depth averaged ice flow speed `ubar` which has the
same flux across elevation band boundaries as the supplied 1D flux
`q1d`.  Only `q1d` is a flow-band variable.  The flux across elevation
bands is calculated with `calc_fluxdir`.

Input:
- q1d: 1d flux (m/a) on elevation bands
- boundaries: from calc_boundaries
- u_trial: 2D velocity field which has the right shape (say cross-valley) but not the right magnitude
- thick: ice thickness
- ux,uy: unsmoothed velocity field (calculated as fall-line of surface DEM with get_ux_uy)
- dx: grid size
- mask: true where glacier is
- bands : elevation bands
- bandi
- lengths: length of each elevation band
- flux_dir_window: as used in calc_fluxdir
- boxcarM or window_frac:
  - pre-calculated boxcar operator
  - window_frac: window over which to smooth as fraction of the maximum(length)
- scale_u (ones(q1d)): if provided scale output by this much -> so make surface speeds

Output:
- ubar2d: the smoothed IV
- ubar: the IV at the elevation-band boundaries
- mask_ubar: true where ubar contains a value
- facs: scaling factor in 1D

Note:
- An assumption about the distribution of the u is needed, which will
  be scaled to conform to mass conservation.  At the moment a function
  `thick[i,j].^u_exp` is used.
"""
function calc_u(q1d, boundaries, u_trial, thick,
                ux, uy, dx, mask, bands, bandi, lengths,
                flux_dir_window, # in [m]
                filter, scale_u=ones(q1d);
                plotyes=false,
                x=nothing, y=nothing)
    u, mask_u, facs = _calc_u(q1d, boundaries, u_trial, thick,
                                    ux, uy, dx, mask, bands,
                                    flux_dir_window, # in [m]
                                    plotyes, x, y, scale_u)
    u_ = copy(u)
    u_[isnan.(u_)] = 0  # remove NaNs to allow filtering

    # type assertion is needed for type-stability
    out::Array{eltype(q1d),2} = if filter isa AbstractMatrix
        VAWTools.apply_boxcar_matrix(filter, u_)
    else
        boxcar(u_, Int((filter*maximum(lengths))÷dx)+1, mask_u, (!).(mask) )
    end

    # Make depth-averaged flow speed into a surface flow speed.
    scale_u2d = map_back_to_2D(size(out), bandi, scale_u)
    out = out .* scale_u2d

    return out, u, mask_u, facs
end
# this helper function is needed for type stability.
function _calc_u(q1d, boundaries, u_trial, thick,
                 ux, uy, dx, mask, bands,
                 flux_dir_window,
                 plotyes,
                 x, y, scale_u) # these are only needed for plotting
    #plotyes && figure()
    dims = size(mask)

    # calculate the u at all elevation band boundaries:
    u2d = zeros(dims)*NaN
    facs = Float64[] # scaling factor in 1D
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
            # This probably needs updating where several elevation bands contribute (tide-water)
            # -> no, only ever the last band does outflow in our 1D model
            #
            # This errors at times, for example on RGI60-14.11814 (which is not sea-terminating)
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
                println("Flow speed<0: $(u[n]) at location ($i,$j).  This should not happen!  Setting to zero.")
                u[n]=0
            end
            # also scale u:
            u2d[i,j] = u[n] * scale_u[ib]
        end
    end
    mask_u = mask .& ((!).(isnan.(u2d))) # location of all cells for which `u` was calculated
    # if plotyes
    #     # imshow(binmat',origin="lower", extent=(x[1],x[end],y[1],y[end]), cmap="flag"); colorbar();
    #     imshow(u',origin="lower", extent=(x[1],x[end],y[1],y[end]),); colorbar(); clim(0,50)
    # end

    return u2d, mask_u, facs
end

"""
    get_iv_boxcar_M(F, dem, mask, bands, bandi, lengths, iv_window_frac)

Precalculate the boxcar operator for the IV calculation (this is the
most expensive part).
"""
function get_iv_boxcar_M(F, dem, mask, bands, bandi, lengths, iv_window_frac)
    ux,uy = get_ux_uy(dem, mask)
    binmat = bins2matrix(dem, bands, bandi)
    boundaries = calc_boundaries(bands,bandi,binmat)
    q1d = ones(bands)
    u_trial = ones(dem.v)
    thick = u_trial
    dx = step(dem.x)
    flux_dir_window = 2
    mask_u = _calc_u(q1d, boundaries, u_trial, thick,
                     ux, uy, dx, mask,
                     bands,
                     flux_dir_window,
                     false,nothing,nothing,ones(q1d))[2]

    return boxcar_matrix(F, Int((iv_window_frac*maximum(lengths))÷dx)+1, mask_u, (!).(mask)),
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
    Main.PyPlot.contourf(dem.x,dem.y,binmat',aspect_ratio=:equal)
end
