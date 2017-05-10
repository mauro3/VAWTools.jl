##############
# Bands
##############
#
# TODO: add extrapolation too (second half of make_bands)

export map_onto_bands, make_1Dglacier, map_back_to_2D, map_back_to_2D!

_binround(binsize_or_bins) = -floor(Int, log10(binsize_or_bins))

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
- bands -- elevation bands
- bandi -- linear index into dem which assigns each cell to a elevation band
- alphas, areas, lengths, widths, x, xmid -- elevation band slope, area, etc
- dem -- the used DEM, a possibly smoothed version of the input DEM
- alpha2d -- slopes at each point of `dem`
"""
function make_1Dglacier(dem::Gridded, binsize_or_bins, glaciermask=BitArray([]);
                        binround=_binround(binsize_or_bins),
                        window_dem_smooth=0.0,
                        window_width_smooth=0.0,
                        alpha_min=deg2rad(0.4),
                        alpha_max=deg2rad(60.0),
                        FILL=-9999999.0)


    # Smooth dem to get smooth alpha, smoother bands.  This is in
    # particular important when there is cross-flow bumpiness, such as
    # on Uaar.  However, it can also be bad.  YMMV, check!
    if window_dem_smooth>0
        dem = deepcopy(dem)
        dx = step(dem.x)
        # TODO: should smoothing only happen over glacier?  But then
        # absslope does not work.
        fillmask = dem.v.!=FILL
        dem.v[:] = boxcar(dem.v, round(Int,window_dem_smooth/dx), fillmask)
        dem.v[!fillmask] = FILL
    end
    # 2D slopes
    alpha2d = absslope(dem)

    bands, bandi = bin_grid(dem, binsize_or_bins, glaciermask, binround=binround)
    bandsize = step(bands)

    nb = length(bands)
    dz = step(bands)
    cellsz = step(dem.x)^2

    totalarea = sum(glaciermask)*cellsz

    malphas, widths, lengths, areas = (zeros(nb) for i=1:4)
    for i=1:nb
        ind = bandi[i]
        # this is the critical step:
        malphas[i] = band_slope(alpha2d[ind], i, alpha_min, alpha_max)
        areas[i] = length(ind)*cellsz
        lengths[i] = dz/tan(malphas[i])
        widths[i] = areas[i]/lengths[i]
    end

    # Smooth the width:
    # Note, this can make the malphas noisy!
    if window_width_smooth>0
        widths = boxcar(widths, round(Int, window_width_smooth/bandsize) )
        for i=1:nb
            lengths[i] = areas[i]/widths[i]
            malphas[i] = atan(dz/lengths[i])
            @assert malphas[i]>=0
        end
    end

    x = vcat(0,cumsum(lengths))
    xmid = x[1:end-1] + diff(x)/2

    # tests
    if abs(totalarea-sum(areas))>1.0
        error("Something's amiss, sum of area of bands $(sum(areas)) not equal total area $totalarea,")
    end
    box_diag = sqrt((dem.x[end]-dem.x[1])^2 + (dem.y[end]-dem.y[1])^2)
    if abs(sum(lengths)-box_diag)/box_diag>0.4
        warn("Glacier length from might be wrong. Band-length: $(sum(lengths)/1e3)km, bounding-box diag: $(box_diag/1e3)km")
    end

    return bands, bandi, malphas, areas, lengths, widths, x, xmid, dem, alpha2d
end


"""
    band_slope(alphas, bandnr, alpha_min=deg2rad(0.4), alpha_max=deg2rad(60.0))

Need to calculate a meaningful mean of the slopes in a elevation bin.
*This is tricky but critical!*

One check can be that all the bin-lengths should add up to the total
glacier length.

Input:
- alphas -- slope angles in one band
- bandnr -- which band those alphas belong to (only used for error message)
- alpha_max, alpha_min -- maximal and minimal allowed slope
"""
function band_slope(alphas, bandnr, alpha_min, alpha_max)
    # parameters
    ratio_fac = 2
    f_q5 = 0.05
    f_q20 = 0.2
    f_q80 = 0.8
    q_band = (0.55, 0.95)

    n = length(alphas)
    if n==0
        error("Band with no element?!")
        return deg2rad(45)
    end
    # magic slope calculation
    sort!(alphas)
    # calculate indices of 5th, 20th and 80th quantiles:
    iq5  = max(1,round(Int,n*f_q5))
    iq20 = max(1,round(Int,n*f_q20))
    iq80 = min(n,round(Int,n*f_q80))
    # angle of those quantiles:
    q5, q20, q80 = map(rad2deg, (alphas[iq5], alphas[iq20], alphas[iq80]))
    # Now some of Matthias' magic:
    a = (q20/q80)*ratio_fac # 2x ratio of angles
    a = min(a, q_band[2])
    a = max(a, q_band[1]) # scaled to lie [0.55, 0.95]
    iq_magic = round(Int,n*a)  # set a "new" magic quantile to that value
    q_magic = rad2deg(alphas[iq_magic])
    # only use indices within those quantiles
    ind = q5 .<= rad2deg(alphas) .< q_magic
    if sum(ind)<2
        println("Magic band-slope calculation done with all cells for band $bandnr")
    end
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
function bin_grid(v::Matrix, binsize_or_bins, mask=BitArray([]); binround=-floor(Int, log10(binsize_or_bins)))
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
function bin_traj(tr::Traj, g::Gridded, binsize_or_bins, mask=trues(size(g.v)); binround=-floor(Int, log10(binsize_or_bins)))
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
    bandi_for_other_grid(binmat, g::Gridded, othergrid::Gridded, othermask=trues(size(othergrid.v))
    bandi_for_other_grid(bands, bandi, g::Gridded, othergrid::Gridded,
                         othermask=trues(size(othergrid.v)))

Returns vector of indices (bandi) to map a different grid onto the
bands encoded in `binmat` (or `bands, bandi`) and grid `g`.  It only
maps points onto bands which are within the mask applied to generate
the bands.  Additionally & optionally, a mask for the othergrid can
also be given.  Return:

    bandi
"""
function bandi_for_other_grid(bands, bandi::Vector{Vector{Int}}, g::Gridded, othergrid::Gridded,
                              othermask=trues(size(othergrid.v)))
    binmat=bins2matrix(g, bands, bandi)
    bandi_for_other_grid(bands, binmat, g, othergrid, othermask)
end
function bandi_for_other_grid(bands, binmat::Matrix{Int}, g::Gridded, othergrid::Gridded, othermask=trues(size(othergrid.v)) )
    og = othergrid
    @assert size(othergrid)==size(othermask)
    if g.x!=og.x || g.y!=og.y
        bandi_ = [Int[] for i=1:length(bands)]
        dims = size(og.v)
        itpm = Interpolations.interpolate((g.x, g.y), binmat, Interpolations.Gridded(Interpolations.Constant()) );
        itpm = Interpolations.extrapolate(itpm, 0);
        for j=1:size(og.v,2)
            for i=1:size(og.v,1)
                if othermask[i,j]
                    ind = itpm[og.x[i], og.y[j]]
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

"""
One cell-edge of a regular gird of cells.  Arguably not the best
datastructure for what is done below.
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
    orientation(e1::Edge, e2::Edge)

Return whether the cells are on left, right or both of the two cells.

"""
function orientation(e1::Edge, e2::Edge)
    n1,n2 = get_nodes(e1)
    m1,m2 = get_nodes(e2)
    if n2==m1
        return right
    elseif n1==m2
        return left
    elseif n1==m1 || n2==m2
        return _noedge
    else
        error("The two edges:\n $e1 $e2 are not connected")
    end
end

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
orientation(l::Line) = orientation(l.edges)
function orientation(l::Vector)
    length(l)==1 && return _noedge # TODO what to do here?
    # establish orientation of first two segments
    ori = orientation(l[1],l[2])
    for i=3:length(l)
        orientation(l[i-1],l[i])!=ori && error("Line not consistently oriented")
    end
    return ori
end

"""
    next_edge!(edges::Set{Edge}, edge::Edge)

Return an adjacent edge of input edge.  Pops returned edge off
`edges`.  If no edge is found, return `Edge(-1,-1,nodege)`.
"""
function next_edge!(edges::Set{Edge}, edge::Edge)
    edge.loc==_noedge && error()
    # lookup table for where to look for next edge.
    tolookat = Dict(
                left =>Dict(
                        (-1,1) => (lower,right),
                        (-1,0) => (lower,upper),
                        (-1,-1)=> (upper,right),
                        (0,1)  => (left,lower),
                        (0,0)  => (upper,lower),
                        (0,-1) => (left,upper)),
                right=>Dict(
                        (0,1)  => (lower,right),
                        (0,0)  => (lower,upper),
                        (0,-1) => (upper,right),
                        (1,1)  => (left,lower),
                        (1,0)  => (upper,lower),
                        (1,-1) => (left,upper)),
                lower=>Dict(
                        (-1,0) => (lower,right),
                        (0,0)  => (left,right),
                        (1,0)  => (left,lower),
                        (-1,-1)=> (upper,right),
                        (0,-1) => (left,right),
                        (1,-1) => (left,upper)),
                upper=>Dict(
                        (-1,1) => (lower,right),
                        (0,1)  => (left,right),
                        (1,1)  => (left,lower),
                        (-1,0) => (upper,right),
                        (0,0)  => (left,right),
                        (1,0)  => (left,upper))
                )
    look = tolookat[edge.loc]
    for ((i,j), locs) in look
        for loc in locs
            test_edge = Edge(edge.i+i,edge.j+j,loc)
            if test_edge in edges
                return pop!(edges, test_edge)
            end
        end
    end
    # nothing found
    return noedge
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
        binmat[landmask.==0 & binmat.==0] = -1
    end

    # The boundaries of all bands:
    # (bandnr,otherbandnr) => Set of edges
    boundaries_ = Dict{Tuple{Int,Int},Set{Edge}}()

    # Loop to identify all boundary cells of each band and determine
    # which cell-edges are on the boundary:
    for (ib,bup) in enumerate(bands)
        for I in bandi[ib]
            i,j = ind2sub(dims, I)
            loc = 1
            # left-right
            for ii=-1:2:1
                iother = binmat[i+ii,j]
                if iother!=ib
                    push!(get!(boundaries_, (ib, iother), Set{Edge}()),
                          Edge(i,j,loc))
                end
                loc+=1
            end
            # lower-upper
            for jj=-1:2:1
                iother = binmat[i,j+jj]
                if iother!=ib
                    push!(get!(boundaries_, (ib, iother), Set{Edge}()),
                          Edge(i,j,loc))
                end
                loc+=1
            end
        end
    end
    # Sort the points and find breaks.

    # bandnr => otherbandnr => Vector{Line}
    boundaries = [Dict{Int,Vector{Line}}() for i=1:length(bands)]

    for ((band,otherband), edges) in boundaries_

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
            for ii=1:2 # both directions
                curedge = firstedge
                while curedge!=noedge
                    push!(out[ii], curedge)
                    curedge = next_edge!(edges, curedge)
                end
            end
            prepend!(out[2], reverse(out[1][2:end]))
            # make sure the elevation band is on the line's left:
            if orientation(out[2]) == right
                push!(output, Line(reverse(out[2])))
            else
                push!(output, Line(out[2]))
            end
        end
        boundaries[band][otherband] = output
    end

    return boundaries
end

"""
    calc_fluxdir_source(l::Line, ux, uy, window) -> fluxdir,ii,jj,fluxdirx,fluxdiry

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
function calc_fluxdir(l::Line, ux, uy, window, dx)
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
                fx<0 ? (i,j,-fx*dx) : (i-1,j,-fx*dx)
            elseif loc==right
                fx<0 ? (i+1,j,fx*dx) : (i,j,fx*dx)
            elseif loc==lower
                fy<0 ? (i,j,-fy*dx) : (i,j-1,-fy*dx)
            else #if loc==upper
                fy<0 ? (i,j+1,fy*dx) : (i,j,fy*dx)
            end
    end
    return fluxdir,ii,jj,fluxdirx,fluxdiry
end

"""
    calc_u(q1d, boundaries, thick, alpha, ux, uy, dx, mask, bands, lengths,
           flux_dir_window, window_frac=1.5,
           x=nothing, y=nothing) # these are only needed for plotting
    -> u2d, u2d_at_bands, scaling-factors-1d, mask_u2d_at_bands

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
                                    ux, uy, dx, mask, bands, lengths,
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
    ubar_, ubar, facs, mask_ubar = _calc_u(q1d, boundaries, u_trial, thick,
                                    ux, uy, dx, mask, bands, lengths,
                                    flux_dir_window, # in [m]
                                    plotyes,x,y)
    return boxcar(ubar_, Int((window_frac*maximum(lengths))÷dx)+1, mask_ubar, !mask), ubar, facs, mask_ubar
end

# this helper function is needed for type stability
function _calc_u(q1d, boundaries, u_trial, thick,
                ux, uy, dx, mask, bands, lengths,
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
            bb = bnd[ib+1]
        else # outflow at terminus
            # TODO: this probably needs updating where several elevation bands contribute (tide-water)
            # -> no only ever the last band does outflow
            bb = get(bnd, -1, bnd[0]) # if sea-terminating use those edges.
        end
        # loop over segments
        ffs = Float64[]
        is = Int[]
        js = Int[]
        for l in bb
            #@assert orientation(l)
            ff, fi, fj, fx, fy = calc_fluxdir(l, ux, uy, Int(flux_dir_window÷dx), dx)
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
                # This should not happen!
                error("This should not happen!")
                u[n]=0
            end
            ubar[i,j] = u[n]
        end
    end
    mask_ubar = mask & (!isnan(ubar)) # location of all cells for which `ubar` was calculated
    ubar_ = copy(ubar)
    ubar_[isnan(ubar_)] = 0
    # if plotyes
    #     # imshow(binmat',origin="lower", extent=(x[1],x[end],y[1],y[end]), cmap="flag"); colorbar();
    #     imshow(ubar',origin="lower", extent=(x[1],x[end],y[1],y[end]),); colorbar(); clim(0,50)
    # end

    return ubar_, ubar, facs, mask_ubar
end
