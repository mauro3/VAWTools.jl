##############
# Bands
##############
#
# TODO: add extrapolation too (second half of make_bands)

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
    bins2matrix(g::Union{Gridded,Matrix}, bands, bandi) = bins2matrix(g.v, bands, bandi)

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
