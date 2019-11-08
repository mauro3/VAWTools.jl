#########
# Binning / elevation-bands
#########
x = -100.3:0.56:-12
y = 1001:2.3:1100
ele = x .+ y'
mask = trues(size(ele))
mask[1:10,:] = false
mask[[1,end],:] = false
mask[:,[1,end]] = false
g = Gridded(x,y,ele)
# negative step
bands, bandi = bin_grid(g, -10.0, mask)
@test length(bands)==length(bandi)
binmat = VAWTools.bins2matrix(g, bands, bandi)
@test all(binmat[1:10,:].==0)
@test all(binmat[:,1].==0)
# positive step
bands, bandi = bin_grid(g, 10.0, mask)
@test length(bands)==length(bandi)
binmat = VAWTools.bins2matrix(g, bands, bandi)
@test all(binmat[1:10,:].==0)
@test all(binmat[:,1].==0)
# adjusting step
bands, bandi = bin_grid(g, 200.0, mask)
@test length(bands)==length(bandi)
binmat = VAWTools.bins2matrix(g, bands, bandi)
@test all(binmat[1:10,:].==0)
@test all(binmat[:,1].==0)
# range steps:
bands_, bandi_ = bin_grid(g, bands, mask)
@test all(bands_.==bands)
@test all(bandi_.==bandi)
# non range steps:
@test_throws MethodError bands_, bandi_ = bin_grid(g, collect(bands), mask)
# @test all(bands_.==bands)
# @test all(bandi_.==bandi)


bands, bandi = bin_grid(g, 10.0, mask)
binmat = VAWTools.bins2matrix(g, bands, bandi)
xx,yy = -50:11:23, 900.0:7:1200
othergrid = Gridded(xx, yy, rand(length(xx),length(yy)))
bandii = VAWTools.bandi_for_other_grid(bands, bandi, binmat, g, othergrid)
bandii2 = VAWTools.bandi_for_other_grid(bands, bandi, g, othergrid)
@test bandii==bandii2

binmatt = VAWTools.bins2matrix(othergrid, bands, bandii)
@test binmatt == [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  6   7   7   8   9  10  10  11  12  12  13  14  14  15  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  7   8   9   9  10  11  11  12  13  13  14  15  15  16  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  8   9  10  10  11  12  12  13  14  14  15  16  17  17  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  9  10  11  11  12  13  14  14  15  16  16  17  18  18  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]

# 1D glacier
bands_g, bandi_g, malphas, areas, lengths, widths, x, xmid, dem, alpha2d =
    make_1Dglacier(g, 10.0, mask)
@test bands_g==bands
@test bandi_g==bandi
bands_g, bandi_g, malphas, areas, lengths, widths, x, xmid, dem, alpha2d =
    make_1Dglacier(g, bands_g, mask)
@test all(bands_g.==bands)
@test all(bandi_g.==bandi)

###############
# Fluxes
###############
tmp       = [14 14 15 15 16 17 18 18 19 0 0 0;
             14 14 15 15 16 17 18 19 19 0 0 0;
             14 15 15 16 16 17 18 19 19 0 0 0;
             14 15 15 16 17 18 18 19 0 0 0 0;
             15 15 16 16 17 18 19 19 0 0 0 0;
             15 16 16 17 17 18 19 20 0 0 0 0;
             15 16 16 17 18 17 20 21 0 0 0 0;
             16 16 17 18 18 17 21 22 23 0 0 0;
             16 17 17 18 19 20 21 23 24 0 0 0;
             16 17 18 19 19 20 22 23 0 0 24 24;
             17 18 19 19 20 21 22 0 23 24 24 24]-13;
elevation = zeros(size(tmp,1)+2,size(tmp,2)+2)-13
elevation[2:end-1,2:end-1] = tmp
mask = elevation.>0
bands, bandi = VAWTools.bin_grid(elevation,1,mask)
binmat = VAWTools.bins2matrix(elevation,bands,bandi)
# using PyPlot
# bm = binmat*1.0
# bm[bm.==0] = NaN
# matshow(bm',cmap="flag", origin="lower", extent=(1,size(bm,1)+1,1,size(bm,2)+1)); colorbar()
# xticks(1:size(bm,1))
# yticks(1:size(bm,2))
# grid()
@test all(elevation[mask].==binmat[mask])
@test all(elevation[(!).(mask)].==-13)
@test all(binmat[(!).(mask)].==0)
edges_at_boundaries = VAWTools.get_cells_on_boundary(bands, bandi, binmat)
@test sort(collect(keys(edges_at_boundaries)))== Tuple{Int64,Int64}[(1,0),(1,2),(2,0),(2,1),(2,3),(3,0),(3,2),(3,4),(4,0),(4,3),(4,5),(4,7),(4,8),(5,0),(5,4),(5,6),(6,0),(6,5),(6,7),(7,0),(7,4),(7,6),(7,8),(7,9),(8,0),(8,4),(8,7),(8,9),(8,10),(9,0),(9,7),(9,8),(9,10),(10,0),(10,8),(10,9),(10,11),(11,0),(11,10)]

boundaries = VAWTools.calc_boundaries(bands, bandi, binmat)
# symmetry
for i=1:length(bands), j=1:length(bands)
    if haskey(boundaries[j],i)
        @test length(boundaries[i][j])==length(boundaries[j][i])
    end
end
@test length(boundaries[2-1][3-1])==1
@test length(boundaries[3-1][2-1])==1
@test length(boundaries[3-1][4-1])==1
@test length(boundaries[4-1][5-1])==1
@test length(boundaries[4][3])==1
@test length(boundaries[4][5])==2
@test length(boundaries[4][7])==2
@test length(boundaries[4][8])==1
@test length(boundaries[5][4])==2
@test length(boundaries[5][6])==2
@test length(boundaries[6][5])==2
@test length(boundaries[8-1][7-1])==2
@test length(boundaries[8-1][9-1])==3
@test length(boundaries[9-1][8-1])==3
@test length(boundaries[9-1][10-1])==3
@test length(boundaries[9-1][11-1])==1
@test length(boundaries[10-1][9-1])==3
@test length(boundaries[10-1][11-1])==2
@test length(boundaries[10-1][1-1])==1
@test length(boundaries[11-1][1-1])==3
@test length(boundaries[11-1][10-1])==2
@test length(boundaries[11-1][12-1])==2
@test length(boundaries[12-1][11-1])==2
@test length(boundaries[12-1][1-1])==2

q1d = bands+1
u_trial = ones(elevation)
thick = sqrt.(elevation+13)
ux,uy = (-).(gradient3by3(1:13,1:14,elevation,mask))
dx = 1.0
lengths = rand(length(bands))
flux_dir_window = 4
window_frac = 0.4
plotyes=false
x=nothing
y=nothing
ubar, mask_ubar, facs = VAWTools._calc_u(q1d, boundaries, u_trial, thick,
                                                ux, uy, dx, mask, bands,
                                                flux_dir_window, # in [m]
                                                plotyes,x,y)


u2d, u2d_at_bands, scaling_factors_1d, mask_u2d_at_bands = VAWTools.calc_u(q1d, boundaries, u_trial, thick,
                                                                           ux, uy, dx, mask, bands, lengths,
                                                                           flux_dir_window, # in [m]
                                                                           window_frac)

# more tests would be good...
@test isequal(u2d_at_bands, ubar)

### Round trip 1D -> 2D -> 1D
#
# Note that 2D -> 1D -> 2D will not round-trip as information
# is lost in the first ->

x = -100.3:0.56:-12
y = 1001:2.3:1100
ele = 0.01*x.^2 .+ y';
ele[:, end÷2:end] += 120; # make a step, just to make it awkward
# contour(x,y,ele')
dem = Gridded(x,y,ele)
mask = sqrt.((x.+60).^2 .+ (y'.-1050).^2) .<45 #; ele[.!mask] = NaN; heatmap(ele)
maxs = [maximum(ele),maximum(ele[mask])]
mins = [minimum(ele),minimum(ele[mask])]

bands = reverse(mean(mins):5:maxs[1]+20)
vals1d = 0.4 * [149, 957, -777, -732, -260, 168, 491, 424, -709, 505, -52, 252, -6, 388, -98, -378, -568, -363, -885, 705, 761, 485, -346, -803, 950, -245, -253, 471, -303, 567, -108,
 773, -545, 23, -57, 524, -730, 521, -232, 600, 311, -862, -284, 151, 694, 541, 565, 707, 184, -189, 192, -471, -283, 293, -669, -244, -808, 155, -553, -392, -769, -940,
 828, -489, 863, 126]
# test that it's good input data
@test length(bands) == length(vals1d)
@test length(unique(vals1d))==length(vals1d)

bands_, bandi = VAWTools.bin_grid(dem, bands, mask, min_bin_number_ends=5);
@test length(bands_)==53
bands_, bandi = VAWTools.bin_grid(dem, bands, mask);
@test bands==bands_

# check that mask works:
@test !any(mask[vcat(bandi...)].==false)
mm = trues(mask)
mm[vcat(bandi...)] = false
@test !any(mask[mm])

vals2d = VAWTools.map_back_to_2D(size(dem.v), bandi, vals1d)

# test that all values are the same in each band
for bi in bandi
    for b in bi
        @test vals2d[bi[1]]==vals2d[b]
    end
end

vals1d_ = VAWTools.map_onto_bands(bandi, vals2d)
@test all(vals1d_[length.(bandi).>0] .≈ vals1d[length.(bandi).>0])
@test all(isnan.(vals1d_[length.(bandi).==0]))
