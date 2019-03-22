using VAWTools
using Base.Test
using SHA, Compat

function sha(fn)
    open(fn) do f
        sha2_256(f)
    end
end

include("smoothing-functions.jl")
include("elevation-bands.jl")

## Gridded
##########
g = Gridded(1.:37, 1:17., rand(37,17))
gc = Gridded(1.:10,1:11., 5*ones(10,11) )

@test g.v==downsample(g,1).v
@test_throws ArgumentError downsample(g,0)

for i=1:9
    @test all(downsample(gc,i).v.==5)
    @test downsample(g,i).v==boxcar(g.v,Int(floor((i)/2)))[1:i:end,1:i:end]
end
@test_throws AssertionError downsample(gc,10)

gc.v[4,5] = NaN
for i=1:9
    ds = downsample(gc,i).v
    @test all((ds.==5) .| (isnan.(ds)))
    @test downsample(g,i).v==boxcar(g.v,Int(floor((i)/2)))[1:i:end,1:i:end]
end
averagemask = trues(size(gc.v))
for i=1:9
    ds = downsample(gc,i,1,true,averagemask).v
    @test all((ds.==5) .| (isnan.(ds)))
    @test downsample(g,i,1,true,trues(size(g.v))).v==boxcar(g.v,Int(floor((i)/2)))[1:i:end,1:i:end]
end
averagemask = trues(size(g.v))
averagemask[2,1] = false
for i=1:9
    ds = downsample(g,i,1,true,averagemask).v
    bc = boxcar(g.v,Int(floor((i)/2)),averagemask)[1:i:end,1:i:end]
    @test bc==ds
end



## File IO
#################
const todelete = []
function tempfn(ext="")
    out = tempname()*ext
    push!(todelete, out)
    out
end
function delfls()
    for f in todelete
        if isfile(f)
            rm(f)
        end
    end
end

## .agr readers:
g1 = VAWTools._read_agr("testfiles/wiki.agr")
g1b = VAWTools._read_agr("testfiles/wiki.bin")
@test g1.v == [-9999 -9999 5 2
                -9999 20 100 36
                3 8 35 10
                32 42 50 6
                88 75 27 9
                13 5 1 -9999]
@test g1.nc==4
@test g1.nr==6
@test g1.xll==0.0
@test g1.yll==-40.0
@test g1.dx==50.0
@test g1.NA==-9999
# write and read back:
tfn1 = tempfn(".agr")
write_agr(g1, tfn1)
g2 = VAWTools._read_agr(tfn1, Float64)
@test g1==g2
tfn2 = tempfn(".agr")
write_agr(g2, tfn2)
@test sha(tfn1)==sha(tfn2)

# reading from io thingies
io = open("testfiles/wiki.agr", "r")
g2 = VAWTools.read_agr(io)
close(io)
g1 = VAWTools.read_agr("testfiles/wiki.agr")
@test g1.x==g2.x
@test g1.y==g2.y
@test g1.midpoint==g2.midpoint

# Matthias abuse of NODATA_value field as UTM zone:
VAWTools.read_agr("testfiles/wiki.agr")


## RasterIO agr reading
g1 = VAWTools.read_agr("testfiles/wiki.agr")
# println("Expected output: ERROR 1: ERROR - failed to load SRS ...")
# g2,proj4 = VAWTools.read_rasterio("testfiles/wiki.agr")
# @test isequal(g1.v,g2.v)
# @test g1.x==g2.x
# @test g1.y==g2.y
# @test g1.midpoint==g2.midpoint
# @test proj4==""

# write and read back binary:
tfn3 = tempfn(".bin")
write_agr(g1, tfn3, NA_agr=-9999.0)
g2 = VAWTools._read_agr(tfn3)
g1_ = VAWTools._read_agr("testfiles/wiki.agr", Float32)
@test g1_==g2
@test sha(tfn3)==sha("testfiles/wiki.bin")

# larger file:
g3 = VAWTools._read_agr("testfiles/t2.bin", Float32)
@test g3.nc==121
@test g3.nr==81
@test g3.xll==678750.0f0
@test g3.yll==140250.0f0
@test g3.dx==25.0f0
@test g3.NA==-9999.0f0

tfn4 = tempfn(".bin")
write_agr(g3, tfn4)
@test sha("testfiles/t2.bin")==sha(tfn4)

# Gridded
g1 = VAWTools._read_agr("testfiles/wiki.agr")
gg = Gridded(g1)
@test isequal(VAWTools.AGR(gg, NA_agr=-9999).v, g1.v)

# geotiff
gt = VAWTools.read_geotiff("testfiles/small_world.tif")
gt_ag = VAWTools.read_geotiff("testfiles/small_world.asc")
@test gt.proj=="+proj=longlat +datum=WGS84 +no_defs"
@test gt.x==gt_ag.x
@test gt.y==gt_ag.y
@test isa(gt.v, Array{Float32,2})
@test gt.v==gt_ag.v

# Trajectory
@test_throws AssertionError Traj(1:5, 6:11, [0.0:5.0;])

tr = Traj(1.0:5, 6.0:10)
@test tr.x==collect(1.0:5)
@test tr.y==collect(6.0:10)
@test tr.splits==[1:5]
@test !VAWTools.haserror(tr)
@test !VAWTools.hasvalues(tr)

tr = Traj(1.0:5, 6.0:10, Void[], Void[], [1:3,4:5])
@test tr.x==collect(1.0:5)
@test tr.y==collect(6.0:10)
@test tr.splits==[1:3,4:5]
@test !VAWTools.haserror(tr)
@test !VAWTools.hasvalues(tr)

tr = Traj(1.0:5, 6.0:10, 1.0:5.0)
@test tr.x==collect(1.0:5)
@test tr.y==collect(6.0:10)
@test tr.v==collect(1.0:5)
@test tr.splits==[1:5]
@test !VAWTools.haserror(tr)
@test VAWTools.hasvalues(tr)

tr = Traj(1.0:5, 6.0:10, 1.0:5.0, Float64[], [1:3,4:5])
@test tr.x==collect(1.0:5)
@test tr.y==collect(6.0:10)
@test tr.v==collect(1.0:5)
@test tr.splits==[1:3,4:5]
@test !VAWTools.haserror(tr)
@test VAWTools.hasvalues(tr)

tr = Traj(1.0:5, 6.0:10, 1.0:5.0, 1.0:5.0)
@test tr.x==collect(1.0:5)
@test tr.y==collect(6.0:10)
@test tr.v==collect(1.0:5)
@test tr.err==collect(1.0:5)
@test tr.splits==[1:5]
@test VAWTools.haserror(tr)
@test VAWTools.hasvalues(tr)

tr = Traj(1.0:5, 6.0:10, 1.0:5.0, 1.0:5.0, [1:3,4:5])
@test tr.x==collect(1.0:5)
@test tr.y==collect(6.0:10)
@test tr.v==collect(1.0:5)
@test tr.err==collect(1.0:5)
@test tr.splits==[1:3,4:5]
@test VAWTools.haserror(tr)
@test VAWTools.hasvalues(tr)
@test all(VAWTools.distances(tr).==sqrt(2))

# split it
tr = Traj([1,2,4,5], 1:4)
@test tr.splits==[1:4]
@test !VAWTools.haserror(tr)
@test !VAWTools.hasvalues(tr)
VAWTools.split_traj!(tr,2)
@test tr.splits==[1:2,3:4]
@test !VAWTools.haserror(tr)
@test !VAWTools.hasvalues(tr)


## .xyn reader
po = read_xyn("testfiles/poly.xyn")
@test length(po)==1
@test size(po[1])==(2,24)
cpo = VAWTools.concat_poly(po)
@test po[1]==cpo[1]
@test po==VAWTools.split_poly(cpo...)


po = read_xyn("testfiles/multipoly.xyzn", hasz=true)
@test length(po)==4
@test size(po[1])==(3,4)
@test size(po[2])==(3,5)
@test size(po[3])==(3,4)
@test size(po[4])==(3,3)
cpo = VAWTools.concat_poly(po)
@test po==VAWTools.split_poly(cpo...)


# tidy up temp-files
delfls()


## IN poly
@test VAWTools.leftorright(0.5,0.5, 1,0,1,1)==-1
@test VAWTools.leftorright(1.5,.5, 1,0,1,1)==1
@test VAWTools.leftorright(1,0.5, 1,0,1,1)==0

poc = VAWTools.concat_poly(po)
@test !inpoly([0.0,0], poc[1][1:2,:])
@test inpoly([6.28793e5, 90129.7+10], poc[1][1:2,:])

# counterclockwise poly
poly = Float64[0  0  1  1  0
               0  1  1  0  0]
# plot with:
# using Plots; plot(poly[1,:], poly[2,:])

p1 = [0.5, 0.5]
p2 = [0.5, 0.99]
p22 = [0.5, 1] # on top edge (per def outside)
p22_ = [1, 0.5] # on right edge (per def outside)
p23 = [0.5, 0] # on bottom edge (per def inside)
p23_ = [0, 0.5] # on left edge (per def inside)
p24 = [0, 0]   # on corner
p25 = [0, .4]   # on edge
p3 = [0.5, 1.1]

@test inpoly(p1, poly)
@test inpoly(p2, poly)
@test !inpoly(p22, poly)
@test !inpoly(p22_, poly)
@test inpoly(p23, poly)
@test inpoly(p23_, poly)
@test inpoly(p24, poly)
@test inpoly(p25, poly)
@test !inpoly(p3, poly)

# clockwise poly
poly = Float64[0  1  1  0  0
               0  0  1  1  0]

@test inpoly(p1, poly)
@test inpoly(p2, poly)
@test !inpoly(p22, poly)
@test !inpoly(p22_, poly)
@test inpoly(p23, poly)
@test inpoly(p23_, poly)
@test inpoly(p24, poly)
@test inpoly(p25, poly)
@test !inpoly(p3, poly)


# cross-over poly
poly = Float64[ 0  1  0  1  0
                0  0  1  1  0]

if VERSION>=v"0.5-"
    eval(:(@test_broken inpoly(p1, poly) )) # should be true
end
@test inpoly(p2, poly)
@test !inpoly(p22, poly)
@test !inpoly(p22_, poly)
@test inpoly(p23, poly)
@test !inpoly(p23_, poly)
@test inpoly(p24, poly)
@test !inpoly(p25, poly) # different
@test !inpoly(p3, poly)


# with interior region
poly = Float64[0 0
               # interior
               0.1 0.1
               0.1 0.6
               0.6 0.6
               0.6 0.1
               0.1 0.1
               # exterior
               0 0
               0 1
               1 1
               1 0
               0 0]'
# inside interior poly: i.e. labeled as outside
@test !inpoly([0.3,0.3], poly)
@test !inpoly([0.3,0.5], poly)

poly = Float64[0 0
               # interior
               0.1 0.1
               0.1 0.6
               0.6 0.6
               # in-interior
               0.4 0.4
               0.4 0.2
               0.2 0.2
               0.2 0.4
               0.4 0.4
               # interior
               0.6 0.6
               0.6 0.1
               0.1 0.1
               # exterior
               0 0
               0 1
               1 1
               1 0
               0 0]'
# inside in-interior poly
@test inpoly([0.3,0.3], poly)
@test !inpoly([0.3,0.5], poly)

poly = Float64[0 0
               # interior
               0.1 0.1
               0.1 0.6
               0.6 0.6
               # in-interior
               0.4 0.4
               0.2 0.4
               0.2 0.2
               0.4 0.2
               0.4 0.4
               # interior
               0.6 0.6
               0.6 0.1
               0.1 0.1
               # exterior
               0 0
               0 1
               1 1
               1 0
               0 0]'
# inside in-interior poly
@test inpoly([0.3,0.3], poly)
@test !inpoly([0.3,0.5], poly)

poly = Float64[0 0
               # interior #1
               0.1 0.1
               0.1 0.6
               0.4 0.6
               0.4 0.6
               0.4 0.1
               0.1 0.1
               0 0
               # interior #2
               0.6 0.4
               0.6 0.6
               0.8 0.6
               0.8 0.4
               0.6 0.4
               0 0
               # exterior
               0 1
               1 1
               1 0
               0 0]'
@test !inpoly([0.2,0.4], poly)
@test !inpoly([0.3,0.15], poly)
@test inpoly([0.5,0.4], poly)
@test inpoly([0.5,0.2], poly)
@test !inpoly([0.7,0.5], poly)

# Test points which are on same y-level as a endpoint of poly edge
poly = Float64[0.1 1 2 1  0.1
               0.1 -1 0 1 0.1];
@test inpoly([1,0.0], poly)
@test inpoly([1,-0.0001], poly)
@test inpoly([1,0.0001], poly)
poly = Float64[0.0 1 2 1  0.0
               0.0 -1 0 1 0.0];
@test inpoly([1,0.0], poly)
@test inpoly([1,-0.0001], poly)
@test inpoly([1,0.0001], poly)


#####
# boxcar
#####
## 1D
orig = [1:10.0;]
window = 3
@inferred VAWTools.boxcar(orig, window)
filt1 = VAWTools.boxcar(orig, window)
@test size(filt1)==size(orig)
@test filt1==[2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 7.5, 8.0, 8.5]

window = 3.3
@inferred VAWTools.boxcar(orig, window)
filt1 = VAWTools.boxcar(orig, window)
@test size(filt1)==size(orig)
@test filt1≈[2.67442, 3.16981, 3.66667, 4.16438, 5.0, 6.0, 6.83562, 7.33333, 7.83019, 8.32558] atol=1e-4

windows = (2.4,3.3)
@inferred VAWTools.boxcar(orig, windows)
filt1 = VAWTools.boxcar(orig, windows)
@test size(filt1)==size(orig)
@test filt1≈[2.67442, 3.16981, 3.66667, 4.44776, 5.44776, 6.44776, 7.28125, 7.77778, 8.27273, 8.76471] atol=1e-4


## 2D
#nr,nc = 4,5
T = Float32
nr,nc = 20,31
#nr,nc = 600,500
window = 3
# orig = (1.0:nr)''*(1:nc)'
# weights = ones(Int,nr,nc)
srand(1)
orig = rand(T,nr,nc)
weights = rand(nr,nc)
weightsb = bitrand(nr,nc)
weightsbb = convert(Matrix{Bool}, weightsb)

dropmask = falses(nr,nc)
dropmask[5,6] = true
# with weights
filt1 = VAWTools.boxcar(orig, window, weights)
@inferred VAWTools.boxcar(orig, window, weights)
M = VAWTools.boxcar_matrix(T, window, weights)
@inferred VAWTools.boxcar_matrix(T, window, weights)
@test size(M,1)==size(M,2)
@test size(M,1)==length(orig)
filt2 = VAWTools.apply_boxcar_matrix(M, orig)
@test eltype(orig)==eltype(filt1)
@test eltype(orig)==eltype(filt2)
for i=eachindex(orig)
    @test filt1[i] ≈ filt2[i]
end
# with weightsb
filt1 = VAWTools.boxcar(orig, window, weightsb)
@inferred VAWTools.boxcar(orig, window, weightsb)
M = VAWTools.boxcar_matrix(T, window, weightsb)
@inferred VAWTools.boxcar_matrix(T, window, weightsb)
@test size(M,1)==size(M,2)
@test size(M,1)==length(orig)
filt2 = VAWTools.apply_boxcar_matrix(M, orig)
@test eltype(orig)==eltype(filt1)
@test eltype(orig)==eltype(filt2)
for i=eachindex(orig)
    @test filt1[i] ≈ filt2[i]
end
# with weightsbb
filt1 = VAWTools.boxcar(orig, window, weightsbb)
@inferred VAWTools.boxcar(orig, window, weightsbb)
M = VAWTools.boxcar_matrix(T, window, weightsbb)
@inferred VAWTools.boxcar_matrix(T, window, weightsbb)
@test size(M,1)==size(M,2)
@test size(M,1)==length(orig)
filt2 = VAWTools.apply_boxcar_matrix(M, orig)
@test eltype(orig)==eltype(filt1)
@test eltype(orig)==eltype(filt2)
for i=eachindex(orig)
    @test filt1[i] ≈ filt2[i]
end

# with weights and dropmask
filt1 = VAWTools.boxcar(orig, window, weightsbb, dropmask)
@inferred VAWTools.boxcar(orig, window, weightsbb, dropmask)
M = VAWTools.boxcar_matrix(T, window, weightsbb, dropmask)
@inferred VAWTools.boxcar_matrix(T, window, weightsbb, dropmask)
@test size(M,1)==size(M,2)
@test size(M,1)==length(orig)
filt2 = VAWTools.apply_boxcar_matrix(M, orig)
@test eltype(orig)==eltype(filt1)
@test eltype(orig)==eltype(filt2)
for i=eachindex(orig)
    @test filt1[i] ≈ filt2[i]
end

# with NaN poisoning
orig[1,1] = NaN
weightsbb[1,1] = 1
filt1 = VAWTools.boxcar(orig, window, weightsbb, dropmask)
M = VAWTools.boxcar_matrix(T, window, weightsbb, dropmask)
@test size(M,1)==size(M,2)
@test size(M,1)==length(orig)
filt2 = VAWTools.apply_boxcar_matrix(M, orig)
@test eltype(orig)==eltype(filt1)
@test eltype(orig)==eltype(filt2)
for i=eachindex(orig)
    if isnan(filt1[i])
        @test isequal(filt1[i], filt2[i])
    else
        @test filt1[i] ≈ filt2[i]
    end
end

# NaN is masked: no poisoning
orig[1,1] = NaN
weightsbb[1,1] = 0
filt1 = VAWTools.boxcar(orig, window, weightsbb, dropmask)
M = VAWTools.boxcar_matrix(T, window, weightsbb, dropmask)
@test size(M,1)==size(M,2)
@test size(M,1)==length(orig)
filt2 = VAWTools.apply_boxcar_matrix(M, orig)
@test eltype(orig)==eltype(filt1)
@test eltype(orig)==eltype(filt2)
for i=eachindex(orig)
    if isnan(filt1[i])
        @test isequal(filt1[i], filt2[i])
    else
        @test filt1[i] ≈ filt2[i]
    end
end

## different left and right windows
windows1 = (3,4)
windows2 = ((2,3), (1,2))

srand(1)
orig = rand(T,nr,nc)

for w in [windows1, windows2]
    @inferred VAWTools.boxcar(orig, w)
    @inferred VAWTools.boxcar(orig, w, weights)
    @inferred VAWTools.boxcar(orig, w, weights, dropmask)
end

@inferred VAWTools.boxcar(orig, (3.3,4.2))

####
# piecewiselinear
####
f = VAWTools.piecewiselinear([78.0], [-12.0])
@test f(-1)==-12
@test f(100)==-12

f = VAWTools.piecewiselinear([0, 1, 2], [-1, 1, -2])
@test f(0)==-1
@test f(0.5)==0
@test f(1)==1
@test f(1.5)==-0.5
@test f(2)==-2


#####
# projections
####
import Proj4
wgs84 = Proj4.Projection("+proj=longlat +datum=WGS84 +no_defs")
utm56 = Proj4.Projection("+proj=utm +zone=56 +south +datum=WGS84 +units=m +no_defs")

pt  = [150.0 -27 110]
@test isapprox(VAWTools.transform_proj(VAWTools.transform_proj(pt, wgs84, utm56), utm56, wgs84), pt)
