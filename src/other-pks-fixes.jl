# Fixes for other packages or Julia-Base

# Make RasterIO conditional as it is not Julia 0.6 compatible
if haskey(Pkg.installed(), "RasterIO")
    eval(quote
         import RasterIO
         # https://github.com/wkearn/RasterIO.jl/pull/26
         function getrasternodatavalue(rasterband::RasterIO.GDALRasterBandH)
         success = Cint[0]
         nodata = RasterIO._getrasternodatavalue(rasterband, pointer(success))
         nodata, Bool(success[])
         end
         end)
end
