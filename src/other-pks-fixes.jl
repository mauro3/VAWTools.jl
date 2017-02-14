# Fixes for other packages or Julia-Base

import RasterIO

# https://github.com/wkearn/RasterIO.jl/pull/26
eval(RasterIO, :(
    function getrasternodatavalue(rasterband::GDALRasterBandH)
    success = Ref{Cint}()
    nodata = _getrasternodatavalue(rasterband, success)
    nodata, Bool(success[])
    end))
