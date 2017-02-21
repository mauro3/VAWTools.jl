# Fixes for other packages or Julia-Base

import RasterIO

# https://github.com/wkearn/RasterIO.jl/pull/26
eval(RasterIO, :(
    function getrasternodatavalue(rasterband::GDALRasterBandH)
    success = Cint[0]
    nodata = _getrasternodatavalue(rasterband, pointer(success))
    nodata, Bool(success[])
    end))
