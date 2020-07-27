# ####################################### LOAD REQUIRED LIBRARIES ############################################# #
import time
import gdal, osr, ogr
import numpy as np
# ####################################### SET TIME-COUNT ###################################################### #
starttime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print("--------------------------------------------------------")
print("Starting process, time: " +  starttime)
print("")
# ####################################### DRIVERS AND FOLDER-PATHS ############################################ #
drvMemR = gdal.GetDriverByName('MEM')
drvR = gdal.GetDriverByName('GTiff')
drvMemV = ogr.GetDriverByName('Memory')
workFolder = "D:/OneDrive - Conservation Biogeography Lab/_RESEARCH/Publications/Publications-submitted-in-review/2020_Christiasen_Causality-conflict-forest-lss/Maps/"
roadSHP = drvMemV.CopyDataSource(ogr.Open(workFolder + "SHP-Files/COL_roads.shp"),'')
FL = gdal.Open(workFolder + "COL_Cummulated_ForestLoss_2001-2018.tif")
# ####################################### PROCESSING ########################################################## #
# (0) Reproject road layer into common projection
road_lyr = roadSHP.GetLayer()
# Build a coordinate transformation
from_SR = road_lyr.GetSpatialRef()
to_SR = osr.SpatialReference()
to_SR.ImportFromWkt(FL.GetProjection())
tr = osr.CoordinateTransformation(from_SR, to_SR)
# Create new layer --> add attribute "Value" and set to 1
road_proj = drvMemV.CreateDataSource('')
road_proj_lyr = road_proj.CreateLayer('', geom_type=ogr.wkbLineString)
road_proj_lyr.CreateField(ogr.FieldDefn("Value", ogr.OFTInteger))
road_proj_lyrDefn = road_proj_lyr.GetLayerDefn()
# Transform and copy each feature into the new layer
inFeat = road_lyr.GetNextFeature()
while inFeat:
    geom = inFeat.GetGeometryRef()
    geom.Transform(tr)
    # create new feature
    outFeat = ogr.Feature(road_proj_lyrDefn)
    # set the geometry and attribute
    outFeat.SetGeometry(geom)
    outFeat.SetField("Value", 1)
    # add the feature to the shapefile
    road_proj_lyr.CreateFeature(outFeat)
    # Get the next feature
    inFeat = road_lyr.GetNextFeature()
# adjust the geotranform to achieve a finer resolution at first
in_gt = FL.GetGeoTransform()
out_gt = [in_gt[0], in_gt[1]/10, 0, in_gt[3], 0, in_gt[5]/10]

# (1) Rasterize the projected road shapefile
road_ras = drvMemR.Create('', FL.RasterXSize*10, FL.RasterYSize*10, 1, gdal.GDT_Byte)
road_ras.SetProjection(FL.GetProjection())
road_ras.SetGeoTransform(out_gt)
gdal.RasterizeLayer(road_ras, [1], road_proj_lyr, options=["ATTRIBUTE=Value"])
# (2) Calculate the distance raster
road_dist = drvMemR.Create('', FL.RasterXSize*10, FL.RasterYSize*10, 1, gdal.GDT_Float32)
road_dist.SetProjection(FL.GetProjection())
road_dist.SetGeoTransform(out_gt)
gdal.ComputeProximity(road_ras.GetRasterBand(1), road_dist.GetRasterBand(1), options=['VALUES=1', 'DISTUNITS=GEO'])
# (4) Copy new raster to disc
drvR.CreateCopy(workFolder + "DistanceToRoad_1km.tif", road_dist)
# ####################################### END TIME-COUNT AND PRINT TIME STATS################################## #
print("")
endtime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print("--------------------------------------------------------")
print("--------------------------------------------------------")
print("start: " + starttime)
print("end: " + endtime)
print("")