# ####################################### LOAD REQUIRED LIBRARIES ############################################# #
import os, csv
import time
import gdal, osr, ogr
import numpy as np
import baumiTools as bt
from joblib import Parallel, delayed
from tqdm import tqdm
# ####################################### SET TIME-COUNT ###################################################### #
if __name__ == '__main__':
    starttime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("--------------------------------------------------------")
    print("Starting process, time:" +  starttime)
    print("")
# ####################################### FOLDER PATHS AND BASIC VARIABLES FOR PROCESSING ##################### #
# Tiles from the global forest loss dataset are being merged into a virtual raster dataset (.vrt) using the gdal_buildvrt function
    rootFolder = "Z:/Warfare/"
    pol_shp = rootFolder + "_SHPs/BIOMES_TropicsSavannas_10kmGrid_polygons_COL.shp"
    out_csv = rootFolder + "_DataSummaries/RoadDist_Summaries.csv"
    roadDist = rootFolder + "_Variables/DistanceToRoad_1km.tif"
    nPackages = 20
    nr_cores = 10
# ####################################### PROCESSING ########################################################## #
# (1) Build job list
    jobList = []
    # Get the number of total features in the shapefile
    eco = ogr.Open(pol_shp)
    ecoLYR = eco.GetLayer()
    nFeat = ecoLYR.GetFeatureCount()
    # Create a list of UIDs and subdivide the into smaller chunks
    #featIDs = list(range(1,nFeat+1, 1))
    featIDs = []
    for feat in ecoLYR:
        featIDs.append(feat.GetField("UniqueID"))
    packageSize = int(nFeat / nPackages)
    #
    IDlist = [featIDs[i * packageSize:(i + 1) * packageSize] for i in range((len(featIDs) + packageSize - 1) // packageSize )]
    # Now build the jobs and append to job list
    for chunk in IDlist:
        job = {'ids': chunk,
               'shp_path': pol_shp,
               'rodDist_raster': roadDist}
        jobList.append(job)
# (2) Build Worker_Function
    def SumFunc(job):
    # Define the drivers that we need for creating the summaries
        drvMemV = ogr.GetDriverByName('Memory')
        drvMemR = gdal.GetDriverByName('MEM')
        # Load the shapefile into mem, get the layer and subset by the IDs that are in the chunk
        shpMem = bt.baumiVT.CopyToMem(job['shp_path'])
        lyr = shpMem.GetLayer()
        idSubs = job['ids']
        lyr.SetAttributeFilter("UniqueID IN {}".format(tuple(idSubs)))
        # Create coordinate transformation rule
        pol_SR = lyr.GetSpatialRef()
        # Define the output-list that we want to return
        outList = []
        # Now loop through the selected features in our lyr
        feat = lyr.GetNextFeature()
        while feat:
    # Get needed properties from the SHP-File, the take geometry, and transform to Target-EPSG
            # UID Info
            UID = feat.GetField("UniqueID")
    # Instantiate output and take the geometry of the feature, transform it to our epsg
            vals = [UID]
            geom = feat.GetGeometryRef()
    # Rasterize the geometry, pixelSize is 30m
        # Create new SHP-file in memory to which we copy the geometry
            geom_shp = drvMemV.CreateDataSource('')
            geom_lyr = geom_shp.CreateLayer('geom', pol_SR, geom_type=ogr.wkbMultiPolygon)
            geom_lyrDefn = geom_lyr.GetLayerDefn()
            geom_feat = ogr.Feature(geom_lyrDefn)
            geom_feat.SetGeometry(geom)
            geom_lyr.CreateFeature(geom_feat)
        # Check if the geometry we are processing is larger than 1x1 pixel
            x_min, x_max, y_min, y_max = geom_lyr.GetExtent()
            x_res = int((x_max - x_min) / 30)
            y_res = int((y_max - y_min) / 30)
            geom_ras = drvMemR.Create('', x_res, y_res, gdal.GDT_Byte)
            geom_ras.SetProjection(pol_SR.ExportToWkt())
            geom_ras.SetGeoTransform((x_min, 30, 0, y_max, 0, -30))
            gdal.RasterizeLayer(geom_ras, [1], geom_lyr, burn_values=[1])
        # Reproject the Hansen-Rasters "into" the geometry-raster
            def ReprojectRaster(valRaster, GEOMraster):
                vasRaster_sub = drvMemR.Create('', GEOMraster.RasterXSize, GEOMraster.RasterYSize, 1, gdal.GDT_Byte)
                vasRaster_sub.SetGeoTransform(GEOMraster.GetGeoTransform())
                vasRaster_sub.SetProjection(GEOMraster.GetProjection())
                gdal.ReprojectImage(valRaster, vasRaster_sub, valRaster.GetProjection(), GEOMraster.GetProjection(), gdal.GRA_NearestNeighbour)
                return vasRaster_sub
            rd = ReprojectRaster(gdal.Open(job['rodDist_raster']), geom_ras)
        # Open all rasters into np-Arrays
            geom_np = geom_ras.GetRasterBand(1).ReadAsArray(0, 0, x_res, y_res)
            roadDist_np = rd.GetRasterBand(1).ReadAsArray(0, 0, x_res, y_res)
        # Now extract the summaries
            rd_m = np.where((geom_np == 1), roadDist_np, 0)
            mean = np.mean(rd_m)
            vals.append(format(mean, '.5f'))
            
        # Append the values to the output-DS, then take the next feature
            outList.append(vals)
            feat = lyr.GetNextFeature()
    # return the outList as output from the function
        return outList
# (3) Execute the Worker_Funtion parallel
    job_results = Parallel(n_jobs=nr_cores)(delayed(SumFunc)(i) for i in tqdm(jobList))

# (4) Merge the different packages back together into one dataset, instantiate colnames first
    print("Merge Outputs")
    outDS = [["UniqueID",
              "RoadDist"]]
    # Now extract the information from all the evaluations
    # 1st loop --> the different chunks
    for result in job_results:
        # 2nd loop --> all outputs in each chunk
        for out in result:
            outDS.append(out)
# (5) Write all outputs to disc
    print("Write output")
    with open(out_csv, "w") as theFile:
        csv.register_dialect("custom", delimiter = ",", skipinitialspace = True, lineterminator = '\n')
        writer = csv.writer(theFile, dialect = "custom")
        for element in outDS:
            writer.writerow(element)
# ####################################### END TIME-COUNT AND PRINT TIME STATS################################## #
    print("")
    endtime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("--------------------------------------------------------")
    print("--------------------------------------------------------")
    print("start: " + starttime)
    print("end: " + endtime)
    print("")