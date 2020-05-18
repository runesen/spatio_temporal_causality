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
    rootFolder = "Z:/Warfare/"
    pol_shp = rootFolder + "_SHPs/BIOMES_TropicsSavannas_10kmGrid_polygons.shp"
    out_csv = rootFolder + "_DataSummaries/GPWv4_summary_20190129.csv"
    pop00 = rootFolder + "_Variables/GPWv4/gpw_v4_population_density_rev10_2000_30_sec.tif"
    pop05 = rootFolder + "_Variables/GPWv4/gpw_v4_population_density_rev10_2005_30_sec.tif"
    pop10 = rootFolder + "_Variables/GPWv4/gpw_v4_population_density_rev10_2010_30_sec.tif"
    pop15 = rootFolder + "_Variables/GPWv4/gpw_v4_population_density_rev10_2015_30_sec.tif"
    nPackages = 300
    nr_cores = 60
# ####################################### PROCESSING ########################################################## #
# (1) Build job list
    jobList = []
    # Get the number of total features in the shapefile
    pol = ogr.Open(pol_shp)
    polLYR = pol.GetLayer()
    nFeat = polLYR.GetFeatureCount()
    # Create a list of UIDs and subdivide the into smaller chunks
    featIDs = list(range(1,nFeat+1, 1))
    packageSize = int(nFeat / nPackages)
    #
    IDlist = [featIDs[i * packageSize:(i + 1) * packageSize] for i in range((len(featIDs) + packageSize - 1) // packageSize )]
    # Now build the jobs and append to job list
    for chunk in IDlist:
        job = {'ids': chunk,
               'shp_path': pol_shp,
               'pop00': pop00,
               'pop05': pop05,
               'pop10': pop10,
               'pop15': pop15}
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
    # Rasterize the geometry, pixelSize is 1000m
        # Create new SHP-file in memory to which we copy the geometry
            geom_shp = drvMemV.CreateDataSource('')
            geom_lyr = geom_shp.CreateLayer('geom', pol_SR, geom_type=ogr.wkbMultiPolygon)
            geom_lyrDefn = geom_lyr.GetLayerDefn()
            geom_feat = ogr.Feature(geom_lyrDefn)
            geom_feat.SetGeometry(geom)
            geom_lyr.CreateFeature(geom_feat)
        # Check if the geometry we are processing is larger than 1x1 pixel
            x_min, x_max, y_min, y_max = geom_lyr.GetExtent()
            x_res = int((x_max - x_min) / 1000)
            y_res = int((y_max - y_min) / 1000)
        # Do the rest of the operation for this polygon only if x_res and y_res are >= 1
            #if x_res > 0 and y_res > 0:
            geom_ras = drvMemR.Create('', x_res, y_res, gdal.GDT_Byte)
            geom_ras.SetProjection(pol_SR.ExportToWkt())
            geom_ras.SetGeoTransform((x_min, 1000, 0, y_max, 0, -1000))
            gdal.RasterizeLayer(geom_ras, [1], geom_lyr, burn_values=[1])
        # Reproject the Hansen-Rasters "into" the geometry-raster
            p00 = gdal.Open(job['pop00'])
            p00rb = p00.GetRasterBand(1)
            noData = p00rb.GetNoDataValue()
            p00 = None
            p00rb = None
            def ReprojectRaster(valRaster, GEOMraster):
                vasRaster_sub = drvMemR.Create('', GEOMraster.RasterXSize, GEOMraster.RasterYSize, 1, gdal.GDT_Float32)
                vasRaster_sub.SetGeoTransform(GEOMraster.GetGeoTransform())
                vasRaster_sub.SetProjection(GEOMraster.GetProjection())
                gdal.ReprojectImage(valRaster, vasRaster_sub, valRaster.GetProjection(), GEOMraster.GetProjection(), gdal.GRA_NearestNeighbour)
                return vasRaster_sub
            p00 = ReprojectRaster(gdal.Open(job['pop00']), geom_ras)
            p05 = ReprojectRaster(gdal.Open(job['pop05']), geom_ras)
            p10 = ReprojectRaster(gdal.Open(job['pop10']), geom_ras)
            p15 = ReprojectRaster(gdal.Open(job['pop15']), geom_ras)
        # Open all rasters into np-Arrays
            p00_np = p00.GetRasterBand(1).ReadAsArray(0, 0, x_res, y_res)
            #print(p00_np.mean())
            p00_masked = np.ma.masked_equal(p00_np, noData)
            #print(np.nanmean(p00_masked))
            p05_np = p05.GetRasterBand(1).ReadAsArray(0, 0, x_res, y_res)
            p05_masked = np.ma.masked_equal(p05_np, noData)
            p10_np = p10.GetRasterBand(1).ReadAsArray(0, 0, x_res, y_res)
            p10_masked = np.ma.masked_equal(p10_np, noData)
            p15_np = p15.GetRasterBand(1).ReadAsArray(0, 0, x_res, y_res)
            p15_masked = np.ma.masked_equal(p15_np, noData)
        # Now extract the summaries
            vals.append(format(np.nanmean(p00_masked), '.5f'))
            vals.append(format(np.nanmean(p05_masked), '.5f'))
            vals.append(format(np.nanmean(p10_masked), '.5f'))
            vals.append(format(np.nanmean(p15_masked), '.5f'))
        # Append the values to the output-DS, then take the next feature
            outList.append(vals)
            feat = lyr.GetNextFeature()
    # return the outList as output from the function
        return outList
# (3) Execute the Worker_Funtion parallel
    job_results = Parallel(n_jobs=nr_cores)(delayed(SumFunc)(i) for i in tqdm(jobList))

# (4) Merge the different packages back together into one dataset, instantiate colnames first
    print("Merge Outputs")
    outDS = [["UniqueID", "PopDens_00", "PopDens_05", "PopDens_10", "PopDens_15"]]
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