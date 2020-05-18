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
    rootFolder = "D:/baumamat/Warfare/"
    pol_shp = rootFolder + "_SHPs/BIOMES_TropicsSavannas_10kmGrid_polygons.shp"
    point_shp = rootFolder + "_Variables/CONFLICT_DATA/PRIO/Disaggregated_Data/UCDP Georeferenced Event Dataset (GED) Global version 18.1 (2018)/ged181.shp"
    out_csv = rootFolder + "_DataSummaries/PRIO-summaries_highEstimate_20191217.csv"
    yrs = range(2000, 2019)
    nPackages = 1000
    nr_cores = 50
# ####################################### PROCESSING ########################################################## #
# (1) Build job list
    jobList = []
    # Get the number of total features in the shapefile
    eco = ogr.Open(pol_shp)
    ecoLYR = eco.GetLayer()
    nFeat = ecoLYR.GetFeatureCount()
    # Create a list of UIDs and subdivide the into smaller chunks
    featIDs = list(range(1, nFeat+1, 1))
    packageSize = int(nFeat / nPackages)
    #
    IDlist = [featIDs[i * packageSize:(i + 1) * packageSize] for i in range((len(featIDs) + packageSize - 1) // packageSize )]
    #
    # Now build the jobs and append to job list
    for chunk in IDlist:
        job = {'ids': chunk,
               'pol_path': pol_shp,
               'point_path': point_shp,
               'years': yrs}
        jobList.append(job)
# (2) Build Worker_Function
    def SumFunc(job):
    # Prepare the stuff we need for the processing of the data
        # Define the drivers that we need for creating the summaries
        drvMemV = ogr.GetDriverByName('Memory')
        # Load the polygon layer into mem, get the layer and subset by the IDs that are in the chunk
        grids = bt.baumiVT.CopyToMem(job['pol_path'])
        pol_lyr = grids.GetLayer()
        grids_copy = bt.baumiVT.CopyToMem(job['pol_path'])
        pol_lyr_copy = grids_copy.GetLayer()
        idSubs = job['ids']
        pol_lyr.SetAttributeFilter("UniqueID IN {}".format(tuple(idSubs)))
        # Load the conflict layer into memory, get the layer
        conflicts = bt.baumiVT.CopyToMem(job['point_path'])
        conf_lyr = conflicts.GetLayer()
        # Create coordinate transformation rule
        grid_SR = pol_lyr.GetSpatialRef()
        conflict_SR = conf_lyr.GetSpatialRef()
        trans = osr.CoordinateTransformation(grid_SR, conflict_SR)
        # Define the output-list that we want to return
        outList = []
    # Now loop through the selected features in our lyr
        feat = pol_lyr.GetNextFeature()
        while feat:
        # Get the UID from the field
            UID = feat.GetField("UniqueID")
    # Get the geometry, and start running the analyses
            geom = feat.GetGeometryRef()
            # Create a clone of the geometry, transform it to the point-CS, set spatial filter on conflict.lyr
            geom_cl = geom.Clone()
            geom_cl.Transform(trans)
            conf_lyr.SetSpatialFilter(geom_cl)
            # Loop through the years; make in addition thematic selection on the precision (where_prec <= 2)
            for yr in job['years']:
                select_statement = '(where_prec < 3) AND (year = ' + str(yr) + ')'
                conf_lyr.SetAttributeFilter(select_statement)
                # Calculate the statistics we want to gather
        # (1) FOR THIS GRID-CELL ONLY
                # (a) Nr. events
                nr_events = conf_lyr.GetFeatureCount()
                # Calculate the rest only if nr_events > 0, otherwise set to zero
                if nr_events > 0:
                # (b) Nr. fatalities
                    nr_fatalities = 0
                    conf_feat = conf_lyr.GetNextFeature()
                    while conf_feat:
                        deaths = conf_feat.GetField('best')
                        nr_fatalities = nr_fatalities + deaths
                        conf_feat = conf_lyr.GetNextFeature()
                # (c) Is it "conflict" by PRIO-definition --> more than 25 deaths per year
                    if nr_fatalities >= 25:
                        conflict_yn = 1
                    else:
                        conflict_yn = 0
                # (d) Average # fatalities per events
                    if nr_events > 0:
                        fat_per_event = nr_fatalities / nr_events
                    else:
                        fat_per_event = 0
                else:
                    nr_fatalities = 0
                    conflict_yn = 0
                    fat_per_event = 0
                # Reset the reading and remove attribute filter
                conf_lyr.ResetReading()
                conf_lyr.SetAttributeFilter(None)
                 # Add to the list we want to store
                vals = [UID, yr, conflict_yn, nr_events, nr_fatalities, format(fat_per_event, '.3f')]
				# Append the values to the output-DS
                outList.append(vals)
            # after all years are processed, take the next feature
            feat = pol_lyr.GetNextFeature()
    # return the outList as output from the function
        return outList
# (3) Execute the Worker_Funtion parallel
    job_results = Parallel(n_jobs=nr_cores)(delayed(SumFunc)(i) for i in tqdm(jobList))
# (4) Merge the different packages back together into one dataset, instantiate colnames first
    print("Merge Outputs")
    outDS = [["PolygonID", "Year", "Conflict_YN", "nr_events", "nr_fatalities", "nr_fatalities_per_event"]]
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