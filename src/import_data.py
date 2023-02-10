
# TODO/QUESTION  - Should AK and HI be removed

from optparse import Values
import os, sys
from select import select
import pathlib
import time
import pandas as pd
import arcpy
import requests
from arcgis.features import GeoAccessor, GeoSeriesAccessor
from arcgis.geometry import Polygon
import math
import numpy as np
from pathlib import Path
# --------------------------------------------------------------------------------------------------

"""
This script loads data from different sources to find a point location of the parking lot and merge with visitation data. 

Create Date: 2022
Created By: The Volpe Center's Geospatial Team 
"""


BASE_DIR = Path(__file__).parents[1].resolve()

INPUT_DIR  = os.path.join(BASE_DIR, "data\\input")
OUTPUT_DIR = os.path.join(BASE_DIR, "data\\output")

#### Spatial Features
DOI_BLM_GDB = os.path.join(INPUT_DIR, r"DOI_BLM\BLM_National_Recreation_Site_pts\recs.gdb")
DOI_BR_GDB  = os.path.join(INPUT_DIR, r"DOI_BR\Jan_2021_FLMA_Submittal_UPDATE2021June22.gdb")
DOI_FWS_SHP = os.path.join(INPUT_DIR, r"DOI_FWS\FWS_Public_Parking_Lot_All\FWS_Public_Parking_Lot_All.shp")
DOI_FWS_SHP_AMEN = os.path.join(INPUT_DIR, r"DOI_FWS\FWS_NWRS_HQ_VisitorServiceAmenities_V2\FWS_NWRS_HQ_VisitorServiceAmenities_V2.shp") #This is used to get the CCCODE
DOI_NPS_SHP = os.path.join(INPUT_DIR, r"DOI_NPS\NPS_-_Parking_Areas_-_Web_Mercator\NPS_-_Parking_Areas_-_Web_Mercator.shp")
DOI_NPS_GDB = os.path.join(INPUT_DIR, r"DOI_BR\Jan_2021_FLMA_Submittal_UPDATE2021June22.gdb")
USDA_FS_XLS = os.path.join(INPUT_DIR, r"USDA_FS\FS_Points_of_Interest_mod.xlsx") #Modified so that it is easier to load e.g. header row
USACE_SHP = os.path.join(INPUT_DIR, r"USACE\USACE_Recreation_Project_Site_Areas_(Points).shp")

PAD_BOUNDARY = os.path.join(INPUT_DIR,r"PAD\PAD_US2_1.gdb\PADUS 2.1 Combined  Fee, Designation, Easement")


STATE_BOUNDARIES = os.path.join(INPUT_DIR, r"Mapping_Info.gdb\all_states_natural_earth")
SDF_STATES = pd.DataFrame.spatial.from_featureclass(STATE_BOUNDARIES)
#### Visitation files

DOI_BLM_BND_GDB = os.path.join(INPUT_DIR, r"VISITATION\DOI_BLM\admu.gdb")
DOI_BR_BND = os.path.join(INPUT_DIR, r"VISITATION\DOI_BOR\Area_Office_Boundaries.shp")
DOI_FWS_BND = os.path.join(INPUT_DIR, r"VISITATION\DOI_FWS\FWSBoundaries.shp")
DOI_NPS_BND = os.path.join(INPUT_DIR, r"VISITATION\DOI_NPS\NPS_-_Land_Resources_Division_Boundary_and_Tract_Data_Service.shp")
USDA_FS_BND = os.path.join(INPUT_DIR, r"VISITATION\USDA_FS\AdministrativeForest.shp") 

DOI_BLM_CENT = "doi_blm_centroid"
DOI_BR_CENT = "doi_bor_centroid"
DOI_FWS_CENT = "doi_fws_centroid"
DOI_NPS_CENT= "doi_nps_centroid"
USDA_FS_CENT = "usda_fs_centroid"

DOI_BLM_VISIT = os.path.join(INPUT_DIR, r"VISITATION\DOI_BLM\Visitation_Cleaned_blm.xlsx")
DOI_BR_VISIT = os.path.join(INPUT_DIR, r"VISITATION\DOI_BOR\Visitation_Cleaned_bor.xlsx")
DOI_FWS_VISIT = os.path.join(INPUT_DIR, r"VISITATION\DOI_FWS\Visitation_Cleaned_fws.xlsx")
DOI_NPS_VISIT = os.path.join(INPUT_DIR, r"VISITATION\DOI_NPS\Visitation_Cleaned_nps.xlsx")
USDA_FS_VISIT = os.path.join(INPUT_DIR, r"VISITATION\USDA_FS\Visitation_Cleaned_usfs.xlsx") 
VISITATION_FILES = [DOI_BLM_VISIT,DOI_BR_VISIT,DOI_FWS_VISIT,DOI_NPS_VISIT,USDA_FS_VISIT]
#Or use an API Key
EV_STATIONS_CSV = os.path.join(INPUT_DIR, r"EV_stations\alt_fuel_stations (Feb 18 2022).csv")

#Units in Meters
MAX_SYMBOL_SIZE = 90000
MIN_SYMBOL_SIZE = 10000 

# --------------------------------------------------------------------------------------------------
def import_doi_blm(parts:list,centroids:list):
    """load data from blm source file
    Parameters:
        parts (list): data will be appended to this list object, Order [shape,orgcode,name,agency,matchid]
        centroids (list): data will be appended to this list object, Order [shape,matchid,name,visitation]
    returns:
        parts (list): with appended data"""

    print('\nImporting DOI BLM ...')

    # 1 information center
    # 8 parking area
    # TODO - seems to me there would be more to include (e..g campgrounds)


    lyr = arcpy.SelectLayerByAttribute_management(
        os.path.join(DOI_BLM_GDB, r"RECS\recs_pt"),
        'NEW_SELECTION',
        "ADMIN_ST NOT IN ('AK') and  FET_TYPE IN (8, 1)"
        )

    arcpy.Project_management(lyr, "doi_blm", arcpy.SpatialReference("WGS 1984"))

    #arcpy.AddField_management("doi_blm", "agency", "TEXT", "", "", 10)
    #arcpy.AddField_management("doi_blm", "visitfld", "TEXT", "", "", 10)
    #arcpy.CalculateField_management("doi_blm", "agency", "'doi_blm'", "PYTHON3", '')
    #arcpy.CalculateField_management("doi_blm", "visitfld", "!ADM_UNIT_CD!", "PYTHON3", '')
    with arcpy.da.SearchCursor("doi_blm",["SHAPE@","ADM_UNIT_CD","FET_NAME"]) as sc:
        [parts.append([row[0],row[1],row[2],"doi_blm",row[1]]) for row in sc]

    cent_fc = project(DOI_BLM_BND_GDB+"\\blm_natl_admu_labels_webpub",DOI_BLM_CENT)

    with arcpy.da.SearchCursor(cent_fc, ["SHAPE@","ADM_UNIT_CD","Label_Full_Name"]) as sc:
        [centroids.append([row[0],row[1],row[2],"doi_blm",-1,row[1]]) for row in sc]

    print('done')
    return parts,centroids

# --------------------------------------------------------------------------------------------------
def import_doi_br(parts:list,centroids:list):
    """load data from blm source file
    Parameters:
        parts (list): data will be appended to this list object, Order [shape,orgcode,name,agency,matchid]
    returns:
        parts (list): with appended data"""
    print('\nImporting DOI BR ...')

    lyr = arcpy.SelectLayerByAttribute_management(
        os.path.join(DOI_BR_GDB, "BOR_2021"),
        'NEW_SELECTION',
        "FLTP IN ('YES FLTP NETWORK')"
        )
        
    domain_codes = {}
    domains = arcpy.da.ListDomains(DOI_BR_GDB)
    for domain in domains:
        if domain.name == "domAreaOffice":
            coded_values = domain.codedValues
            for val, desc in coded_values.items():
                domain_codes[val]=desc

    # TODO figure out how to trim down these vertices
    arcpy.FeatureVerticesToPoints_management(
        lyr,
        "doi_bor_road_vertices",
        point_location="BOTH_ENDS"
        )

    with arcpy.da.SearchCursor("doi_bor_road_vertices",["SHAPE@",'FLMA_UNIT_ID', 'RTE_NAME']) as sc:
        [parts.append([row[0],row[1],row[2],"doi_bor",domain_codes[row[1]]]) for row in sc]

    cent_fc = project_and_create_centroid_fc(DOI_BR_BND,DOI_BR_CENT)
    with arcpy.da.SearchCursor(cent_fc, ["SHAPE@","AreaOffice"]) as sc:
        [centroids.append([row[0],row[1],row[1],"doi_bor",-1,row[1]]) for row in sc]
    #arcpy.AddField_management("doi_bor_road_vertices", "agency", "TEXT", "", "", 10)

    #arcpy.CalculateField_management("doi_bor_road_vertices", "agency", "'doi_bor'", "PYTHON3", '')

    print('done')
    return parts,centroids

# --------------------------------------------------------------------------------------------------
def import_doi_fws(parts:list,centroids:list):
    """load data from fws source file
    Parameters:
        parts (list): data will be appended to this list object, Order [shape,orgcode,name,agency,matchid]
    returns:
        parts (list): with appended data"""

    print('\nImporting DOI FWS ...')

    arcpy.management.FeatureToPoint(
        DOI_FWS_SHP,
        "doi_fws",
        "INSIDE")

    #arcpy.AddField_management("doi_fws", "agency", "TEXT", "", "", 10)

    #arcpy.CalculateField_management("doi_fws", "agency", "'doi_fws'", "PYTHON3", '')
    with arcpy.da.SearchCursor(DOI_FWS_SHP_AMEN,["OrgCode","OrgName","CCCODE"]) as sc:
        codes = {str(row[0]):(row[1],row[2]) for row in sc}

    with arcpy.da.SearchCursor("doi_fws",["SHAPE@",'Org_Code',"Route_Name"]) as sc:
        for row in sc:
            try:
                cc = codes[row[1]][1]
                lbl = codes[row[1]][0]
            except:
                cc = None
                lbl = row[2]
            parts.append([row[0],row[1],lbl,"doi_fws",cc])

    cent_fc = project_and_create_centroid_fc(DOI_FWS_BND,DOI_FWS_CENT)

    with arcpy.da.SearchCursor(cent_fc, ["SHAPE@","CostCenter","ORGNAME","ORGCODE"]) as sc:
        [centroids.append([row[0],row[1],row[2],"doi_fws",-1,row[3]]) for row in sc]

    print('done')
    return parts,centroids

# --------------------------------------------------------------------------------------------------
def import_doi_nps(parts:list,centroids:list):
    """load data from nps source file
    Parameters:
        parts (list): data will be appended to this list object, Order [shape,orgcode,name,agency,matchid]
    returns:
        parts (list): with appended data"""

    print('\nImporting DOI NPS ...')

    lyr = arcpy.SelectLayerByAttribute_management(
        os.path.join(DOI_NPS_GDB, "NPS_PAVED_2021"),
        'NEW_SELECTION',
        "FLTP IN ('YES FLTP NETWORK')"
        )
    arcpy.FeatureVerticesToPoints_management(
        lyr,
        "doi_nps_road_vertices",
        point_location="BOTH_ENDS"
        )

    arcpy.management.FeatureToPoint(
        DOI_NPS_SHP,
        "doi_nps_tmp",
        "INSIDE")

    DOI_NPS_GDB

    arcpy.Project_management("doi_nps_tmp", "doi_nps", arcpy.SpatialReference("WGS 1984"))

    arcpy.Delete_management("doi_nps_tmp")


    #arcpy.AddField_management("doi_nps", "agency", "TEXT", "", "", 10)
    #arcpy.CalculateField_management("doi_nps", "agency", "'doi_nps'", "PYTHON3", '')
    with arcpy.da.SearchCursor("doi_nps",["SHAPE@",'UNITCODE', 'LOTNAME']) as sc:
        [parts.append([row[0],row[1],row[2],"doi_nps",row[1]]) for row in sc]

    with arcpy.da.SearchCursor("doi_nps_road_vertices",["SHAPE@",'FLMA_UNIT_ID', 'RTE_NAME']) as sc:
        [parts.append([row[0],row[1],row[2],"doi_nps",row[1]]) for row in sc]

    cent_fc = project_and_create_centroid_fc(DOI_NPS_BND,DOI_NPS_CENT)
    with arcpy.da.SearchCursor(cent_fc, ["SHAPE@","UNIT_CODE","UNIT_NAME"]) as sc:
        [centroids.append([row[0],row[1],row[2],"doi_nps",-1,row[1]]) for row in sc]

    print('done')
    return parts,centroids
# --------------------------------------------------------------------------------------------------
def import_usda_fs(parts:list,centroids:list):
    """load data from usda fs source file
    Parameters:
        parts (list): data will be appended to this list object, Order [shape,orgcode,name,agency,matchid]
    returns:
        parts (list): with appended data"""

    print('\nImporting USDA FS ...')

    df = pd.read_excel(USDA_FS_XLS, sheet_name='Details',converters={"SECURITY_ID":str,"MANAGING_ORG":str})
    df = df[df["EV_POTENTIAL"] == "Y"]

    cpfc = create_point_feature_class("usda_fs")

    if cpfc == False:
        raise Exception("Error creating the point feature class for the usda_fs")

    icursor = arcpy.da.InsertCursor("usda_fs", ["SHAPE@XY", "org", "name","agency","matchid"])

    for index, row in df[~(pd.isnull(df["LATITUDE"]))].iterrows():
        icursor.insertRow([(row['LONGITUDE'], row['LATITUDE']), row["SECURITY_ID"], row['REC_SITE_or_OFFICE_NAME'], "usda_fs", row["SECURITY_ID"]])
        #icursor.insertRow([(row['LONGITUDE'], row['LATITUDE']), row['MANAGING_ORG'], row['REC_SITE_or_OFFICE_NAME'], "usda_fs", row["SECURITY_ID"]])

    del icursor

    with arcpy.da.SearchCursor("usda_fs", ["SHAPE@", "org", "name","agency","matchid"]) as sc:
        [parts.append([row[0],row[1],row[2],row[3],row[4]]) for row in sc]

    cent_fc = project_and_create_centroid_fc(USDA_FS_BND,USDA_FS_CENT)
    with arcpy.da.SearchCursor(cent_fc, ["SHAPE@","FORESTORGC","FORESTNAME"]) as sc:
        [centroids.append([row[0],row[1],row[2],"usda_fs",-1,row[1]]) for row in sc]


    print('done')
    return parts,centroids


def import_usace(parts:list):
    """load data from usace source file
    Parameters:
        parts (list): data will be appended to this list object, Order [shape,orgcode,name,agency,matchid]
    returns:
        parts (list): with appended data"""

    print('\nImporting USACE ...')

    arcpy.Project_management(USACE_SHP, "usace", arcpy.SpatialReference("WGS 1984"))

    
    with arcpy.da.SearchCursor("usace", ["SHAPE@", "recproje_3", "featurenam"]) as sc:
        [parts.append([row[0],row[1],row[2],"usace",""]) for row in sc ]




    print('done')
    return parts
# --------------------------------------------------------------------------------------------------




def merge(parts:list,collapse_feet:int):
    """
        combines all point information within the parts list of lists.
        Parameters:
            parts (list): [shape, org, name, agency, matchid]
            collapse_feet (int): used for minimum distance to consider points as unique.
        Returns:
            None
    """
    print('\nMerging flma point layers ...')

    # USDA_FS has the desired schema, start with that
    cpfc = create_point_feature_class("flma_points")
    if cpfc == False:
        raise Exception("Error creating the point feature class for the usda_fs")
    
    with arcpy.da.InsertCursor("flma_points",["SHAPE@", "org", "name", "agency","matchid"]) as ic:
        for p in parts:
            if p[0]:
                #print(p[2])
                try:
                    ic.insertRow(p)
                except:
                    pass
        
    print("\nCollapsing Points using {0} feet".format(collapse_feet))
    arcpy.management.DeleteIdentical("flma_points",["SHAPE","matchid","agency"],"{0} Feet".format(collapse_feet))#


    
    


    print('done')
# --------------------------------------------------------------------------------------------------

def import_ev_stations():

    print('\ncreating ev station layer')

    use_cols = ['Latitude', 'Longitude', 'EV DC Fast Count', 'Access Code']

    df = pd.read_csv(EV_STATIONS_CSV, usecols=use_cols, low_memory = False)

    arcpy.CreateFeatureclass_management(arcpy.env.workspace, "charging_stations", "POINT",
            "", "DISABLED", "DISABLED", arcpy.SpatialReference("WGS 1984"))

    arcpy.AddField_management("charging_stations", "ev_dc_fast_count", 'long')
    arcpy.AddField_management("charging_stations", "access_code", "TEXT", "", "", 10)


    icursor = arcpy.da.InsertCursor("charging_stations", ["SHAPE@XY", "ev_dc_fast_count", "access_code"])

    for index, row in df.iterrows():
        icursor.insertRow([(row['Longitude'], row['Latitude']), row['EV DC Fast Count'], row['Access Code']])

    del icursor
    print('done')


def import_ev_stations_from_api(key_str:str):
    """Load Charging stations using the API. Filters for connection type and creates Level 2 and DCFC classes. Writes to geodatabase.
    Parameters:
        key_str (str): personal API key requested from https://developer.nrel.gov/signup/
    Returns:
        Nothing"""
    print('\nRequesting EV station data from API')

    connector_types = ["J1772","J1772COMBO","CHADEMO"]
    api_url = "https://developer.nrel.gov/api/alt-fuel-stations/v1.json?fuel_type=ELEC&access=public&status=all&ev_connector_type={1}&api_key={0}".format(key_str,",".join(connector_types))

    data = requests.get(api_url).json()["fuel_stations"]
    df = pd.DataFrame.from_dict(data)
    df["ev_connector_types_level"] = df["ev_connector_types"].apply(lambda x: check_ct(x))
    df["ev_connector_types_splt"] = df["ev_connector_types"].apply(lambda x: "|".join(x))
    use_cols = ['latitude', 'longitude', 'ev_dc_fast_num', 'access_code','ev_connector_types_level','ev_connector_types_splt']
    print("\nwriting station layer")
    

    arcpy.CreateFeatureclass_management(arcpy.env.workspace, "charging_stations", "POINT",
            "", "DISABLED", "DISABLED", arcpy.SpatialReference("WGS 1984"))

    arcpy.AddField_management("charging_stations", "ev_dc_fast_count", 'long')
    arcpy.AddField_management("charging_stations", "access_code", "TEXT", "", "", 10)
    arcpy.AddField_management("charging_stations", "ev_connector_types_lev", "TEXT", "", "", 50)
    arcpy.AddField_management("charging_stations", "ev_connector_types", "TEXT", "", "", 100)
    with arcpy.da.InsertCursor("charging_stations", ["SHAPE@XY", "ev_dc_fast_count", "access_code","ev_connector_types","ev_connector_types_lev"]) as ic:
        for index, row in df.iterrows():
            ic.insertRow([(row['longitude'], row['latitude']), row['ev_dc_fast_num'],  row['access_code'], row['ev_connector_types_splt'], row['ev_connector_types_level']])

def match_visitation_data_with_centroids(centroids):
    print("\nCombining visitation data...")
    sr=arcpy.SpatialReference(102003)
    arcpy.management.CreateTable(arcpy.env.workspace,"visitation_table")
    arcpy.management.AddField("visitation_table","matchid","TEXT",field_length=100)
    arcpy.management.AddField("visitation_table","matchsites","TEXT",field_length=100)
    arcpy.management.AddField("visitation_table","label","TEXT",field_length=150)
    arcpy.management.AddField("visitation_table","agency","TEXT",field_length=10)
    arcpy.management.AddField("visitation_table","visitation","LONG")

    with arcpy.da.InsertCursor("visitation_table",["matchid","label","agency","visitation","matchsites"]) as ic:
        for tbl in VISITATION_FILES:
            temp = pd.read_excel(tbl,sheet_name="Sheet1",
                    converters={'MATCHID':str,'VISITATION':int,'MATCHSITES':str})
            temp["VISITATION"] = temp["VISITATION"].fillna(0)
            for idx,row in temp.iterrows():
                ic.insertRow([row["MATCHID"],row["LABEL"],row["AGENCY"],row["VISITATION"],row['MATCHSITES']])
    cdf = pd.DataFrame(centroids,columns=["SHAPE","matchid","label","agency","visitation","matchsites"])
    #cdf.to_csv("c:/temp/centroids.csv")
    sdf = pd.DataFrame.spatial.from_df(cdf,geometry_column="SHAPE")
    sdf["sqrt"] = 0
    sdf["matched"] = 0
    #print(sdf.columns)
    visits = pd.DataFrame.spatial.from_table(arcpy.env.workspace + "\\visitation_table")
    for idx,row in visits.iterrows():
        sqroot = row['visitation']**(1/2)
        sdf.loc[(sdf["agency"]==row["agency"])&(sdf["matchid"]==row["matchid"]),"matched"] = 1
        sdf.loc[(sdf["agency"]==row["agency"])&(sdf["matchid"]==row["matchid"]),"visitation"] = row['visitation']
        sdf.loc[(sdf["agency"]==row["agency"])&(sdf["matchid"]==row["matchid"]),"sqrt"] = sqroot
    #print(sdf.columns)
    #print(sdf.head())
    #sdf.spatial.set_geometry("geometry")
    #sdf = sdf.reset_index()
    sdf.spatial.to_featureclass(arcpy.env.workspace + "\\flma_points_visits")

    try:
        project("flma_points_visits","flma_points_visits_prj")
    except:
        print("flma_points_visits_prj alreayd exists...")
    
    print("joining state information")
    fieldmappings = arcpy.FieldMappings()
    fieldmappings.addTable("flma_points_visits_prj")
    fieldmappings2 = arcpy.FieldMappings()
    fieldmappings2.addTable(STATE_BOUNDARIES)


    fieldmappings.addFieldMap(fieldmappings2.getFieldMap(fieldmappings2.findFieldMapIndex("name")))
    fieldmappings.addFieldMap(fieldmappings2.getFieldMap(fieldmappings2.findFieldMapIndex("postal")))
    fieldmappings.addFieldMap(fieldmappings2.getFieldMap(fieldmappings2.findFieldMapIndex("scale")))
    fieldmappings.addFieldMap(fieldmappings2.getFieldMap(fieldmappings2.findFieldMapIndex("maxsymb")))
    fieldmappings.addFieldMap(fieldmappings2.getFieldMap(fieldmappings2.findFieldMapIndex("minsymb")))
    flma_pnts = arcpy.SpatialJoin_analysis("flma_points_visits_prj",STATE_BOUNDARIES,"flma_points_visits_states","JOIN_ONE_TO_ONE","KEEP_ALL",fieldmappings,match_option="CLOSEST").getOutput(0)
    return flma_pnts

def buffer_by_agency_state(flma_pnts):
    print("Buffering for symbols...")
    sr=arcpy.SpatialReference(102003)
    sdf = pd.DataFrame.spatial.from_featureclass(flma_pnts)
    sdf['bufmeters'] = 0
    sdf.loc[sdf["sqrt"]>0,'bufmeters'] = MAX_SYMBOL_SIZE*sdf[sdf["sqrt"]>0]["sqrt"] / sdf[sdf["sqrt"]>0]["sqrt"].max()
    cutoff_value = MIN_SYMBOL_SIZE/MAX_SYMBOL_SIZE * sdf[sdf["sqrt"]>0]["sqrt"].max()
    sdf.loc[sdf["sqrt"]<=cutoff_value,'bufmeters'] = MIN_SYMBOL_SIZE

    sdf["bufbyagency"] = 0
    sdf["bufall"] = 0
    legend_data = {"SHAPE":[],"agency":[],"value":[],"state":[]}
    for state in sdf["name"].unique():
        temp = sdf[sdf["name"]==state]
        max_symb_size = temp['maxsymb'].max()
        min_symb_size = temp['minsymb'].min()
        
        for agency in temp["agency"].unique():

            sqrt_max = sdf[(sdf["agency"]==agency)&(sdf["sqrt"]>0)&(sdf['name']==state)]['sqrt'].max()
            values = sdf[(sdf["agency"]==agency)&(sdf["sqrt"]>0)&(sdf['name']==state)]['sqrt'].values
            sdf.loc[(sdf["agency"]==agency)&(sdf["sqrt"]>0)&(sdf['name']==state),'bufbyagency'] = max_symb_size * values  / sqrt_max
            cutoff_value = min_symb_size/max_symb_size* sqrt_max
            sdf.loc[(sdf["agency"]==agency)&(sdf["sqrt"]<=cutoff_value)&(sdf['name']==state),'bufbyagency'] = min_symb_size
            bvals = []
            value = sqrt_max * sqrt_max
            minvalue = cutoff_value*cutoff_value
            cnt = 0
            while cnt <2:
                for x in [10000000,1000000,100000,10000,1000]:
                    if value>=x:
                        rp = round_place(value,x)
                        bvals.append(rp)
                        break

                value = value *.25
                if value<minvalue:
                    break
                cnt +=1
            for x in [1000000,100000,10000,1000]:
                if cutoff_value*cutoff_value>=x:
                    bvals.append(round_place(cutoff_value*cutoff_value,x))
                    break

            bds = [max_symb_size * math.sqrt(x) / sqrt_max for x in bvals]

            for bd,x in zip(bds,bvals):
                legend_data["agency"].append(agency)
                legend_data["value"].append("{:,}".format(x))
                geom = arcpy.PointGeometry(arcpy.Point(0,bd+max(bds)), sr).buffer(bd)
                legend_data["SHAPE"].append(geom)
                legend_data["state"].append(state)
            if "<=" not in legend_data["value"][-1]:
                legend_data["value"][-1] = "<=" + legend_data["value"][-1] 

        sqrt_max = sdf[(sdf["sqrt"]>0)&(sdf['name']==state)]['sqrt'].max()
        values = sdf[(sdf["sqrt"]>0)&(sdf['name']==state)]['sqrt'].values
        sdf.loc[(sdf["sqrt"]>0)&(sdf['name']==state),'bufall'] = max_symb_size * values  / sqrt_max
        cutoff_value = min_symb_size/max_symb_size* sqrt_max
        sdf.loc[(sdf["sqrt"]<=cutoff_value)&(sdf['name']==state),'bufall'] = min_symb_size
        bvals = []
        value = sqrt_max * sqrt_max
        minvalue = cutoff_value*cutoff_value
        cnt = 0
        while cnt <2:
            for x in [10000000,1000000,100000,10000,1000]:
                if value>=x:
                    rp = round_place(value,x)
                    bvals.append(rp)
                    break

            value = value *.25
            if value<minvalue:
                break
            cnt +=1
        for x in [1000000,100000,10000,1000]:
            if cutoff_value*cutoff_value>=x:
                bvals.append(round_place(cutoff_value*cutoff_value,x))
                break

        bds = [max_symb_size * math.sqrt(x) / sqrt_max for x in bvals]

        for bd,x in zip(bds,bvals):
            legend_data["agency"].append("ALL")
            legend_data["value"].append("{:,}".format(x))
            geom = arcpy.PointGeometry(arcpy.Point(0,bd+max(bds)), sr).buffer(bd)
            legend_data["SHAPE"].append(geom)
            legend_data["state"].append(state)
        if "<=" not in legend_data["value"][-1]:
            legend_data["value"][-1] = "<=" + legend_data["value"][-1] 

    legend = pd.DataFrame.spatial.from_df(pd.DataFrame(legend_data),geometry_column="SHAPE")
    
    legend.spatial.to_featureclass(arcpy.env.workspace + "\\legend_data")
    arcpy.Delete_management(arcpy.env.workspace + "\\flma_points_visits")
    
    sdf.spatial.to_featureclass(arcpy.env.workspace + "\\flma_points_visits_final")
    arcpy.analysis.Buffer(arcpy.env.workspace + "\\flma_points_visits_final",arcpy.env.workspace + "\\flma_points_prj_charge_stations_buff","bufall")
    arcpy.analysis.Buffer(arcpy.env.workspace + "\\flma_points_visits_final",arcpy.env.workspace + "\\flma_points_prj_charge_stations_buffagen","bufbyagency")


def project_envelopes(pth,newname):
    project(pth,newname,4326)
    res = arcpy.FeatureEnvelopeToPolygon_management(newname,newname+"_ENV").getOutput(0)
    return res
    

def set_values(agency,lblfield,lblorgcode,sdf_env,flma_pnts,check_id=True):
    print("set values")
    for org in flma_pnts[flma_pnts["agency"]==agency]["org"].unique():
        sdf_temp = flma_pnts.loc[flma_pnts["org"]==org].copy().reset_index()
        pntcnt = len(sdf_temp)
        name = ""
        poss_name = sdf_temp["name"].unique()[0]
        for idx,row in sdf_env.spatial.relationship(sdf_temp,"contains").iterrows():
            selcnt = len(sdf_temp[sdf_temp["SHAPE"].geom.within(row["SHAPE"])])
            if selcnt == pntcnt:
                if check_id == True:
                    if row[lblorgcode] == org:
                        name = row[lblfield]
                        break
                    else:
                        poss_name = row[lblfield]
                else:
                    name = row[lblfield]
                    break
        if name == "":
            name = poss_name
        flma_pnts.loc[flma_pnts["org"]==org,"label"] = name

def get_labeling():
    flma_pnts = pd.DataFrame.spatial.from_featureclass(arcpy.env.workspace + "\\flma_points")

    #BLM
    print("Getting BLM Labels")
    envelope = project_envelopes(os.path.join(DOI_BLM_BND_GDB,r"blm_natl_admu_poly_webpub\blm_natl_admu_field_poly_webpub"),"blm_boundaries")
    sdf_env = pd.DataFrame.spatial.from_featureclass(envelope)
    lbl_field = "ADMU_NAME"
    code_field = "ADM_UNIT_CD"
    set_values("doi_blm",lbl_field,code_field,sdf_env,flma_pnts)
    
    
    #BOR
    print("Getting BOR Labels")
    envelope = project_envelopes(DOI_BR_BND,"bor_boundaries")
    sdf_env = pd.DataFrame.spatial.from_featureclass(envelope)
    lbl_field = "AreaOffice"
    code_field = "AreaOffice"
    set_values("doi_bor",lbl_field,code_field,sdf_env,flma_pnts,False)

    #FWS
    print("Getting FWS Labels")
    envelope = project_envelopes(DOI_FWS_BND,"fws_boundaries")
    sdf_env = pd.DataFrame.spatial.from_featureclass(envelope)
    lbl_field = "ORGNAME"
    code_field = "CostCenter"
    set_values("doi_fws",lbl_field,code_field,sdf_env,flma_pnts,False)


    #NPS
    print("Getting NPS Labels")
    envelope = project_envelopes(DOI_NPS_BND,"nps_boundaries")
    sdf_env = pd.DataFrame.spatial.from_featureclass(envelope)
    lbl_field = "UNIT_NAME"
    code_field = "CostCenter"
    set_values("doi_nps",lbl_field,code_field,sdf_env,flma_pnts,False)

    #FS
    print("Getting FS Labels")
    envelope = project_envelopes(USDA_FS_BND,"usfs_boundaries")
    sdf_env = pd.DataFrame.spatial.from_featureclass(envelope)
    lbl_field = "FORESTNAME"
    code_field = "CostCenter"
    set_values("usda_fs",lbl_field,code_field,sdf_env,flma_pnts,False)


    flma_pnts.spatial.to_featureclass(arcpy.env.workspace + "\\flma_points")

    fieldmappings = arcpy.FieldMappings()
    fieldmappings.addTable("flma_points")
    fieldmappings2 = arcpy.FieldMappings()
    fieldmappings2.addTable(STATE_BOUNDARIES)

    fieldmappings.addFieldMap(fieldmappings2.getFieldMap(fieldmappings2.findFieldMapIndex("postal")))
    flma_pnts = arcpy.SpatialJoin_analysis("flma_points",STATE_BOUNDARIES,"flma_points_states","JOIN_ONE_TO_ONE","KEEP_ALL",fieldmappings,match_option="CLOSEST").getOutput(0)
    return flma_pnts


def round_place(n,place=1000):
    return int(math.ceil(n / place)) * place

def project(in_fc:str,out_fc:str,epsg=102003):
    """Project the feature class to a new one
    Parameters:
        in_fc (str): input feature class path
        out_fc (str):  output feature class path
        epsg (int): ouput spatial reference system
    Returns:
        path to new projected fc"""

    pth = arcpy.management.Project(in_fc, out_fc, epsg).getOutput(0)
    return pth


def get_medoid(vX:np.array):
    try:
        vMean = np.mean(vX, axis=0)
        medoid = vX[np.argmin([sum((x - vMean)**2) for x in vX])]
        return medoid
    except:
        return [np.nan,np.nan]


def check_ct(l:str):
    """Checks for the type of charging station
    Parameters:
        string input
    Returns:
        DCFC or Level 2
    """
    if "J1772COMBO" in l:
        if "CHADEMO" in l:
            return "DCFC"
        elif "J1772" in l:
            return "Level 2"
    elif "J1772" in l:
        return "Level 2"
    return None

def create_point_feature_class(name:str):
    """ Create the template point feature class of a given name within the current workspace. WGS 1984 projection
    Parameters:
        name (str): Name of the output feature class.
    Returns:
        True if successful
        False if an error occured
    """
    try:
        arcpy.CreateFeatureclass_management(
                arcpy.env.workspace,
                name,
                "POINT",
                "", "DISABLED", "DISABLED",
                arcpy.SpatialReference("WGS 1984")
                )

        arcpy.AddField_management(name, "org", "TEXT", "", "", 100)
        arcpy.AddField_management(name, "name", "TEXT", "", "", 250)  # width driven by largest across layers
        arcpy.AddField_management(name, "matchid", "TEXT", "", "", 250)
        arcpy.AddField_management(name, "agency", "TEXT", "", "", 10)
        arcpy.AddField_management(name, "label", "TEXT", "", "", 250)
        #arcpy.AddField_management(name, "state", "TEXT", "", "", 10)
        return True
    except:
        return False
# --------------------------------------------------------------------------------------------------

def project_and_create_centroid_fc(polygon_fc:str,outpoint_fc:str):
    """
    Description:
        Project Polygon Feature Class and convert to point featureclass. Workspace should be set beforehand.
    Parameters:
        polygon_fc (str): path to polygon feature class
        outpoint_fc (str): name of output centroid featureclass
    Returns:
        Path to projected centroids
    """

    pth = arcpy.FeatureToPoint_management(polygon_fc, "in_memory//centroid_fc", "INSIDE").getOutput(0)
    res = project(pth,outpoint_fc)
    arcpy.Delete_management("in_memory//centroid_fc")
    return res


def main():



    gdb_name = 'ev_gap_analysis.gdb'
    fp_to_gdb = os.path.join(OUTPUT_DIR, gdb_name)

    if arcpy.Exists(fp_to_gdb):
        arcpy.Delete_management(fp_to_gdb)
        time.sleep(1)

    arcpy.CreateFileGDB_management(OUTPUT_DIR, gdb_name)
    arcpy.env.workspace = fp_to_gdb
    parts = []
    centroids = []
    parts,centroids = import_doi_br(parts,centroids)
    parts,centroids = import_doi_blm(parts,centroids)
    
    parts,centroids = import_doi_fws(parts,centroids)
    parts,centroids = import_doi_nps(parts,centroids)
    parts,centroids = import_usda_fs(parts,centroids)
    parts = import_usace(parts)
    merge(parts,2500)
    #match_visitation_data()
    pth = match_visitation_data_with_centroids(centroids)
    buffer_by_agency_state(pth)
    
    get_labeling()
    #
    #Get an API KEY from : https://developer.nrel.gov/signup/
    key_str = ""
    if key_str == "":
        import_ev_stations()
    else:
        print("\nUsing NREL api to get EV stations...")
        import_ev_stations_from_api(key_str)


# --------------------------------------------------------------------------------------------------
if __name__ == "__main__":


    main()

    print('\ndone')





