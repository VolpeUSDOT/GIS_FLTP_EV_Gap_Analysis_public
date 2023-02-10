import os, sys
from pathlib import Path

import pandas as pd
from arcgis.features import GeoAccessor, GeoSeriesAccessor

import arcpy

import numpy as np
from sklearn.neighbors import KDTree


"""
This script finds charging stations within 50 miles of the parking lot. 

Create Date: 2022
Created By: The Volpe Center's Geospatial Team 
"""
# --------------------------------------------------------------------------------------------------

BASE_DIR = Path(__file__).parents[1].resolve()

INPUT_DIR  = os.path.join(BASE_DIR, "data\\input")
OUTPUT_DIR = os.path.join(BASE_DIR, "data\\output")

DISTANCE_TO_DCFC = 80467.2 #50 miles in meters
DISTANCE_TO_LEVEL2 = 16093.4 #10 miles in meters


def project(in_fc,out_fc,epsg=102003):
    """Project the feature class to a new one
    Parameters:
        in_fc (str): input feature class path
        out_fc (str):  output feature class path
        epsg (int): ouput spatial reference system
    Returns:
        Nothing"""

    arcpy.management.Project(in_fc, out_fc, epsg) 


def find_charging_stations(gdb_path,distance_metric="manhattan"):
    """Calculates the number of charging stations within specified distances
    Parameters:
        gdb_path (str): path to geodatabase where flma and charging station point feature classes are located
    Returns:
        (str) path to new feature class storing the flma data"""

    cs_df = pd.DataFrame.spatial.from_featureclass(gdb_path + "\\charging_stations_prj")
    flma_df = pd.DataFrame.spatial.from_featureclass(gdb_path + "\\flma_points_prj")
    flma_df["DCFC_CNT"]=0
    flma_df["LEV2_CNT"]=0
    flma_df["XC"] = flma_df["SHAPE"].apply(lambda s: s.x)
    flma_df["YC"] = flma_df["SHAPE"].apply(lambda s: s.y)
    cs_df["XC"] = cs_df["SHAPE"].apply(lambda s: s.x)
    cs_df["YC"] = cs_df["SHAPE"].apply(lambda s: s.y)

    for cs,r,fld in [("DCFC",DISTANCE_TO_DCFC,"DCFC_CNT"),("Level 2",DISTANCE_TO_LEVEL2,"LEV2_CNT")]:
        cs_tree = KDTree(cs_df[cs_df["ev_connector_types_lev"]==cs][["XC","YC"]].values, metric=distance_metric)
        res = cs_tree.query_radius(flma_df[["XC","YC"]].values,r=r,return_distance=True)
        flma_df[fld]= [len(x) for x in res[0]]

    flma_df["ACCESS"] = flma_df["DCFC_CNT"] + flma_df["LEV2_CNT"]
    flma_df["CLASS"] = ""
    flma_df["aclass"] = ""
    flma_df.loc[(flma_df["DCFC_CNT"]==0)&(flma_df["LEV2_CNT"]==0),"CLASS"] = "No Charging Stations Within Ranges"
    flma_df.loc[(flma_df["DCFC_CNT"]>0)&(flma_df["LEV2_CNT"]==0),"CLASS"] = "DCFC Charging Stations Within Range"
    flma_df.loc[(flma_df["DCFC_CNT"]>0)&(flma_df["LEV2_CNT"]>0),"CLASS"] = "DCFC and Level 2 Charging Stations Within Ranges"
    flma_df.loc[(flma_df["DCFC_CNT"]==0)&(flma_df["LEV2_CNT"]>0),"CLASS"] = "Level 2 Charging Stations Within Range"
    flma_df.loc[(flma_df["CLASS"] == "No Charging Stations Within Ranges"),"aclass"]=1
    flma_df.loc[(flma_df["CLASS"] == "DCFC Charging Stations Within Range"),"aclass"]=2
    flma_df.loc[(flma_df["CLASS"] == "DCFC and Level 2 Charging Stations Within Ranges"),"aclass"]=3
    flma_df.loc[(flma_df["CLASS"] == "Level 2 Charging Stations Within Range"),"aclass"]=4
    flma_df.spatial.to_featureclass(gdb_path + "\\flma_points_prj_charge_stations")
    
    return gdb_path + "\\flma_points_prj_charge_stations"

def main():


    gdb_name = 'ev_gap_analysis.gdb'
    fp_to_gdb = os.path.join(OUTPUT_DIR, gdb_name)
    print(fp_to_gdb)
   

    if arcpy.Exists(fp_to_gdb) == False:
        raise Exception("Missing geodatabase....")

    arcpy.env.workspace = fp_to_gdb
    arcpy.env.overwriteOutput = True
    print("Project feature classes...")

    try:
        project("flma_points_states","flma_points_prj")
    except:
        print("flma_points_prj already exists...")
    
    try:
        project("charging_stations","charging_stations_prj")
    except:
        print("charging_stations_prj alreayd exists...")
    


    print("Find closest Charging Stations....)")
    out_fc = find_charging_stations(fp_to_gdb)
    #print("Convert to Raster Layer....")
    #arcpy.conversion.PointToRaster(out_fc,"ACCESS","cs_access",cell_assignment="MEAN",cellsize=4000)

# --------------------------------------------------------------------------------------------------
if __name__ == "__main__":


    main()

    print('\ndone')