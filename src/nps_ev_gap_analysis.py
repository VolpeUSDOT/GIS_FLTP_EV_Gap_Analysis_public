import os
import arcpy
import datetime

"""
This script takes origin-destination pairs and splits segments of road by EV stations that are 
within 3 miles of the road. Using the GAP_SIZE threshold, these segments are identified as areas 
with fuel station gaps and are subset into a different feature class. Gaps that are shared by 
multiple routes are also subset into separate feature class, and are referred to as "shared gaps". 

Create Date: 6/20/2019
Created By: The Volpe Center's Geospatial Team 
"""


# UPDATE USER ARGUMENTS
# ---------------------
PROGRAM_DIR = r"H:\nps_ev\program"
ALT_FUEL_CSV_NAME = "alt_fuel_stations (Jun 14 2019).csv"
ORIGINS_XLS_NAME = "origins.xls"
DESTINATION_XLS_NAME = "destinations.xls"
GAP_SIZE = 50  # Minimum distance between stations in Miles
MAX_STATION_DIST_FROM_ROAD = 3  # Max distance station can be from road in Miles

# ---------------------------------------------------------------------------


def load_and_update_table(csv):

    print("Table to Table Conversion")
    arcpy.TableToTable_conversion(csv, PROGRAM_GDB, "alt_fuel_stations")

    arcpy.AddField_management("alt_fuel_stations", "CONNECTOR_TYPE", "TEXT")

    arcpy.MakeTableView_management("alt_fuel_stations", "alt_fuel_stations_tbl_view")

    # CATEGORIZE CHADEMO STATIONS
    # ---------------------------
    print("Categorizing CHAdeMO Stations")
    arcpy.SelectLayerByAttribute_management("alt_fuel_stations_tbl_view", "NEW_SELECTION",
                                            "EV_Connector_Types = 'J1772 CHADEMO' "
                                            "OR EV_Connector_Types = 'J1772 CHADEMO TESLA' "
                                            "OR EV_Connector_Types = 'J1772 NEMA1450 CHADEMO' "
                                            "OR EV_Connector_Types = 'J1772 NEMA515 CHADEMO' "
                                            "OR EV_Connector_Types = 'NEMA520 CHADEMO J1772' "
                                            "OR EV_Connector_Types = 'NEMA520 J1772 CHADEMO' "
                                            "OR EV_Connector_Types = 'TESLA J1772 CHADEMO' "
                                            "OR EV_Connector_Types = 'J1772 CHADEMO NEMA520 TESLA'")

    arcpy.CalculateField_management("alt_fuel_stations_tbl_view", "CONNECTOR_TYPE", '"CHAdeMO"')

    # CATEGORIZE J1772 AND CHADEMO STATIONS
    # -------------------------------------
    print("Categorizing J1772 and CHAdeMO stations")
    arcpy.SelectLayerByAttribute_management("alt_fuel_stations_tbl_view", "NEW_SELECTION",
                                            "EV_Connector_Types = 'CHADEMO J1772' "
                                            "OR EV_Connector_Types = 'J1772 CHADEMO' "
                                            "OR EV_Connector_Types = 'J1772 CHADEMO TESLA' "
                                            "OR EV_Connector_Types = 'J1772 NEMA1450 CHADEMO' "
                                            "OR EV_Connector_Types = 'J1772 NEMA515 CHADEMO' "
                                            "OR EV_Connector_Types = 'NEMA520 CHADEMO J1772' "
                                            "OR EV_Connector_Types = 'NEMA520 J1772 CHADEMO' "
                                            "OR EV_Connector_Types = 'TESLA J1772 CHADEMO' "
                                            "OR EV_Connector_Types = 'J1772 CHADEMO NEMA520 TESLA'")

    arcpy.CalculateField_management("alt_fuel_stations_tbl_view", "CONNECTOR_TYPE",
                                    '"J1772 & CHAdeMO"')

    # CATEGORIZE CHADEMO AND SAE STATIONS
    # -----------------------------------
    print("Categorizing CHAdeMO and SAE Stations")
    arcpy.SelectLayerByAttribute_management("alt_fuel_stations_tbl_view", "NEW_SELECTION",
                                            "EV_Connector_Types = 'CHADEMO J1772COMBO' "
                                            "OR EV_Connector_Types = 'CHADEMO J1772COMBO TESLA' "
                                            "OR EV_Connector_Types = 'J1772COMBO CHADEMO' "
                                            "OR EV_Connector_Types = 'NEMA520 CHADEMO J1772COMBO' "
                                            "OR EV_Connector_Types = 'TESLA CHADEMO J1772COMBO'")

    arcpy.CalculateField_management("alt_fuel_stations_tbl_view", "CONNECTOR_TYPE",
                                    '"CHAdeMO & SAE"')

    # CATEGORIZE SAE STATIONS
    # -----------------------
    print("Categorizing SAE stations")
    arcpy.SelectLayerByAttribute_management("alt_fuel_stations_tbl_view", "NEW_SELECTION",
                                            "EV_Connector_Types = 'J1772COMBO'")
    arcpy.CalculateField_management("alt_fuel_stations_tbl_view", "CONNECTOR_TYPE", '"SAE"')

    # CATEGORIZE J1772 AND SAE STATIONS
    # ---------------------------------
    print("Categorizing J1772 and SAE Stations")
    arcpy.SelectLayerByAttribute_management("alt_fuel_stations_tbl_view", "NEW_SELECTION",
                                            "EV_Connector_Types = 'J1772 J1772COMBO' "
                                            "OR EV_Connector_Types = 'J1772COMBO J1772'")
    arcpy.CalculateField_management("alt_fuel_stations_tbl_view", "CONNECTOR_TYPE", '"J1772 & SAE"')

    # CATEGORIZE J1772, CHADEMO, AND SAE STATIONS
    # -------------------------------------------
    print("Categorizing J1772, CHAdeMO, and SAE Stations")
    arcpy.SelectLayerByAttribute_management("alt_fuel_stations_tbl_view", "NEW_SELECTION",
                                            "EV_Connector_Types = 'CHADEMO J1772 J1772COMBO' "
                                            "OR EV_Connector_Types = 'CHADEMO J1772COMBO J1772' "
                                            "OR EV_Connector_Types = 'CHADEMO J1772COMBO J1772 "
                                            "TESLA' "
                                            "OR EV_Connector_Types = 'J1772 CHADEMO J1772COMBO' "
                                            "OR EV_Connector_Types = 'J1772 J1772COMBO CHADEMO' "
                                            "OR EV_Connector_Types = 'J1772 TESLA CHADEMO "
                                            "J1772COMBO' "
                                            "OR EV_Connector_Types = 'J1772COMBO CHADEMO J1772'")

    arcpy.CalculateField_management("alt_fuel_stations_tbl_view", "CONNECTOR_TYPE",
                                    '"J1772 & CHAdeMO & SAE"')

    # CATEGORIZE J1772 STATIONS
    # -------------------------
    print("Categorizing J1772 Stations")
    arcpy.SelectLayerByAttribute_management("alt_fuel_stations_tbl_view", "NEW_SELECTION",
                                            "EV_Connector_Types = 'J1772' "
                                            "OR EV_Connector_Types = 'J1772 NEMA1450' "
                                            "OR EV_Connector_Types = 'J1772 NEMA1450 TESLA' "
                                            "OR EV_Connector_Types = 'J1772 NEMA515' "
                                            "OR EV_Connector_Types = 'J1772 NEMA520' "
                                            "OR EV_Connector_Types = 'J1772 TESLA' "
                                            "OR EV_Connector_Types = 'J1772 TESLA NEMA515' "
                                            "OR EV_Connector_Types = 'NEMA1450 J1772' "
                                            "OR EV_Connector_Types = 'NEMA515 J1772' "
                                            "OR EV_Connector_Types = 'NEMA520 J1772' "
                                            "OR EV_Connector_Types = 'NEMA520 J1772 TESLA' "
                                            "OR EV_Connector_Types = 'TESLA J1772'")

    arcpy.CalculateField_management("alt_fuel_stations_tbl_view", "CONNECTOR_TYPE", '"J1772"')

    # # CATEGORIZE TESLA STATIONS
    # # -------------------------
    # print("Categorizing Tesla Stations")
    # arcpy.SelectLayerByAttribute_management("alt_fuel_stations_tbl_view", "NEW_SELECTION",
    #                                         "EV_Connector_Types = 'TESLA' "
    #                                         "OR EV_Connector_Types = 'NEMA515 TESLA' "
    #                                         "OR EV_Connector_Types = 'TESLA NEMA1450 NEMA520' "
    #                                         "OR EV_Connector_Types = 'TESLA NEMA520'")

    # DELETE UNUSED STATION TYPES
    # ---------------------------
    arcpy.SelectLayerByAttribute_management("alt_fuel_stations_tbl_view", "NEW_SELECTION",
                                            "CONNECTOR_TYPE IS NULL")
    arcpy.DeleteRows_management("alt_fuel_stations_tbl_view")


def create_routes(origins_tbl, destinations_tbl, output_fc):

    # CREATE ROUTE LAYER AND SOLVE ROUTES
    # ------------------------------------
    print("Create route layer and solve routes")

    route_lyr = arcpy.MakeRouteAnalysisLayer_na("https://www.arcgis.com/", "route_lyr",
                                                "Driving Time", "USE_CURRENT_ORDER", "6/8/2019",
                                                "LOCAL_TIME_AT_LOCATIONS", "ALONG_NETWORK",
                                                None)

    origins_lyr = arcpy.MakeXYEventLayer_management(origins_tbl, "LON", "LAT", "origins_lyr")
    destinations_lyr = arcpy.MakeXYEventLayer_management(destinations_tbl, "LON", "LAT",
                                                         "destinations_lyr")

    fields = ["OBJECTID", "NAME"]
    origin_cursor = arcpy.da.SearchCursor(origins_lyr, fields)
    destination_cursor = arcpy.da.SearchCursor(destinations_lyr, fields)
    for origin_row in origin_cursor:
        origin_oid, origin_name = origin_row
        arcpy.SelectLayerByAttribute_management(origins_lyr, "NEW_SELECTION",
                                                "OBJECTID = {}".format(origin_oid))

        destination_cursor.reset()

        for destination_row in destination_cursor:
            destination_oid, destination_name = destination_row
            # if ((origin_name == "Salt Lake City" and destination_name == "Arches")
            #     or (origin_name == "Phoenix" and destination_name == "Arches")
            #     or (origin_name == "Las Vegas" and destination_name == "Bryce Canyon")
            #     or (origin_name == "Denver" and destination_name == "Bryce Canyon")
            #    ):
            #     continue

            print("Routing {} to {}".format(origin_name, destination_name))

            # CLEAR STOPS SUBLAYER AND ADD ORIGIN
            # -----------------------------------
            arcpy.AddLocations_na(route_lyr, "Stops", origins_lyr, "", "", "", "", "", "CLEAR")

            # APPEND DESTINATION TO STOPS SUBLAYER
            # ------------------------------------
            arcpy.SelectLayerByAttribute_management(destinations_lyr, "NEW_SELECTION",
                                                    "OBJECTID = {}".format(destination_oid))
            arcpy.AddLocations_na(route_lyr, "Stops", destinations_lyr, "", "", "", "", "",
                                  "APPEND")

            # MAKE ROUTE AND APPEND TO ROUTE FC
            # ---------------------------------
            arcpy.Solve_na(route_lyr)

            arcpy.CalculateField_management("route_lyr\Routes", "Name",
                                            "'{}' + ' - ' + '{}'".format(origin_name,
                                                                         destination_name))

            origin_name_no_spaces = origin_name.replace(" ", "_")
            destination_name_no_spaces = destination_name.replace(" ", "_")
            arcpy.CopyFeatures_management("route_lyr\Routes",
                                          os.path.join(ROUTES_GDB,
                                                       "{}_to_{}".format(origin_name_no_spaces,
                                                                         destination_name_no_spaces)
                                                       ))

    return output_fc


def segment_routes():

    print("Running function: segment_routes ...")

    # MAKE BASELINE STATION LAYER
    # ---------------------------
    print("Make baseline station layer")
    arcpy.MakeTableView_management("alt_fuel_stations", "alt_fuel_stations_view",
                                   "CONNECTOR_TYPE = 'CHAdeMO & SAE' "
                                   "Or CONNECTOR_TYPE = 'J1772' "
                                   "Or CONNECTOR_TYPE = 'J1772 & CHAdeMO' "
                                   "Or CONNECTOR_TYPE = 'J1772 & CHAdeMO & SAE' "
                                   "Or CONNECTOR_TYPE = 'J1772 & SAE'")
    arcpy.MakeXYEventLayer_management("alt_fuel_stations_view", "Longitude", "Latitude",
                                      "stations_layer")
    #                                   , "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',\
    # SPHEROID['WGS_1984',6378137.0,298.257223563]],\
    # PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]];\
    # -400 -400 1000000000;-100000 10000;-100000 10000;8.98315284119522E-09;0.001;0.001;\
    # IsHighPrecision", None)

    # SPLIT ROUTE BY STATIONS WITHIN 3 MILES TO CREATE ROUTE SEGMENTS
    # ---------------------------------------------------------------
    print("Segment route by stations...")
    arcpy.LocateFeaturesAlongRoutes_lr("stations_layer", "in_memory\\route", "NAME",
                                       "{} Miles".format(MAX_STATION_DIST_FROM_ROAD),
                                       "in_memory\\located_stations", "RID Point METERS",
                                       "FIRST", "NO_DISTANCE", "ZERO", "FIELDS", "M_DIRECTON")

    arcpy.DissolveRouteEvents_lr("in_memory\\located_stations", "RID Point METERS",
                                 "Station_Name;Street_Address;ZIP;CONNECTOR_TYPE",
                                 "in_memory\\located_stations_dissolve",
                                 "RID Point METERS", "DISSOLVE", "INDEX")

    arcpy.MakeRouteEventLayer_lr("in_memory\\route", "NAME", "in_memory\\located_stations_dissolve",
                                 "RID Point METERS", "located_stations_events", None,
                                 "NO_ERROR_FIELD", "NO_ANGLE_FIELD", "NORMAL", "ANGLE", "LEFT",
                                 "POINT")

    arcpy.SplitLineAtPoint_management("in_memory\\route", "located_stations_events",
                                      "in_memory\\route_split_by_stations", "1 Meters")

    arcpy.SpatialJoin_analysis("in_memory\\route_split_by_stations", "located_stations_events",
                               "in_memory\\route_with_station_data", "JOIN_ONE_TO_ONE", "KEEP_ALL",
                               r'NAME "NAME" true true false 260 Text 0 0,First,#,'
                               r'route_split_by_stations,NAME,0,255;'
                               r'METERS "METERS" true true false 8 Double 0 0,First,#,'
                               r'located_stations_events,METERS,-1,-1;'
                               r'Station_Name "Station Name" true true false 30000 Text 0 0,'
                               r'Join,"/",located_stations_events,Station_Name,0,255;'
                               r'Street_Address "Street Address" true true false 30000 Text 0 0,'
                               r'Join,"/",located_stations_events,Street_Address,0,255;'
                               r'ZIP "ZIP" true true false 1000 Text 0 0,Join,"/",'
                               r'located_stations_events,ZIP,0,255;'
                               r'CONNECTOR_TYPE "CONNECTOR_TYPE" true true false 255 Text 0 0,'
                               r'Join,"/",located_stations_events,CONNECTOR_TYPE,0,255',
                               "INTERSECT", "1 Meters", None)

    arcpy.DeleteField_management("in_memory\\route_with_station_data", "Join_Count;TARGET_FID")

    arcpy.Append_management("in_memory\\route_with_station_data", "all_routes_with_station_data",
                            "TEST", None, None)


def find_gaps():

    print('\n')
    print("Running function: find_gaps ...")

    # CALCULATE SEGMENTS TO FIND GAPS
    # -------------------------------
    print("Calculate segments...")
    arcpy.AddField_management("all_routes_with_station_data", "FULL_GAP_LENGTH", "DOUBLE")
    arcpy.CalculateField_management("all_routes_with_station_data", "FULL_GAP_LENGTH",
                                    "!shape.length@miles!")

    arcpy.MakeFeatureLayer_management("all_routes_with_station_data",
                                      "all_routes_with_station_data_lyr")
    arcpy.SelectLayerByAttribute_management("all_routes_with_station_data_lyr", "NEW_SELECTION",
                                            "FULL_GAP_LENGTH >= " + str(GAP_SIZE), None)
    arcpy.CopyFeatures_management("all_routes_with_station_data_lyr", "routes_with_gaps")
    arcpy.SelectLayerByAttribute_management("all_routes_with_station_data_lyr", "CLEAR_SELECTION")

    # FIND SHARED GAPS - GAPS THAT ARE USED BY MULTIPLE ROUTES
    # --------------------------------------------------------
    print("Find shared gaps")
    arcpy.Intersect_analysis("routes_with_gaps #", "shared_gap_segments", "ALL", None, "INPUT")
    arcpy.AddField_management("shared_gap_segments", "SHARED_GAP_LENGTH", "DOUBLE")


if __name__ == "__main__":

    start_time = datetime.datetime.now()

    PROGRAM_GDB = os.path.join(PROGRAM_DIR, "program.gdb")
    if not arcpy.Exists(PROGRAM_GDB):
        arcpy.CreateFileGDB_management(PROGRAM_DIR, "program.gdb")

    arcpy.env.workspace = PROGRAM_GDB
    arcpy.env.overwriteOutput = True

    sr = arcpy.SpatialReference("WGS 1984")
    arcpy.env.outputCoordinateSystem = sr

    # STEP 1 - CATEGORIZE STATIONS
    # ----------------------------
    print("STEP 1: CATEGORIZE STATIONS")

    alt_fuel_station_csv = os.path.join(PROGRAM_DIR, ALT_FUEL_CSV_NAME)
    load_and_update_table(alt_fuel_station_csv)

    # STEP 2 - CREATE ROUTES
    # ----------------------
    print("STEP 2: CREATE ROUTES")

    origins = os.path.join(PROGRAM_DIR, ORIGINS_XLS_NAME)
    destinations = os.path.join(PROGRAM_DIR, DESTINATION_XLS_NAME)

    ROUTES_GDB = os.path.join(PROGRAM_DIR, "routes.gdb")
    if not arcpy.Exists(ROUTES_GDB):
        arcpy.CreateFileGDB_management(PROGRAM_DIR, "routes.gdb")

    # Importing Origins and Destinations
    print("Importing Origins and Destinations")
    arcpy.ExcelToTable_conversion(origins, "origins")
    arcpy.ExcelToTable_conversion(destinations, "destinations")

    # Create Route Output FC
    print("Creating Route Feature Class")
    all_routes = arcpy.CreateFeatureclass_management(PROGRAM_GDB, "routes", "POLYLINE", "", "", "",
                                                     sr)
    arcpy.AddFields_management(all_routes, [["Name", "TEXT"], ["Total_TravelTime", "DOUBLE"],
                                            ["Total_Miles", "DOUBLE"]])

    create_routes("origins", "destinations", all_routes)

    # STEP 3 - ANALYSIS
    # -----------------
    print("STEP 3: ANALYSIS")

    # Change GDB workspace just to collect route list
    arcpy.env.workspace = ROUTES_GDB
    ROUTE_GDB_LIST = arcpy.ListFeatureClasses("*")

    # Reset workspace to the main program GDB
    arcpy.env.workspace = PROGRAM_GDB

    arcpy.CreateFeatureclass_management(PROGRAM_GDB, "all_routes_with_station_data", "POLYLINE")
    arcpy.AddFields_management("all_routes_with_station_data",
                               [["NAME", "TEXT"],
                                ["METERS", "DOUBLE"],
                                ["Station_Name", "TEXT", "", 30000],
                                ["Street_Address", "TEXT", "", 30000],
                                ["ZIP", "TEXT", "", 1000],
                                ["CONNECTOR_TYPE", "TEXT", "", 1000]
                                ])

    for route_name in ROUTE_GDB_LIST:
        print("\n")
        print("Current Route:", route_name)

        arcpy.CreateRoutes_lr(os.path.join(ROUTES_GDB, route_name), "Name", "in_memory\\route",
                              "LENGTH")

        segment_routes()

    find_gaps()

    print("H:M:S " + str(datetime.datetime.now() - start_time))
