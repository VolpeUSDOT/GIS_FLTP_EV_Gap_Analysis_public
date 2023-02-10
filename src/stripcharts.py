import stripcharts_data
from pathlib import Path
import arcpy
import requests
import plot_tools



"""
Automates the process of creating the strip charts. 

Create Date: 2022
Created By: The Volpe Center's Geospatial Team 
"""

DISSOLVE_FIELDS = None#"AADT;Route_ID;URBAN_CODE"
CHARGING_STATIONS = stripcharts_data.OUTPUTS / "ev_gap_analysis.gdb" / "charging_stations_prj"
PLACES = stripcharts_data.INPUTS / "us_places.shp"
#templates
MAP_PROJ = stripcharts_data.INPUTS / "templates"  / "EmptyProj.aprx"
LAYOUT_TEMP = stripcharts_data.INPUTS / "templates"  / "LayoutTemplate.pagx"
STATIONS_LAYER = stripcharts_data.INPUTS / "templates"  / "StationsLayer.lyrx"
CENTERLINE_TEMP = stripcharts_data.INPUTS / "templates"  / "centerline2.lyrx"

#INPUT INFORMATION
ATTRIBUTE_DATA =[{"field_name":"AADT",
                  "label":"AADT",
                  "color":"#1E88E5"},
                #   {"field_name":"VMT",
                #   "label":"VMT",
                #   "color":"#E69F00"},
                {"field_name":"ELEVFT",
                  "label":"Elevation\n(Feet)",
                  "color":"#004D40"},
                {"field_name":"URBANSTR",
                "label":"Urban Areas",
                  "color":"#E69F00"}
                  
 ]


NEAR_POINT_DATA=[{"source":CHARGING_STATIONS,"label":"Charging Station","attribute":"ev_connector_types_lev",
                "color":{"DCFC":"#994455","Level 2":"#AAAA00"},"layer_file":stripcharts_data.INPUTS/"templates"/"Charging_Station.lyrx","search_distance":26400}]

POLYGON_DATA = [{"source":PLACES,"label":"Populated Places","attribute":"NAME","color":"#828282","layer_file":stripcharts_data.INPUTS/"templates"/"us_places.lyrx"}]

STATION_LENGTHS =  8000 #feet
STATION_SPACING = 20 #miles 100000
FILL_SPACE = True

NEW_FGDB = "strip_chart_data.gdb"
DEM_FOLDER = stripcharts_data.OUTPUTS / "dem_images"
PROJECTS_FOLDER = stripcharts_data.OUTPUTS / "strip_chart_projs"
MAP_NAME = "Basemap"
DEM_FOLDER.mkdir(exist_ok=True)
PROJECTS_FOLDER.mkdir(exist_ok=True)

#API Data for Elevations
ELEVATION_DATASET_NAME = "National%20Elevation%20Dataset%20(NED)%201/3%20arc-second"
TNM_API_URL = "https://tnmaccess.nationalmap.gov/api/v1/products?"

def create_fgdb(dir:Path, fgdb:str,delete=False)->Path:
    """
    Creates a file geodatabase in the supplied directory.
    Deletes if existing.
    Args:
        dir (Path): path to workspace folder
        fgdb (str): name with gdb extension
    Returns:
        Path: to file geodatabase
    """
    if delete:
        if arcpy.Exists(str(dir/fgdb)):
            arcpy.Delete_management(str(dir/fgdb))
        res = arcpy.CreateFileGDB_management(str(dir),fgdb).getOutput(0)
    else:
        if arcpy.Exists(str(dir/fgdb)):
            return dir/fgdb
        res = arcpy.CreateFileGDB_management(str(dir),fgdb).getOutput(0)
    return Path(res)


def create_index(workspace:Path,delete=False)->Path:
    """
    Creates a file geodatabase in the supplied directory.
    Deletes if existing.
    Args:
        dir (Path): path to workspace folder
        fgdb (str): name with gdb extension
    Returns:
        Path: to file geodatabase
    """
    nm = workspace/"index_routes"
    if delete:
        if arcpy.Exists(str(nm)):
            arcpy.Delete_management(str(nm ))
    else:
        if arcpy.Exists(str(nm)):
            return nm

    res = arcpy.CreateFeatureclass_management(str(workspace),"index_routes","POLYGON",spatial_reference=arcpy.SpatialReference(4326)).getOutput(0)
    arcpy.AddField_management(res,"NAME","TEXT")
    arcpy.AddField_management(res,"LABEL","TEXT")
    return Path(res)

def extract_and_project(fc:Path,wc:str,epsg:int,outfile:Path)->Path:
    """
    Selects and projects.
    Args:
        fc (Path): path to feature class
        wc (str): where clause
        epsg (int): coordinate reference system code
        outfile (Path): path and name of the output projected feature class
    Returns:
        Path: to projected feature class
    """
    outsr = arcpy.SpatialReference(epsg)
    print(outfile)
    selection = arcpy.SelectLayerByAttribute_management(str(fc),"NEW_SELECTION",wc).getOutput(0)
    projected = arcpy.Project_management(selection,str(outfile),outsr).getOutput(0)
    return Path(projected)

def extent_as_string_latlong_url(fc:Path)->str:
    wgs84 = arcpy.SpatialReference(4326)
    desc = arcpy.Describe(str(fc))
    ext = desc.extent
    ll = arcpy.PointGeometry(arcpy.Point(ext.XMin,ext.YMin),desc.spatialReference)
    ur = arcpy.PointGeometry(arcpy.Point(ext.XMax,ext.YMax),desc.spatialReference)
    llprj = ll.projectAs(wgs84)
    urprj = ur.projectAs(wgs84)
    return "{},{},{},{}".format(llprj.firstPoint.X,llprj.firstPoint.Y,urprj.firstPoint.X,urprj.firstPoint.Y)

def extent_as_string(fc:Path)->str:
    wgs84 = arcpy.SpatialReference(4326)
    desc = arcpy.Describe(str(fc))
    rect = "{} {} {} {}".format(desc.extent.XMin,desc.extent.YMin,desc.extent.XMax,desc.extent.YMax)
    return rect

def convert_urban(fc:Path):
    try:
        arcpy.AddField_management(str(fc),"URBANSTR","TEXT",field_length="50")
    except:
        pass
    with arcpy.da.UpdateCursor(str(fc),["URBAN_CODE","URBANSTR"]) as uc:
        for row in uc:
            if row[0] == 99999:
                row[1] = "Rural"
            else:
                row[1] = "Urban"
            uc.updateRow(row)

def add_vmt(fc:Path):
    units = arcpy.Describe(str(fc)).spatialReference.linearUnitName
    unitconv = 1
    if units == "Meter":
        unitconv = 0.000621371
    elif units == "Foot_US":
        unitconv = 0.000189394
    try:
        arcpy.AddField_management(str(fc),"VMT","DOUBLE")
    except:
        pass
    with arcpy.da.UpdateCursor(str(fc),["SHAPE@","AADT","VMT"]) as uc:
        for row in uc:
            if row[0]:
                miles = row[0].length * unitconv
                row[2] = miles * row[1]
            else:
                row[2] = 0
            uc.updateRow(row)

def get_elevation_data(base_url:str,bbox:str,datasets:str,outdir:Path,timeout=60):
    filelist = []
    api_url = "{}datasets={}&bbox={}".format(base_url,datasets,bbox)
    response = requests.get(api_url,timeout=timeout)
    result = response.json()
    if 'message' in result:
        if result['message'] == 'Endpoint request timed out':
            print("timed out, trying again...")
            response = requests.get(api_url,timeout=timeout+timeout)
            result = response.json()
        #raise Exception("Error getting the DEM from the api.")
    if 'total' not in result:
        return None
        #raise Exception("Error getting the DEM from the api.")
    print("Total DEM Files that may need to be downloaded: {}".format(result["total"]))

    for item in result["items"]:
        downloadurl = item["downloadURL"]
        name = downloadurl.split("/")[-1]
        if ".tif" in name:
            outfile = outdir / name
            if outfile.exists() == False:
                imgresponse = requests.get(downloadurl)
                with open(outfile,"wb") as file:
                    file.write(imgresponse.content)
            filelist.append(str(outfile))
        else:
            print("Results are tif files")
            return filelist
    return filelist

def mosaic_and_append(rasterFiles:list,outws:Path,mosaic_name:str,linefc:Path):
    """
    mosaics raster dataset, clips to extent, and appends to line feature class

    Args:
        rasterFiles (list): list of raster paths
        outws (Path): output workspace

    """
    ext = "Spatial"
    if arcpy.CheckExtension("Spatial") == "Available":
        arcpy.CheckOutExtension("Spatial")
    else:
        if arcpy.CheckExtension("3D") == "Available":
            arcpy.CheckOutExtension("3D")
            ext = "3D"
        else:
            print("No spatial or 3d analyst extension available")
            raise Exception()


    wstype = arcpy.Describe(str(outws)).workspaceType
    rect = extent_as_string(linefc)
    desc = arcpy.Describe(str(linefc))
    clip_name = mosaic_name + "clp"
    rect = "{} {} {} {}".format(desc.extent.XMin,desc.extent.YMin,desc.extent.XMax,desc.extent.YMax)
    if wstype == "FileSystem":
        clip_name+=".tif"
        mosaic_name+=".tif"
    clipout = outws / clip_name
    if arcpy.Exists(str(clipout)):
        arcpy.Delete_management(str(clipout))
    with arcpy.EnvManager(scratchWorkspace=str(outws), workspace=str(outws)):
        print("Mosaicing Rasters")
        res = arcpy.management.MosaicToNewRaster(rasterFiles, str(outws), mosaic_name , desc.spatialReference, "32_BIT_FLOAT", None, 1, "LAST", "FIRST").getOutput(0)
        print("Clipping to study area Rasters")
        cres = arcpy.management.Clip(res, rect, out_raster=str(clipout), in_template_dataset=str(linefc), clipping_geometry= "NONE", maintain_clipping_extent="NO_MAINTAIN_EXTENT").getOutput(0)
        print("Adding Surface Information.")
        if ext == "Spatial":
            arcpy.sa.AddSurfaceInformation(str(linefc), cres, "Z_MEAN")
        else:
            arcpy.ddd.AddSurfaceInformation(str(linefc), cres, "Z_MEAN")
        try:
            arcpy.AddField_management(str(linefc),"ELEVFT","DOUBLE")
        except:
            pass
        with arcpy.da.UpdateCursor(str(linefc),["Z_Mean","ELEVFT"]) as uc:
            for row in uc:
                if row[0]:
                    row[1] = row[0]*3.28084
                uc.updateRow(row)

def main(centerlines):
    index = [("State_Code","hpms_state_code"),("ROUTE_NUMB","hpms_route_num"),("ROUTE_NAME","hpms_route_name"),("Route_ID","hpms_route_id")]
    ws = create_fgdb(stripcharts_data.OUTPUTS,NEW_FGDB)
    working_folder = PROJECTS_FOLDER / "scratch_folder_for_charts"
    working_folder.mkdir(parents=True,exist_ok=True,)
    index_map = create_index(ws,True)
    for d in centerlines:
            
        for f,indnm in index:
            try:
                arcpy.AddIndex_management(d["nhs_hpms"],f,indnm)
            except:
                pass

        name_no_space = d["name"].replace(" ","_")
        print(name_no_space)
        proj = PROJECTS_FOLDER / "{}.aprx".format(name_no_space)
        if proj.exists():
            proj.unlink()
        #shutil.copy(str(MAP_PROJ),str(proj))
        if arcpy.Exists(str(ws/name_no_space)):
            arcpy.Delete_management(str(ws/name_no_space))
        fcPath = extract_and_project(d["nhs_hpms"],d["wc"],d["epsg"],ws/name_no_space)
        arcpy.RepairGeometry_management(str(fcPath),True)
        if DISSOLVE_FIELDS != None:
            outdiss = ws/"{}_diss".format(name_no_space)
            if arcpy.Exists(str(outdiss)):
                arcpy.Delete_management(str(outdiss))
            fcPath = arcpy.management.Dissolve(str(fcPath),str(outdiss), DISSOLVE_FIELDS, None, "SINGLE_PART", "DISSOLVE_LINES").getOutput(0)
            fcPath = Path(fcPath)
        convert_urban(fcPath)
        add_vmt(fcPath)
        bbox = extent_as_string_latlong_url(fcPath)
        rasterFiles = get_elevation_data(TNM_API_URL,bbox,ELEVATION_DATASET_NAME,DEM_FOLDER)
        if rasterFiles:
            rastname = name_no_space + "_m"
            if arcpy.Exists(str(ws/rastname)):
                arcpy.Delete_management(str(ws/rastname))
            mosaic_and_append(rasterFiles,ws,rastname,fcPath)
        desc = arcpy.Describe(str(fcPath))
        ext = desc.extent
        ext_poly = ext.polygon
        ext_poly = ext_poly.projectAs(arcpy.SpatialReference(4326))
        with arcpy.da.InsertCursor(str(index_map),["SHAPE@","NAME","LABEL"]) as ic:
            ic.insertRow([ext_poly,name_no_space,d["name"]])
        #aadt_layer = arcpy.MakeFeatureLayer_management(str(fcPath),d["name"]).getOutput(0)
        aprx_temp = arcpy.mp.ArcGISProject(str(MAP_PROJ))
        aprx_temp.saveACopy(str(proj))
        del aprx_temp
        aprx_proj = arcpy.mp.ArcGISProject(str(proj))

        #try:
        aprx_map = aprx_proj.listMaps(MAP_NAME)[0]
        aprx_map.name = name_no_space
        aprx_proj.defaultGeodatabase = str(ws)
        aprx_proj.homeFolder = str(PROJECTS_FOLDER)
        #res = aprx_map.addLayer(aadt_layer)
        #if len(res)==1:
            #aadt_layer = res[0]

        lyrx = arcpy.mp.LayerFile(str(CENTERLINE_TEMP))
        print(str(CENTERLINE_TEMP))
        aadt_layer = aprx_map.addLayer(lyrx)[0]
        cp = {'dataset': name_no_space, 'workspace_factory': 'File Geodatabase', 'connection_info': {'database': str(ws)}}
        aadt_layer.name = d["name"]
        aadt_layer.updateConnectionProperties(aadt_layer.connectionProperties,cp)
        aprx_map.spatialReference = arcpy.SpatialReference(d['epsg'])
        aprx_proj.save()
        #aprx_map.defaultView.camera.setExtent(aprx_map.defaultView.getLayerExtent(aadt_layer, False, True))
        #aprx_map.defaultView.zoomToAllLayers()
        data_layers = []
        print("Adding attributes")
        for data in ATTRIBUTE_DATA:
            ndata = plot_tools.plot_data()
            ndata.SOURCE = aadt_layer
            ndata.DATA_TYPE = plot_tools.plot_dt.ATTRIBUTE_DATA
            ndata.LABEL = data["label"]
            if data["color"] != '':
                ndata.add_color(data["color"])
            print(data["field_name"])
            for f in arcpy.ListFields(aadt_layer):
                if f.name == data["field_name"]:
                    if f.type == "String":
                        ndata.add_field(f.name,plot_tools.plot_dt.FIELD_CLASS)
                    elif f.type == "Double" or f.type == "Integer":
                        ndata.add_field(f.name,plot_tools.plot_dt.FIELD_NUMERIC)
                    print("found")
                    break  
            data_layers.append(ndata)

        for data in POLYGON_DATA:
            outfin= ws / "{}_{}".format(name_no_space,data["label"].replace(" ","_"))
            print(outfin)
            if arcpy.Exists(str(outfin)):
                arcpy.Delete_management(str(outfin))
            selection = arcpy.SelectLayerByLocation_management(str(data["source"]),"WITHIN_A_DISTANCE",aadt_layer,STATION_LENGTHS).getOutput(0)
            copied = arcpy.CopyFeatures_management(selection,"in_memory//"+data["label"].replace(" ","_")+"mem").getOutput(0)
            projected = arcpy.Project_management(copied,str(outfin),arcpy.SpatialReference(d["epsg"])).getOutput(0)
            arcpy.Delete_management(copied)
            if data["layer_file"]!=None or data["layer_file"]!="":
                lyrx2 = arcpy.mp.LayerFile(str(data["layer_file"]))
                returned = aprx_map.addLayer(lyrx2)
                srclyr = returned[0]
                cp = {'dataset': "{}_{}".format(name_no_space,data["label"].replace(" ","_")), 'workspace_factory': 'File Geodatabase', 'connection_info': {'database': str(ws)}}
                srclyr.updateConnectionProperties(srclyr.connectionProperties,cp)
                srclyr.name = data["label"].replace(" ","_")
            else:
                pnt_layer = arcpy.MakeFeatureLayer_management(projected,data["label"].replace(" ","_")).getOutput(0)
                srclyr = aprx_map.addLayer(pnt_layer)[0]
            aprx_proj.save()
            print("Selected and projected")
            dt = plot_tools.plot_dt.POLYGON_DATA
            ndata = plot_tools.plot_data()
            ndata.SOURCE = srclyr
            ndata.DATA_TYPE = dt
            ndata.LABEL = data["label"]
            ndata.add_color(data["color"])
            if data["attribute"] != "":
                fld = None
                for f in arcpy.ListFields(srclyr):
                    if f.name == data["attribute"]:
                        if f.type == "String":
                            ndata.add_field(f.name,plot_tools.plot_dt.FIELD_CLASS)
                        elif f.type == "Double" or f.type == "Integer":
                            ndata.add_field(f.name,plot_tools.plot_dt.FIELD_NUMERIC)
                        
                        break
            data_layers.append(ndata)



        print("Adding point layers")
        for data in NEAR_POINT_DATA:
            outfin= ws / "{}_{}".format(name_no_space,data["label"].replace(" ","_"))
            print(outfin)
            if arcpy.Exists(str(outfin)):
                arcpy.Delete_management(str(outfin))
            selection = arcpy.SelectLayerByLocation_management(str(data["source"]),"WITHIN_A_DISTANCE",aadt_layer,STATION_LENGTHS).getOutput(0)
            copied = arcpy.CopyFeatures_management(selection,"in_memory//"+data["label"].replace(" ","_")+"mem").getOutput(0)
            projected = arcpy.Project_management(copied,str(outfin),arcpy.SpatialReference(d["epsg"])).getOutput(0)
            arcpy.Delete_management(copied)
            if data["layer_file"]!=None or data["layer_file"]!="":
                lyrx2 = arcpy.mp.LayerFile(str(data["layer_file"]))
                returned = aprx_map.addLayer(lyrx2)
                srclyr = returned[0]
                cp = {'dataset': "{}_{}".format(name_no_space,data["label"].replace(" ","_")), 'workspace_factory': 'File Geodatabase', 'connection_info': {'database': str(ws)}}
                srclyr.updateConnectionProperties(srclyr.connectionProperties,cp)
                srclyr.name = data["label"].replace(" ","_")
            else:
                pnt_layer = arcpy.MakeFeatureLayer_management(projected,data["label"].replace(" ","_")).getOutput(0)
                srclyr = aprx_map.addLayer(pnt_layer)[0]
            aprx_proj.save()
            print("Selected and projected")
            dt = plot_tools.plot_dt.POINT_DATA
            ndata = plot_tools.plot_data()
            ndata.SOURCE = srclyr
            ndata.DATA_TYPE = dt
            ndata.LABEL = data["label"]
            ndata.SELECTION_DISTANCE = data["search_distance"]
            ndata.add_color(data["color"])
            if data["attribute"] != "":
                fld = None
                for f in arcpy.ListFields(srclyr):
                    if f.name == data["attribute"]:
                        if f.type == "String":
                            ndata.add_field(f.name,plot_tools.plot_dt.FIELD_CLASS)
                        elif f.type == "Double" or f.type == "Integer":
                            ndata.add_field(f.name,plot_tools.plot_dt.FIELD_NUMERIC)
                        
                        break
            data_layers.append(ndata)

        

        pt = plot_tools.plotter()
        pt.INPUT_DATA = data_layers
        pt.LAYOUT_TEMPLATE = LAYOUT_TEMP
        pt.STATIONS_LAYER = STATIONS_LAYER
        pt.CENTERLINE_SOURCE = aadt_layer
        pt.MAP_PROJ = aprx_proj
        pt.BASE_MAP = aprx_map
        pt.STATION_DISTANCES = STATION_SPACING
        pt.TANGENT_LINE_LENGTH = STATION_LENGTHS
        pt.SCRATCH_FOLDER = working_folder
        pt.SCRATCH_GDB = ws
        pt.SPATIAL_REF = arcpy.SpatialReference(d["epsg"])
        pt.FILL_SPACE = FILL_SPACE
        pt.logprint = True
        pt.FIG_W = 8.1
        pt.FIG_H = 6.9
        pt.TEXT_ELEMENTS = {"TitleText":d["name"]}
        pt.CONVERT_TO_MILES = True
        print("Adding Charts")
        pt.add_charts()

        #except:
            #print("Error adding layer")
        aprx_proj.save()
        del aprx_proj

# --------------------------------------------------------------------------------------------------
if __name__ == "__main__":


    main(stripcharts_data.centerlines)

    print('\ndone')
