import arcpy
from scipy.spatial import KDTree
import numpy as np
import datetime
import time
from pathlib import Path
import math

"""
This script calculates the network distance between the parking lot sites and the charging station location. 

Create Date: 2022
Created By: The Volpe Center's Geospatial Team 
"""


################## CLASSES  ##################################################################

class GeographicCoordinateSystemError(Exception):
    def __init__(self,message="Need a projected coordinate system to proceed"):
        self.message = message
        super().__init__(self.message)



class helper_functions(object):

    @staticmethod
    def round_place(n,place=1000):
        return int(math.ceil(n / place)) * place

    @staticmethod
    def drop_add_field(featureClass,field_name,field_type,field_length=None,field_precision=None,field_scale=None):
        try:
            arcpy.AddField_management(featureClass,field_name,field_type,field_precision,field_scale,field_length)
        except:
            arcpy.DeleteField_management(featureClass,field_name)
            arcpy.AddField_management(featureClass,field_name,field_type,field_precision,field_scale,field_length)
        return True

    @staticmethod
    def drop_add_featureclass(workspace:Path, featureClassName:str,featureClassType:str,sr:arcpy.SpatialReference):
        if arcpy.Exists(str(workspace/featureClassName)):
            arcpy.Delete_management(str(workspace/featureClassName))
        result = arcpy.CreateFeatureclass_management(str(workspace),featureClassName,featureClassType,spatial_reference=sr).getOutput(0)
        return result


    @staticmethod
    def drop_add_table(workspace:Path, tableName:str):
        if arcpy.Exists(str(workspace/tableName)):
            arcpy.Delete_management(str(workspace/tableName))
        result = arcpy.CreateTable_management(str(workspace),tableName).getOutput(0)
        return result

    @staticmethod
    def drop_add_fgdb(workspace:Path,fgdbName:str):
        """
            Description:
                Deletes file geodatabase if exists and creates it.
            Parameters:
                workspace (pathlib.Path): folder
                fgdbName (str): name with gdb extension
            returns:
                fgdb path as Path object
        """
        if arcpy.Exists(str(workspace/fgdbName)):
            arcpy.Delete_management(str(workspace/fgdbName))
        result = arcpy.CreateFileGDB_management(str(workspace),fgdbName).getOutput(0)
        return Path(result)


    @staticmethod
    def list_fields_for_cursor(featureClass:str,startFields=["SHAPE@","OID@"]):
        fields = arcpy.ListFields(featureClass)

        if startFields == None:
            startFields = []
        
        startFields += [f.name for f in fields]
        startFields.remove("OBJECTID")
        startFields.remove("Shape")
        startFields.remove("Shape_Length")
        return startFields


class point_index_ids(object):
    WGS84_EPSG = 4326
    
    def __init__(self,dataSource):
        """
        Description:
            Creates a kdtree to search for nearest points along with their oids
        Parameters:
            dataSource (str) path to point feature class TODO: Test for point feature class first
        Methods:
            project_to_wgs84 - projects all points to latitude and longitude
            get_closest_point_fc - finds nearest points. returns in memory feature class
            get_closest_points_ids - finds nearest points to one point based on distance value, returns the points plus their oid
            get_closest_points_ids - finds nearest points to array of coordinates. Returns unique coordinates and oid
            get_closest_kpoint_ids - finds k nearest points to one point.
        """
        self.points = []
        self.ids = []
        with arcpy.da.SearchCursor(dataSource,["SHAPE@XY","OID@"]) as sc:
            for row in sc:
                self.points.append([row[0][0],row[0][1]])
                self.ids.append(row[1])
        self.points = np.array(self.points)
        self.points_lat_long = []
        self.ids = np.array(self.ids)
        self.sr = arcpy.Describe(dataSource).spatialReference
        self.tree = KDTree(self.points)
        self.project_to_wgs84()
        
    def project_to_wgs84(self):
        for p in self.points:
            pg = arcpy.PointGeometry(arcpy.Point(p[0],p[1]),spatial_reference=self.sr)
            pglat = pg.projectAs(arcpy.SpatialReference(self.WGS84_EPSG))
            centroid = pglat.centroid
            self.points_lat_long.append([centroid.X,centroid.Y])
        self.points_lat_long = np.array(self.points_lat_long)
        
    def get_closest_point_fc(self,xcoord,ycoord,name,maxdist=80467.2):
        ii = self.tree.query_ball_point([xcoord,ycoord],r=maxdist)
        in_mem_fc = arcpy.CreateFeatureclass_management("in_memory",name,"POINT",spatial_reference=self.sr)[0]

        with arcpy.da.InsertCursor(in_mem_fc,["SHAPE@"]) as ic:
            
            for v in self.points[ii]:
                pnt = arcpy.PointGeometry(arcpy.Point(v[0],v[1]),spatial_reference=self.sr)
        return in_mem_fc

    def get_closest_points_ids(self,xcoord,ycoord,name,maxdist=80467.2,latlong=True):
        ii = self.tree.query_ball_point([xcoord,ycoord],r=maxdist)
        if latlong:
            return self.points_lat_long[ii],self.ids[ii]
        return self.points[ii],self.ids[ii]
    
    def get_closest_points_ids_all(self,coords,maxdist=161000,latlong=True):
        ii = self.tree.query_ball_point(coords,r=maxdist)
        indices = []
        for i in ii:
            for j in i:
                indices.append(j)
        indices = np.array(list(set(indices)))
        if len(indices)==0:
            print("No charging stations within max distance")
            return [],[]
        if latlong:
            return self.points_lat_long[indices],self.ids[indices]
        return self.points[indices],self.ids[indices]
    
    def get_closest_kpoint_ids(self,xcoord,ycoord,k=1,latlong=True):
        dd,ii= self.tree.query([xcoord,ycoord],k=k)
        if latlong:
            return self.points_lat_long[ii],self.ids[ii]
        return self.points[ii],self.ids[ii]

    def get_closest_kpoint_ids_all(self,coords,k=1,latlong=True):
        dd,ii= self.tree.query(coords,k=k)
        indices = []
        for i in ii:
            for j in i:
                indices.append(j)
        indices = np.array(list(set(indices)))
        if len(indices)==0:
            print("No charging stations within max distance")
            return [],[]
        if latlong:
            return self.points_lat_long[indices],self.ids[indices]
        return self.points[indices],self.ids[indices]

def remove_characters(value):
    value = value.replace(" ","_")
    value = value.replace("&","_")
    value = value.replace("-","_")
    value = value.replace(".","")
    return value


def get_org_agency(sites_pth:Path,agencies:list)->dict:
    sites_pth = str(sites_pth)
    org_ag = {}
    with arcpy.da.SearchCursor(sites_pth,["org","agency","matchid","label","postal"]) as sc:
        for row in sc:
            if row[1] in agencies:
                if row[0]==' ':
                    org_ag[" |{}|{}".format(row[1],row[4])]=[row[2],row[3],row[4]]
                    print(" |{}".format(row[1]))
                else:
                    org_ag["{}|{}|{}".format(row[0],row[1],row[4])]=[row[2],row[3],row[4]]
    return org_ag



INPUTS = Path("C:\\Users\\David.Lamb\\OneDrive - DOT OST\\Documents\\GitHub\\GIS_FLTP_EV_Gap_Analysis\\data\\input")
OUTPUTS = Path("C:\\Users\\David.Lamb\\OneDrive - DOT OST\\Documents\\GitHub\\GIS_FLTP_EV_Gap_Analysis\\data\\output")
OD_FGDB = "od_pairs_dcfc_all.gdb"

IN_GDB = OUTPUTS / "ev_gap_analysis.gdb"
SITES_PTH = IN_GDB / "flma_points_prj_charge_stations"
STATIONS_PTH = IN_GDB / "charging_stations_prj"

DCFC_STATION = "DCFC"
LEV2_STATION = "Level 2"
STATION_TYPE = DCFC_STATION#LEV2_STATION#

APPEND_TO_EXISTING = False #True #

AGENCIES = ['doi_blm', 'doi_bor', 'doi_fws', 'doi_nps', 'usda_fs', 'usace']
#AGENCIES = ['doi_nps']
#AGENCIES = ['doi_bor']
#AGENCIES = ['usda_fs', 'usace']
ROUTE_ID_START = 0 #0
if __name__ == "__main__":
    
    if APPEND_TO_EXISTING == False:
        fgdb = helper_functions.drop_add_fgdb(OUTPUTS,OD_FGDB)
        tbl = helper_functions.drop_add_table(fgdb,"org_agency_route")
        helper_functions.drop_add_field(tbl,"org","TEXT")
        helper_functions.drop_add_field(tbl,"agency","TEXT")
        helper_functions.drop_add_field(tbl,"routefcid","TEXT")
        helper_functions.drop_add_field(tbl,"matchid","TEXT")
        helper_functions.drop_add_field(tbl,"labeltxt","TEXT")
        helper_functions.drop_add_field(tbl,"postal","TEXT")
    else:
        fgdb = OUTPUTS / OD_FGDB
        tbl = fgdb / "org_agency_route"
        tbl = str(tbl)

    org_ag = get_org_agency(SITES_PTH,AGENCIES)
    #print(list(org_ag.keys()))
    wc = f"ev_connector_types_lev = '{STATION_TYPE}'"
    tblFields = ["org","agency","routefcid","matchid","labeltxt","postal"]

    lbl = remove_characters(STATION_TYPE)
    nm = f"charging_stations_{lbl}"
    newStations = fgdb / nm
    if APPEND_TO_EXISTING == False:
        dcfc_pth = arcpy.Select_analysis(str(STATIONS_PTH),str(newStations),wc).getOutput(0)
    else:
        dcfc_pth = str(newStations)


    qf = point_index_ids(dcfc_pth)


    route_id = ROUTE_ID_START
    route_links = []
    for k,v in org_ag.items():
        mid,lbl,postal = v
        print(k)
        org = k.split("|")[0].replace("'","''")
        agency = k.split("|")[1]
        #postal = k.split("|")[2]
        wc = f"org = '{org}' and agency='{agency}' and postal='{postal}'"
        route_name = f"routes_dcfc_{agency}_{route_id}"
        #print(wc)
        
        #Get the origin points from the sites based on org and agency
        origins_xy = []
        with arcpy.da.SearchCursor(str(SITES_PTH),["SHAPE@XY","OID@","org","agency","label","postal"],where_clause=wc) as sc:
            origins = [[row[0],row[1],remove_characters(row[2]),row[3],row[4],row[5]] for row in sc]
            origin_xy = np.array([[row[0][0],row[0][1]] for row in origins])
        #print(len(origins))

        print("Create the Closest Facility object using arcgis.com")
        cf = arcpy.nax.ClosestFacility("https://www.arcgis.com/")
        cnter = 0
        #Add all the origin points as incidents to the CF object
        for o in origins:
            with cf.insertCursor(arcpy.nax.ClosestFacilityInputDataType.Incidents,["SHAPE@","Name"]) as ic:
                pg = arcpy.PointGeometry(arcpy.Point(o[0][0],o[0][1]),spatial_reference=arcpy.SpatialReference(102003))
                pglat = pg.projectAs(arcpy.SpatialReference(qf.WGS84_EPSG))
                ic.insertRow([pg,str(o[1])])


            #pnt,ids = qf.get_closest_points_ids(o[0][0],o[0][1],o[2].replace(" ","_"))
        print("Getting the closest charging stations within 100,000 meters")
        pnts,ids = qf.get_closest_points_ids_all(origin_xy)

        if len(pnts)==0:
            pnts,ids = qf.get_closest_kpoint_ids_all(origin_xy,k=3)
        if len(pnts)>5000:
            pnts = pnts[:5000-len(origins)]
            ids = ids[:5000-len(origins)]
        if len(pnts)>0:
            with cf.insertCursor(arcpy.nax.ClosestFacilityInputDataType.Facilities,["SHAPE@","Name"]) as ic:
                for pnt,idv in zip(pnts,ids):
                    pntg = arcpy.PointGeometry(arcpy.Point(pnt[0],pnt[1]),spatial_reference=arcpy.SpatialReference(4326))
                    ic.insertRow([pntg,str(idv)])

            print("Solve Closest Facility")
            result = cf.solve()
            if not result.solveSucceeded:
                print("Solve Failed")
                print(result.solverMessages(arcpy.nax.MessageSeverity.All))
                route_links.append([org,agency,None,mid,lbl,postal])
            else:
                #try:
                output_path = fgdb / route_name
                result.export(arcpy.nax.ClosestFacilityOutputDataType.Routes,str(output_path))
                with arcpy.da.InsertCursor(tbl,tblFields) as ic:
                    ic.insertRow([org,agency,route_name,mid,lbl,postal])
                #route_links.append([org,agency,route_name,mid])
                #except:
                    #route_links.append([org,agency,None])                                                  

                
        
            del result
            del cf
        else:
            route_links.append([org,agency,None,mid,lbl,postal])

        time.sleep(1)
        route_id+=1
        #if route_id >20:
            #break
    with arcpy.da.InsertCursor(tbl,tblFields) as ic:
        for rl in route_links:
            if len(rl)!=len(tblFields):
                print(rl)
            else:
                ic.insertRow(rl)

                

