import numpy as np
from scipy.spatial import KDTree
import arcpy
from pathlib import Path
import math
import re

from .special_errors import GeographicCoordinateSystemError
from .special_errors import FeatureClassDoesNotExistError
from .special_errors import WorkspaceDoesNotExistError

class helper_functions(object):

    @staticmethod
    def round_place(n,place=1000):
        """
        Rounds number to place.
        Args:
            n (float): number to change
            place (int): what place to round to
        Returns:
            float: rounded numbers
        """
        return int(math.ceil(n / place)) * place

    @staticmethod
    def drop_add_field(featureClass:Path,field_name:str,field_type:str,field_length=None,field_precision=None,field_scale=None)->bool:
        """
        Adds a field. If it exists, deletes the field and then adds it.
        Args:
            featureClass (pathlib.Path): path to feature class
            field properties

        Returns:
            bool: True
        """
        try:
            arcpy.AddField_management(str(featureClass),field_name,field_type,field_precision,field_scale,field_length)
        except:
                arcpy.DeleteField_management(str(featureClass),field_name)
                arcpy.AddField_management(str(featureClass),field_name,field_type,field_precision,field_scale,field_length)

        return True

    @staticmethod
    def drop_add_featureclass(workspace:Path, featureClassName:str,featureClassType:str,sr:arcpy.SpatialReference)->Path:
        """
        Deletes feature class if exists and creates it.
        Args:
            workspace (pathlib.Path): workspace path
            ffeatureClassName (str): name of the feature class
            featureClassType (str): POINT, POLYLINE, POLYGON
            sr (SpatialReference): spatial reference for the feature class

        Returns:
            Path: feature class path
        """
        if arcpy.Exists(str(workspace/featureClassName)):
            arcpy.Delete_management(str(workspace/featureClassName))
        result = arcpy.CreateFeatureclass_management(str(workspace),featureClassName,featureClassType,spatial_reference=sr).getOutput(0)
        return Path(result)

    @staticmethod
    def drop_add_fgdb_table(workspace_fgdb:Path, tableName:str)->Path:
        """
        Deletes feature class if exists and creates it.
        Args:
            workspace (pathlib.Path):file geodatabase workspace
            tableName (str): name of the feature class
        Returns:
            Path: table path
        """
        if arcpy.Exists(str(workspace/tableName)):
            arcpy.Delete_management(str(workspace/tableName))
        result = arcpy.management.CreateTable(str(workspace),tableName).getOutput(0)
        return Path(result)

    @staticmethod
    def drop_add_fgdb(workspace:Path,fgdbName:str)->Path:
        """
        Deletes file geodatabase if exists and creates it.
        Args:
            workspace (pathlib.Path): folder
            fgdbName (str): name with gdb extension
        Returns:
            Path: fgdb path
        """
        if arcpy.Exists(str(workspace/fgdbName)):
            arcpy.Delete_management(str(workspace/fgdbName))
        result = arcpy.CreateFileGDB_management(str(workspace),fgdbName).getOutput(0)
        return Path(result)


    @staticmethod
    def list_fields_for_cursor(featureClass:Path,startFields=["SHAPE@","OID@"])->list:
        """
        Creates a list of fields without OBJECTID, Shape, or Shape_Length. Use startFields to place field names in the beginning.
        Args:
            featureClass (Path): path to feature class
            startFields (list): fields to include
        Returns:
            list: fields in feature class.
        """
        fields = arcpy.ListFields(str(featureClass))

        if startFields == None:
            startFields = []
        
        startFields += [f.name for f in fields]
        startFields.remove("OBJECTID")
        startFields.remove("Shape")
        startFields.remove("Shape_Length")
        return startFields


    @staticmethod
    def check_for_GCS(featureClass:Path)->bool:
        """
        Check if the feature class's coordinate system is projected or geographic
        Args:
            featureClass (Path): path to feature class
        Returns:
            bool: True - if coordinate reference system is geographic; False - if not geographic
        
        """
        if arcpy.Exists(str(featureClass)):
            sr = arcpy.Describe(str(featureClass)).spatialReference
            if sr.type == "Geographic":
                return True
            else:
                return False
        else:
            False
    
    @staticmethod
    def check_for_sameproj(featureClass1:Path,featureClass2:Path)->bool:
        """
        Quick check if to feature classes are in the same projection.
        Args:
            featureClass1 (Path): path to first feature class
            featureClass2 (Path): path to second feature class
        Returns:
            bool: True - if they are the same; False - if they are not the same
        Raises:
            FeatureClassDoesNotExistError: if feature class does not exist.

        """

        if arcpy.Exists(str(featureClass1)):
            if arcpy.Exists(str(featureClass2)):
                sr1 = arcpy.Describe(str(featureClass1)).spatialReference
                sr2 = arcpy.Describe(str(featureClass2)).spatialReference
                if sr1.factoryCode == sr2.factoryCode:
                    return True
                else:
                    return False
            else:
                raise FeatureClassDoesNotExistError()
        else:
            raise FeatureClassDoesNotExistError()



    @staticmethod
    def match_crs(featureClassBase:Path,featureClass:Path,workspace:Path)->str:
        """
        Project feature class to same projection as base.
        Args:
            featureClassBase (Path): this will provide the coordinate reference system.
        Returns:
            str: path to projected feature class
        Raises:
            FeatureClassDoesNotExistError: if feature class does not exist.
            WorkspaceDoesNotExistError: if workspace does not exist.
        """
        fn = workspace / featureClass.stem
        
        if arcpy.Exists(str(workspace)):
            if arcpy.Exists(str(featureClassBase)):
                sr1 = arcpy.Describe(str(featureClassBase)).spatialReference
                #sr2 = arcpy.Describe(str(featureClass2)).spatialReference
                res = arcpy.management.Project(str(featureClass),str(fn),sr1).getOutput(0)
                return res
            else:
                raise FeatureClassDoesNotExistError()
        else:
            raise WorkspaceDoesNotExistError()


    @staticmethod
    def direction(compass_angle) -> str:
        if compass_angle < 22.5:
            return "North"
        if compass_angle >= 22.5 and compass_angle < 67.5:
            return "Northeast"
        if compass_angle >= 67.5 and compass_angle < 112.5:
            return "East"
        if compass_angle >= 112.5 and compass_angle < 157.5:
            return "Southeast"
        if compass_angle >= 157.5 and compass_angle < 202.5:
            return "South"
        if compass_angle >= 202.5 and compass_angle < 247.5:
            return "Southwest"
        if compass_angle >= 247.5 and compass_angle < 292.5:
            return "West"
        if compass_angle >= 292.5 and compass_angle < 337.5:
            return "Northwest"
        if compass_angle >= 337.5:
            return "North"

    @staticmethod
    def linear_angle_cardinal(featureClass:str)->tuple:
        """
        Description custom version of Esri's linear directional mean:
        https://pro.arcgis.com/en/pro-app/2.8/tool-reference/spatial-statistics/h-how-linear-directional-mean-spatial-statistics-w.htm
        
        Args:
            featureClass (str): path to polyline feature class
        
        Returns:
            tuple: List of compass angles, list of cardinal directions
        Raises:
            GeographicCoordinateSystemError: if feature class uses GCS
        """
        if arcpy.Exists(featureClass):
            isgcs = helper_functions.check_for_GCS(featureClass)
            if isgcs == True:
                raise GeographicCoordinateSystemError
            with arcpy.da.SearchCursor(featureClass,["SHAPE@"]) as sc:
                lines = [row[0] for row in sc]
            compass_angles = []
            compass_direc = []
            ldm = []
            for cl in lines:
                seg_order = []
                for part in cl:
                    for pnt in part:
                        seg_order.append([pnt.X,pnt.Y])
                seg_order = np.array(seg_order)
                seg_delta = np.diff(seg_order,axis=0)
                seg_theta = np.arctan2(seg_delta[:,0],seg_delta[:,1])
                #print(seg_theta)
                seg_theta_deg = np.mean(np.rad2deg(seg_theta))
                #print(seg_theta_deg)
                if seg_theta_deg>=0:
                    compass_angles.append(seg_theta_deg)
                elif seg_theta_deg <0:
                    compass_angles.append(360 + seg_theta_deg)
            #     seg_theta = np.arctan2(seg_delta[:,0],seg_delta[:,1])
            #     seg_theta = np.arctan(seg_delta[:,0]/seg_delta[:,1])
            #     #compass_angles.append(np.degrees(seg_theta)[0])
            #     seg_cos = np.cos(np.arctan2(seg_delta[:,0],seg_delta[:,1]))
            #     seg_sin = np.sin(np.arctan2(seg_delta[:,0],seg_delta[:,1]))
            #     seg_ldm = np.arctan(np.sum(seg_sin)/np.sum(seg_cos))
            #     #seg_ldm = np.rad2deg(seg_ldm)
            #     ldm.append(seg_ldm)
                
            #     if seg_sin>= 0 and seg_cos > 0:
            #         compass_angles.append(np.rad2deg(seg_ldm))
            #     elif np.sum(seg_sin)>= 0 and np.sum(seg_cos) < 0:
            #         compass_angles.append(180 - np.rad2deg(seg_ldm))
            #     elif np.sum(seg_sin)< 0 and np.sum(seg_cos) > 0:
            #         compass_angles.append(360 - np.rad2deg(seg_ldm))
            #     elif np.sum(seg_sin)< 0 and np.sum(seg_cos) < 0:
            #         compass_angles.append(180 + np.rad2deg(seg_ldm))
            #     else:
            #        compass_angles.append(None)
            compass_direc = [helper_functions.direction(ca) for ca in compass_angles]
            return compass_angles,compass_direc

    @staticmethod
    def overall_orientation_cardinal(featureClass:str)->tuple:
        """
        Description custom version of Esri's linear directional mean:
        https://pro.arcgis.com/en/pro-app/2.8/tool-reference/spatial-statistics/h-how-linear-directional-mean-spatial-statistics-w.htm
        
        Args:
            featureClass (str): path to polyline feature class
        Returns:
            tuple: List of compass angles, list of cardinal directions
        Raises:
            GeographicCoordinateSystemError: if feature class uses GCS
        """
        if arcpy.Exists(featureClass):
            isgcs = helper_functions.check_for_GCS(featureClass)
            if isgcs == True:
                raise GeographicCoordinateSystemError
            with arcpy.da.SearchCursor(featureClass,["SHAPE@"]) as sc:
                lines = [row[0] for row in sc]
            compass_angles = []
            compass_direc = []
            for cl in lines:
                seg_order = []
                for pnt in [cl.firstPoint,cl.lastPoint]:
                    seg_order.append([pnt.X,pnt.Y])
                seg_order = np.array(seg_order)
                seg_delta = np.diff(seg_order,axis=0)
                seg_theta = np.arctan2(seg_delta[:,0],seg_delta[:,1])
                #print(seg_theta)
                seg_theta_deg = np.mean(np.rad2deg(seg_theta))
                #print(seg_theta_deg)
                if seg_theta_deg>=0:
                    compass_angles.append(seg_theta_deg)
                elif seg_theta_deg <0:
                    compass_angles.append(360 + seg_theta_deg)
            compass_direc = [helper_functions.direction(ca) for ca in compass_angles]
            return compass_angles,compass_direc


    @staticmethod
    def point_angles_degrees(center_point:np.array,points:np.array)->list:
        """
        Description custom version of Esri's linear directional mean:
        https://pro.arcgis.com/en/pro-app/2.8/tool-reference/spatial-statistics/h-how-linear-directional-mean-spatial-statistics-w.htm
        
        Args:
            center_point (tuple): x,y
            points (np.array): array of coordinates
        
        Returns:
            list: compass angles

        """
        delta = center_point - points
        sigma = np.rad2deg(np.arctan2(delta[:,0],delta[:,1]))
        compass_angles = list(np.where(sigma<0,sigma+360,sigma))
        return compass_angles


    @staticmethod
    def inbound_outbound(inboundPoint:arcpy.Point,featureClass:str):
        """
        Todo:
            Not implemented. Needs testing.
        Raises:
            GeographicCoordinateSystemError: if feature class uses GCS
        """
        if arcpy.Exists(featureClass):

            isgcs = helper_functions.check_for_GCS(featureClass)
            if isgcs == True:
                raise GeographicCoordinateSystemError

            with arcpy.da.SearchCursor(featureClass,["SHAPE@"]) as sc:
                lines = [row[0] for row in sc]
            in_out = []
            for cl in lines:
                start_dist = inboundPoint.distanceTo(cl.firstPoint)
                end_dist = inboundPoint.distanceTo(cl.lastPoint)
                delta_y = cl.firstPoint.Y - cl.lastPoint.Y
                delta = start_dist / end_dist
                #if delta> -100 and delta < 100:
                    #in_out.append("UNK")
                if start_dist < end_dist:
                    in_out.append("OUTBOUND")
                elif start_dist > end_dist:
                    in_out.append("INBOUND")
            if len(lines)==len(in_out):
                return in_out
            else:
                print("problem")
    @staticmethod
    def convert_rows_columns_to_coordinates(xmin:float,ymax:float,cw:float,ch:float,row:int,column:int,dx=0,dy=0)->tuple:
        """
        Converts the column and row to an x and y coordinate
        Args:
            xmin (float): the minimum x value of the raster (upper left)
            ymax (float): the maximum y value of the raster (upper left)
            cw (float): cell width 
            ch (float): cell height
            row (int): index of the row
            column (int): index of the column
            dx (float): offset if needed in the x direction
            dy (float): offset if needed in the y direction
        Returns:
            tuple: the x and y coordinates"""

        return (xmin + (.5 * cw) + (column * cw)+dx,ymax - (.5 * ch) - (row * ch)+dy)

    @staticmethod
    def process_raster_layer_values_and_coords(raster_pth:Path)->tuple:
        """
            Process one raster layer to retreive values and coords with data. Skips No data.
        Args:
            raster_pth (Path): path to the raster layer
        Returns:
            tuple: Value of each cell (numpy array), Coords of each cell (numpy array)
        Raises:
            GeographicCoordinateSystemError: if feature class uses GCS
        """
        if helper_functions.check_for_GCS(raster_pth) == True:
            raise GeographicCoordinateSystemError()

        r_c = arcpy.Raster(str(raster_pth))

        sr_c = r_c.spatialReference.factoryCode
        
        xmin = r_c.extent.XMin
        ymax = r_c.extent.YMax
        cw = r_c.meanCellWidth
        ch = r_c.meanCellHeight

        #load into numpy arrays for analysis

        #create raster cell iterator, much faster than raster to numpy function, skip the nodata cells.
        with arcpy.sa.RasterCellIterator({'rasters':[r_c],'skipNoData':[r_c]}) as rci_skip:
            values = np.array([ rci_skip[i, j] for i, j in rci_skip])
            coords = np.array([helper_functions.convert_rows_columns_to_coordinates(xmin,ymax,cw,ch,i,j) for i, j in rci_skip])

        return values, coords


    @staticmethod
    def dms_to_dd(str_dms:str)->tuple:
        """
        Converts a string of latitude and longitude in degrees minutes seconds to decimal degrees.
        Args:
            str_dms (str): Uses the string format 42° 11' 3.171" N / 73° 24' 11.666" W
        Returns:
            tuple: latitude and longitide
        """
        
        lat_base = str_dms.split(" / ")[0]
        long_base = str_dms.split(" / ")[1]
        output = []
        for x in [lat_base,long_base]:
            re_pattern = r"(\d+)\s?\°\s?(\d+)\s?\'\s?(\d{1,}\.?\,?\d{0,}?)\"\s?(N|W|S|E)"
            result = re.search(re_pattern,x)
            minutes = float(result.group(2)) / 60
            seconds = float(result.group(3)) / 3600
            degrees =float(result.group(1))+minutes+seconds
            if result.group(4) in ("W","S"):
                degrees*=-1
            output.append(degrees)
        return output


    @staticmethod
    def geodesic(lat1:float,lon1:float,lat2:float,lon2:float)->float:
        """
        Geodesic distance
        Args:
            lat1 (float)
            lon1 (float)
            lat2 (float)
            lon2 (float)

        returns
            float: distance in meters.
        """
        import math
        R = 6371000
        phi1 = lat1 * math.pi / 180
        phi2 = lat2 * math.pi / 180
        deltaphi = (lat2 - lat1) * math.pi / 180
        deltalambda = (lon2 - lon1) * math.pi / 180
        a = math.sin(deltaphi / 2) * math.sin(deltaphi/2) + math.cos(phi1) * math.cos(phi2) * math.sin(deltalambda/2) * math.sin(deltalambda/2)
        c = 2 * math.atan2(math.sqrt(a),math.sqrt(1-a))
        d = R * c
        return d
