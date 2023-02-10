
import numpy as np
from scipy.spatial import KDTree
import arcpy

class polyline_kdtree_index(object):
    def __init__(self,base_data:list,shape_index=0):
        """
        Small class that creates a spatial index of vertices on a polyline, using a list of polylines. Works with arcpy.
        Call `build_tree()` after initializing a class.

        Args:
            base_data (list): list of lists.
            shape_index (int): index position of the polyline
        Properties:
            tree (KDTree): spatial index of points of polylines.
            point_to_polyline_index (np.array): index values corresponding points to polyline base data index.
        Methods:
            build_tree: call to build the KDTree

        Examples:
            all_lines is a list of lists. It should at the least contain a polyline object.
            all_lines_index refers to the position of the input data all_lines.
            
            :param base_data: list
            :param shape_index: int

            >>> all_lines = [[arcpy.Polyline,"Road"],[arcpy.Polyline,"Road"]]
            >>> pki= polyline_kdtree_index(all_lines,shape_index=0)
            >>> pki.build_tree()
            >>> a_point = [1,0]
            >>> all_lines_index = pki.nearest_polyline_single_point((xcoord,ycoord))
            
            
        """
        self.base_data = base_data
        self.tree = None
        self.point_to_polyline_index = np.empty(1)
        self.all_coords = np.empty(1)
        self.shape_index = shape_index
        
    def grab_polyline_points(self,index:int,polyline:arcpy.Polyline):
        """
        Gets polyline's vertices
        Args:
            index (int): index of the polyline to associate the coordinates with
            polyline (arcpy.Polyline): polyline object, arcpy geometry.
        Returns:
            tuple: list of index, list of points as tuples

        """
        index_result = []
        points = []
        for part in polyline:
            # Step through each vertex in the feature
            for pnt in part:
                if pnt:
                    points.append([pnt.X, pnt.Y])
                    index_result.append(index)
        return index_result,points
    
    def build_tree(self)->None:
        """
        Builds KDTree for the points. Should be run before querying the index.
        Args:
            None
        Returns:
            None
        """
        indices = []
        coords = []
        for i,r in enumerate(self.base_data):
            indx,shps = self.grab_polyline_points(i,r[self.shape_index])
            indices+=indx
            coords+= shps
        #self._raw_coords = coords
        self.point_to_polyline_index = np.array(indices)
        self.all_coords = np.array(coords)
        self.tree = KDTree(self.all_coords)
    
    def nearest_polyline_single_point(self,pnt=None)->int:
        """
        Gets the index of the polyline that is closest to the input point.
        
        Args:
            pnt (tuple or arcpy.Point): the query point
        Returns:
            int: index of the polyline in the input data
        """
        try:
            pnt = [pnt.X,pnt.Y]
        except:
            try:
                pnt = [pnt.centroid.X,pnt.centroid.Y]
            except:
                pass
        
        distance, idx_nearest_point = self.tree.query(pnt,k=1)
        idx_nearest_line = self.point_to_polyline_index[idx_nearest_point]
        return idx_nearest_line
        
    def nearest_polylines_from_points(self,pnts:np.array,cutoff_dist:float)->list:
        """
        Retrieves index of nearest points using a cutoff distance.
        Parameters:
            pnts (np.array): array of xy tuples.
        returns:
            np.array: Array of indices.
        """
        dists, idxs= self.tree.query(pnts,k=1)
        idx_nearest_lines = [self.point_to_polyline_index[i] if d<=cutoff_dist else None for i,d in zip(idxs,dists)]
        return idx_nearest_lines
    
    def rows_from_points(self,pnts:np.array,cutoff_dist:float):
        idx_nearest_lines = self.nearest_polylines_from_points(pnts,cutoff_dist)
        return [self.base_data[row] if row!=None else None for row in idx_nearest_lines]

    def get_nearest_point_on_line(self,pnt=None):
        """
        Gets the nearest point on the nearest line
        Args:
            pnt (np.array): tuple of xy coordinates
        Returns:
            tuple: (point,distance from start of line, distance to line from in point, right side)
        Todo:
            Needs to be tested, not sure if this is implemented correctly.
        """
        if type(pnt)=="arcpy.arcobjects.arcobjects.Point":
            qpnt = [pnt.X,pnt.Y]
        elif type(pnt)=="arcpy.arcobjects.arcobjects.PointGeometry":
            qpnt = [pnt.centroid.X,pnt.centroid.Y]
        else:
            qpnt = pnt
            pnt = arcpy.Point(pnt[0],pnt[1])
        nearest_line_index = self.nearest_polyline_single_point(qpnt)
        polyline = self.base_data[nearest_line_index]
        nearest_point_distance = polyline.queryPointAndDistance(pnt)
        return nearest_point_distance
    
    def append_lines(self,rows_to_add:list):
        """
        Recreate the list of lines by appending more lines. Builds the tree.
        Args:
            rows_to_add (list): should match the input data.
        Returns:
            bool: True once complete.
        Todo:
            Needs to be tested, not sure if this is implemented correctly.
        """
        coords = []
        indices = []

        for i,r in enumerate(rows_to_add):
            self.base_data.append(r)
            indx,shps = self.grab_polyline_points(i+len(self.base_data)+1,r[self.shape_index])
            indices+=indx
            coords+= shps


        np.append(self.point_to_polyline_index,indices)
        np.append(self.all_coords,coords)
        del self.tree
        self.tree = KDTree(self.all_coords)
        return True

    def find_indices_by_field(self,field_index,matching_values):
        """
        Todo:
            For testing. Not implemented.
        """
        indices = []
        for i,row in self.base_data:
            if row[field_index] in matching_values:
                indices.append(i)
        if len(indices) == 0:
            return None
        else:
            return indices


    def remove_lines_rebuild_tree(self,indices):
        """
        Todo:
            For testing. Not implemented.
        """
        for i in indices:
            del self.base_data[i]
        
        np.delete(self.all_coords,indices)
        del self.tree
        self.tree = KDTree(self.all_coords)

        return True