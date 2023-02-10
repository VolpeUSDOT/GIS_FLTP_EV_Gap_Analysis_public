
from pathlib import Path
from enum import Enum
import numpy as np
import arcpy
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import FormatStrFormatter, StrMethodFormatter
import matplotlib.patches as mpatches
import arcgis_tools.indextools as index_tools
import arcgis_tools.static_tools as static_tools

import matplotlib.transforms as transforms

"""
This tool creates strip chart plots using ArcGIS Pro and Matplotlib. 

Create Date: 2022
Created By: The Volpe Center's Geospatial Team 
"""

class plot_dt(Enum):
    """
    Data formats for use in the plot_data class.
    LINE_DATA - Polyline feature class
    POINT_DATA - Point feature class
    POLYGON_DATA - Polygon feature class
    ATTRIBUTE_DATA - Field in the centerline feature class
    FIELD_CLASS - type of field is a string
    FIELD_NUMERIC - type of field is integer or float
    """
    STRIP_LINE_DATA = 1
    LINE_DATA = 2
    POINT_DATA = 3
    POLYGON_DATA = 4
    ATTRIBUTE_DATA = 5
    RASTER_DATA = 6
    FIELD_CLASS = 7
    FIELD_NUMERIC = 8



class plot_data(object):
    """
    Data object read by the plotter. Sets up the type of data, color, and other information.
    Attributes:
        FIELD_NAMES (list): name of the fields if used.
        FIELD_TYPES (list): type of fields (should be plot_dt.FIELD_CLASS or plot_dt.FIELD_NUMERIC)
        DATA_TYPE (list): should be LINE_DATA, POLYGON_DATA, POINT_DATA or ATTRIBUTE_DATA from plot_dt
        SOURCE (str or Path): the data source for the data
        SR (arcpy.SpatialReference): spatial refrence for the source layer
        LABEL (str): The axis label for the chart being produced
        COLORS (str or dict): colors for the chart. Can be a string (colors separated by commas) or dictionary. Dictionary should be for class values.
        MARKERS (str or dict): TO DO

    """
    
    def __init__(self):
        self.FIELD_NAMES = []
        self.FIELD_TYPES = []
        self.DATA_TYPE = None
        self.SOURCE = None
        self.SR = None
        self.LABEL = "LABEL FOR CHART"
        self.COLORS = {"BaseColor":"Black"}
        self.MARKERS = {"BaseMarker":"."}
        self._values = {}
        self._classes = {}
        self._colors = {}
        self._markers = {}
        self._markers_base = [".","o","v","^","<",">","s","p","P","*"]
        self._create_legend = False
        self.SELECTION_DISTANCE = None

    def add_field(self,name:str,fieldtype):
        """
        Add a field to the plot_data type.
        Args:
            name (str): name of the field
            fieldtype (plot_dt): data type from plot_dt. Either FIELD_CLASS or FIELD_NUMERIC
        Returns:
            None
        """
        if fieldtype in [plot_dt.FIELD_CLASS,plot_dt.FIELD_NUMERIC]:
            self.FIELD_NAMES.append(name)
            self.FIELD_TYPES.append(fieldtype)
            if fieldtype == plot_dt.FIELD_CLASS:
                if self.DATA_TYPE == plot_dt.POINT_DATA:
                    self.create_classes_color(name)


    def add_color(self,colors):
        """
        Add colors to the plot_data.
        Args:
            colors (str or dict): color, or colors separated by comma. Dict is for matching color to a field value. Should be hex colors.
        Returns:
            None
        """
        if type(colors) == str or type(colors)==list:
            if type(colors) == str:
                colors = colors.split(",")
            for i,c in enumerate(colors):
                if c[0]!="#":
                    self.COLORS[i] = "#"+c
                else:
                    self.COLORS[i] = c
        elif type(colors) == dict:
            for k,v in colors.items():
                if v[0]!="#":
                    self.COLORS[k] =  "#"+v
                else:
                    self.COLORS[k] = v
                
        else:
            pass


    def create_classes_color(self,fname):
        """
        Match colors to different classes in the field.
        Args:
            fname (str): name of the field.
        Returns:
            None
        """
        arcpy.AddMessage(self.COLORS)
        arcpy.AddMessage(self._colors)
        if self.SOURCE != None:
            arcpy.AddMessage("Source is not equal to None")
            with arcpy.da.SearchCursor(str(self.SOURCE),[fname]) as sc:
                self._values[fname] =[]
                for row in sc:
                    if row[0]:
                        self._values[fname].append(row[0])

            self._classes[fname] = list(set(self._values[fname]))
            if self._classes[fname]:
                self._classes[fname].sort()
                #self._values[fname] = [self._classes[fname].index(x) for x in self._values[fname]]
                self._colors[fname] = []
                #print(self._classes)
                #print(self._values)
                for x in self._values[fname]:
                    added = False
                    for i,c in enumerate(self._classes[fname]):
                        #print("matching value to colors: {} to {}".format(x,c))
                        if x == c:
                            if c in self.COLORS:
                                colormatch = self.COLORS[c]
                            elif i in self.COLORS:
                                colormatch = self.COLORS[i]
                            else:
                                colormatch = self.COLORS["BaseColor"]
                            self._colors[fname].append(colormatch)
                            added = True
                            break
                    if added == False:
                        colormatch = self.COLORS["BaseColor"]
                        self._colors[fname].append(colormatch)
                #arcpy.AddMessage(self._colors)
                self._patches = []
                for i,c in enumerate(self._classes[fname]):
                    if c in self.COLORS:
                        colormatch = self.COLORS[c]
                    elif i in self.COLORS:
                        colormatch = self.COLORS[i]
                    else:
                        colormatch = self.COLORS["BaseColor"]
                    self._patches.append(mpatches.Patch(color=colormatch, label=c))
                
                self._create_legend = True
            #arcpy.AddMessage("processed classes colors")
            #arcpy.AddMessage(self._colors)
            # else:
            #     self._markers[fname] = []
            #     for x in self._values[fname]:
            #         for i,c in enumerate(self._classes[fname]):
            #             if x == c:
            #                 self._markers[fname].append(self._markers_base[i])
            #     self._patches = []
            #     for i,c in enumerate(self._classes[fname]):
            #         self._patches.append(mpatches.Patch(marker=self._markers_base, label=c))
            #     self._create_legend = True
            


class MultipartGeometryError(Exception):
    def __init__(self,message="Multipart geometry created for the main path. Should be single line."):
        self.message = message
        super().__init__(self.message)



class plotter(object):

    """
    Creates the strip chart plot for an ArcGIS Project.

    Attributes:
        MAP_PROJ (arcpy.mp.ArcGISProject): map project object.
        BASE_MAP (arcpy.mp.Map): basemap map object. New layers will be added here.
        CENTERLINE_SOURCE (str or Path): the line for the strip chart. This will be oriented from left to right, and stations will be added.
        LAOUT_TEMPLATE (str or Path): pagx file. Needs a map frame and an image. The image dimensions should match FIG_W and FIG_H properties.
        STATIONS_LAYER (str or Path): lyrx file. For the stations data.
        SCRATCH_FOLDER (Path): where the images will be saved.
        SCRATCH_GDB (Path): where the new layers will be created.
        SPATIAL_REF (arcpy.SpatialRefrence): base spatial reference from the centerline source.
        TANGENT_LINE_LENGTH (int): length of the station on one side of the line.
        STATION_DISTANCES (int): spacing of the stations.
        FIG_W (float): dimensions of the plot width.
        FIG_H (float): dimensions of the plot height.
        FILL_SPACE (bool): True - fill the space with the plot, False - divide space by 5.
        TEXT_ELEMENTS (dict): name of the text element on the layout, value is the text value.
        FONT_SIZE (int): font size for all parts of the plots.
    """

    INPUT_DATA = []
    
    MAP_PROJ = None
    BASE_MAP = None
    CENTERLINE_SOURCE = None
    LAYOUT_TEMPLATE = None
    STATIONS_LAYER = None
    SCRATCH_FOLDER = None
    SCRATCH_GDB = None
    SPATIAL_REF = None

    TANGENT_LINE_LENGTH = 100
    STATION_DISTANCES = 500
    
    FIG_W = 8.1
    FIG_H = 7.4

    CONVERT_TO_MILES = False

    #Highlight Ranges
    HIGHLIGHT_RANGES = {"spans":[],"colors":[],"labels":[]}

    #Plot Properties
    Y_LABEL_RIGHT = True
    

    def __init__(self):
        """
        Args:
            max_charts_per_page (int): default is 5 charts per page.
            logprint (bool): True print out the messages, false use arcpy.AddMessage()
        """
        self.__version__ = .1
        self._single_path =None
        self._left_point = None
        self._right_point = None
        self._current_layout = None
        self._current_map_frame = None
        
        self.max_charts_per_page = 5
        self.current_charts_on_page = 5
        self.chart_counter = 0
        self.total_length = 0
        self.total_length_inches = 0
        self.FILL_SPACE = False
        self.FONT_SIZE = 7
        self.logprint = False
        self.TEXT_ELEMENTS = None
        self.LEGENDS = False

    def message(self, statement:str,error=False,warning=False):
        """
        Prints message or sends to arcpy.AddMessage()
        
        """
        if self.logprint == True:
            print(statement)
        else:
            if error:
                arcpy.AddError(statement)
            elif warning:
                arcpy.AddWarning(statement)
            else:
                arcpy.AddMessage(statement)

    def add_charts(self):
        """
        Run add_charts to create the plots and add to the layout. Multiple pages created if needed.
        """
        mapbasename = self.BASE_MAP.name
        mapbasename = mapbasename.replace(" ","_")


        self._lines = [[row[0],row[1]] for row in arcpy.da.SearchCursor(self.CENTERLINE_SOURCE,["SHAPE@","OID@"])]
        pnts = self.lines_to_points()
        _,_,ulline = self.bounding_box(pnts)
        sorted_lines = self.sort_lines(ulline)
        self.message("Sorted Lines")
        start_end_points = np.array([[self._left_point.X,self._left_point.Y],[self._right_point.X,self._right_point.Y]])
        rot = self.rotation(start_end_points)
        
        self.stations_fc = static_tools.helper_functions.drop_add_featureclass(self.SCRATCH_FOLDER,"{}_stations.shp".format(mapbasename),"POLYLINE",self.SPATIAL_REF)
        self.message("Created Stations")
        static_tools.helper_functions.drop_add_field(self.stations_fc,"Station","DOUBLE")
        self.create_stations(sorted_lines)


        self.message("Adding stations layer")
        lyrx = arcpy.mp.LayerFile(str(self.STATIONS_LAYER))
        addedlyr = self.BASE_MAP.addLayer(lyrx)[0]
        
        cp = {'dataset': "{}_stations.shp".format(mapbasename), 'workspace_factory': 'Shape File', 'connection_info': {'database': str(self.SCRATCH_FOLDER)}}

        addedlyr.updateConnectionProperties(addedlyr.connectionProperties,cp)

        #how many charts do we need:
        self.MAP_PROJ.save()
        pages = int(len(self.INPUT_DATA) / self.max_charts_per_page)+1
        if len(self.INPUT_DATA) == 5:
            pages = 1
        #Start by figuring out
        chart_running_count = 0
        for page in range(pages):
            self.message("Page {}".format(page))
            leg_patches = []
            if self.import_stripchart_layout():
                lytnm = "{}_stripchart_{}".format(mapbasename,page)
                self._current_layout.name = lytnm
                self.update_map(rot)
                running_ax = (len(self.INPUT_DATA)-chart_running_count)
                if running_ax<self.max_charts_per_page:
                    if self.FILL_SPACE == False:
                        tempheight = self.FIG_H / 5
                        newheight = np.sum([tempheight for x in range(running_ax)])
                    else:
                        newheight = self.FIG_H
                    gridspec = {"width_ratios":[1],'height_ratios':[1 for x in range(running_ax)]}
                    fig, axes = plt.subplots(ncols=1, nrows=running_ax, figsize=(self.FIG_W,newheight),gridspec_kw=gridspec,dpi=300,sharex='all')
                    if running_ax == 1:
                        axes=[axes]
                else:
                    gridspec = {"width_ratios":[1],'height_ratios':[1,1,1,1,1]}
                    fig, axes = plt.subplots(ncols=1, nrows=self.max_charts_per_page, figsize=(self.FIG_W,self.FIG_H),gridspec_kw=gridspec,dpi=300,sharex='all')
                #fig.subplots_adjust(0.05,0.05,0.95,0.95, wspace=0.00, hspace=.05)
                
                for ax in axes:
                    data_info = self.INPUT_DATA[chart_running_count]
                    self.message(data_info.DATA_TYPE)
                    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
                    if self.Y_LABEL_RIGHT:
                        ax.yaxis.set_label_position("right")
                    if data_info.DATA_TYPE == plot_dt.ATTRIBUTE_DATA:
                        if len(data_info.FIELD_NAMES)==1:
                            fldname = data_info.FIELD_NAMES[0]
                            chrttype = data_info.FIELD_TYPES[0]
                            values = [row[0] for row in arcpy.da.SearchCursor(self.CENTERLINE_SOURCE,[fldname])]
                            colormatch = "Black"
                            if 0 in data_info.COLORS:
                                colormatch = data_info.COLORS[0]
                            if chrttype == plot_dt.FIELD_CLASS:
                                yord = list(set(values))
                                yord.sort()
                                values_n = [yord.index(x) for x in values]
                                self.create_att_line_plots(values_n,sorted_lines,ax,colormatch,add_yticklabels=True)
                                
                                if len(yord) >= 1:
                                    ax.set_ylim([-.1,len(yord)-1+.1])
                                    ax.set_yticks([i for i in range(0,len(yord))])
                                    ax.set_yticklabels(yord)

                            if chrttype == plot_dt.FIELD_NUMERIC:
                                self.create_att_line_plots(values,sorted_lines,ax,colormatch)
                                ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
                            ax.set_ylabel(data_info.LABEL)
                    if data_info.DATA_TYPE == plot_dt.POINT_DATA:
                        if data_info.SELECTION_DISTANCE == None:
                                data_info.SELECTION_DISTANCE = self.TANGENT_LINE_LENGTH * 2
                        if len(data_info.FIELD_NAMES) == 0:
                            colormatch = ["Black"]
                            if 0 in data_info.COLORS:
                                colormatch = data_info.COLORS[0]
                            
                            self.create_plot_from_points(data_info.SOURCE,sorted_lines,ax,data_info.SELECTION_DISTANCE,colormatch)
                            ax.set_ylim((0,data_info.SELECTION_DISTANCE+1))
                            ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
                            ax.set_ylabel(data_info.LABEL)
                        elif len(data_info.FIELD_NAMES) == 1:
                            if data_info.FIELD_TYPES[0] == plot_dt.FIELD_CLASS:
                                self.message("Class Point")
                                fname = data_info.FIELD_NAMES[0]
                                try:
                                    colors = data_info._colors[fname]
                                except:
                                    colors = ["Black"]

                                markers = 'o'
                                self.message(data_info._patches)
                                self.create_plot_from_points(data_info.SOURCE,sorted_lines,ax,data_info.SELECTION_DISTANCE,color=colors,markers=markers)
                                leg_patches +=data_info._patches
                                
                                ax2 = ax.twinx()
                                ax2.set_yticks([])
                                ax2.set_ylabel(data_info.LABEL)
                                ax2.yaxis.set_label_position("right")
                                ax.yaxis.set_label_position("left")
                                #ax.set_ylim((0,data_info.SELECTION_DISTANCE+1))
                                if self.CONVERT_TO_MILES:
                                    ax.set_ylim((-.5,self.convert_length_to_miles(data_info.SELECTION_DISTANCE)+1))
                                    yticks = [0,data_info.SELECTION_DISTANCE*.5,data_info.SELECTION_DISTANCE]
                                    yticks = [self.convert_length_to_miles(l) for l in yticks]
                                    ax.set_yticks(yticks)
                                    #nlbls = []
                                    #for l in yticks:
                                    #    print(l)
                                    #    nlbls.append(str(int(self.convert_length_to_miles(l))))
                                    #ax.set_yticklabels(nlbls)
                                    ax.set_ylabel("Distance to Centerline\n(Miles)")
                                    #print(nlbls)
                                else:
                                    ax.set_ylabel("Distance to Centerline\n({})".format(self.unit_label))
                                
                                ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
                                ax.margins(.2) 
                                ax2.margins(.2) 
                                self.set_font_size(ax2,self.FONT_SIZE)
                            else:
                                self.message("Handling point data with a numeric attribut is not implemented at the moment...",warning=True)
                                colormatch = ["Black"]
                                if 0 in data_info.COLORS:
                                   colormatch =  data_info.COLORS[0]
                                self.create_plot_from_points(data_info.SOURCE,sorted_lines,ax,data_info.SELECTION_DISTANCE,colors)
                                ax.set_ylim((0,data_info.SELECTION_DISTANCE+1))
                                ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
                                ax.set_ylabel(data_info.LABEL)
                    if data_info.DATA_TYPE == plot_dt.POLYGON_DATA:
                        self.message(data_info.FIELD_NAMES)
                        if len(data_info.FIELD_NAMES) == 1:
                            if data_info.FIELD_TYPES[0] == plot_dt.FIELD_CLASS:
                            
                                colormatch = ["Black"]
                                fname = data_info.FIELD_NAMES[0]
                                if 0 in data_info.COLORS:
                                    colormatch = data_info.COLORS[0]
                                self.create_plot_from_polygons(data_info.SOURCE,sorted_lines,ax,fname,color=colormatch)
                                
                                ax.set_ylabel(data_info.LABEL)
                            else:
                                self.message("A polygon feature class with numeric variable is not supported at the moment...")

                    self.message("Chart running count: {}".format(chart_running_count))
                    chart_running_count+=1
                    self.set_font_size(ax,self.FONT_SIZE)
                    xlim = ax.get_xlim()
                    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
                    for span,spanColor,spanLabel in zip(self.HIGHLIGHT_RANGES["spans"],self.HIGHLIGHT_RANGES["colors"],self.HIGHLIGHT_RANGES["labels"]):
                        ax.axvspan(span[0], span[1], color=spanColor, alpha=0.4)

                self.message("Finished adding plots for page {}".format(page))
                if self.CONVERT_TO_MILES:
                    ax.set_xlabel("Distance (Miles)")
                else:
                    ax.set_xlabel("Distance ({})".format(self.unit_label))
                for span,spanColor,spanLabel in zip(self.HIGHLIGHT_RANGES["spans"],self.HIGHLIGHT_RANGES["colors"],self.HIGHLIGHT_RANGES["labels"]):
                        leg_patches.append(mpatches.Patch(color=spanColor, label=spanLabel,alpha=.4))

                if len(leg_patches) > 0 and self.LEGENDS == True:

                    #box = ax.get_position()
                    #ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
                    ncols = int(len(leg_patches)/2)+1
                    # Put a legend below current axis
                    plt.legend(handles=leg_patches,loc='upper center', bbox_to_anchor=(0.5, -0.45),
                            fancybox=True, shadow=True, ncol=ncols,borderpad=1.0,fontsize=self.FONT_SIZE)
                    
                
                fig.align_labels()
                self.message("Aligned Labels on Left.")
                plt.tight_layout()
                plt.savefig(str(self.SCRATCH_FOLDER/"{}_plots_{}.png".format(mapbasename,page)))
                plt.close()
                
                pel = self._current_layout.listElements('PICTURE_ELEMENT')[0]

                if self.TEXT_ELEMENTS != None:
                    for tel in self._current_layout.listElements('TEXT_ELEMENT'):
                        if tel.name in self.TEXT_ELEMENTS:
                            tel.text = self.TEXT_ELEMENTS[tel.name]


                pel.sourceImage = str(self.SCRATCH_FOLDER/"{}_plots_{}.png".format(mapbasename,page))
            self.message(self.total_length_inches / bbox.width)
            
            straight_line_length = self.convert_length_to_inches(self.distance_to(self._left_point,self._right_point))

            self._current_map_frame.camera.scale = straight_line_length / bbox.width
            #self._current_map_frame.camera.scale = self.total_length_inches / bbox.width

    def rotation(self,points,minus=90):
        """
        Calculate the rotation for the map.
        """
        delta = points[-1] - points[0]
        sigma = np.rad2deg(np.arctan2(delta[0],delta[1]))
        return sigma-minus

    def lines_to_points(self):
        """
        Lines to all the points.
        """
        points = []
        for v in self._lines:
            l=v[0]
            if l:
                for part in l:
                    for pnt in part:
                        points.append([pnt.X,pnt.Y])
        return np.array(points)


    def bounding_box(self,pnts):
        """
        Calculate the bounding box to find the upper left side of the lines.
        """
        xmin = np.min(pnts[:,0])
        xmax = np.max(pnts[:,0])
        ymin = np.min(pnts[:,1])
        ymax = np.max(pnts[:,1])
        dy = np.abs(ymax-ymin)
        dx = np.abs(xmax-xmin)
        if dy > dx:
            half = dx*.5+xmin
            halfline = arcpy.Polyline(arcpy.Array([arcpy.Point(half,ymin),arcpy.Point(half,ymax)]))
            ulline = arcpy.Polyline(arcpy.Array([arcpy.Point(xmin,ymax),arcpy.Point(xmax,ymax)]))
        else:
            half = dy*.5+ymin
            halfline = arcpy.Polyline(arcpy.Array([arcpy.Point(xmin,half),arcpy.Point(xmax,half)]))
            ulline = arcpy.Polyline(arcpy.Array([arcpy.Point(xmin,ymin),arcpy.Point(xmin,ymax)]))
        bbox = arcpy.Polyline(arcpy.Array([arcpy.Point(xmin,ymin),arcpy.Point(xmin,ymax),arcpy.Point(xmax,ymax),arcpy.Point(xmax,ymin)]))
        return bbox,halfline,ulline


    def sort_lines(self,ulline):
        """Arrange the lines in order from the left to the right of the rotation. Lines should have no gaps."""
        self._index = index_tools.polyline_kdtree_index(self._lines)
        self._index.build_tree()
        self.orig_index_to_sort_index = {}
        closest = []
        newlines = []
        for i,v in enumerate(self._lines):
            ln=v[0]
            res = self._index.tree.query_ball_point([ln.firstPoint.X,ln.firstPoint.Y],0)
            if len(res)==1:
                closest.append((i,False))
            res = self._index.tree.query_ball_point([ln.lastPoint.X,ln.lastPoint.Y],0)
            if len(res)==1:
                closest.append((i,True))
        mindist = 10000^10
        minindx = -1
        rev = False
        for i,r in closest:
            cp = self._lines[i][0]
            _,_,fpd,_ = ulline.queryPointAndDistance(cp.firstPoint)
            _,_,lpd,_= ulline.queryPointAndDistance(cp.lastPoint)
            if fpd < lpd:
                if fpd<mindist:
                    mindist = fpd
                    minindx = i
                    rev=r
            else:
                if lpd<mindist:
                    mindist = lpd
                    minindx = i
                    rev=r
        start_line = self.copy_line(self._lines[minindx][0],rev)
        newidx = 0
        newlines.append([start_line,0,minindx])
        self.orig_index_to_sort_index[minindx] = newidx
        current_line = start_line
        self.total_length = current_line.length
        current_index = minindx
        while len(newlines)<=len(self._lines):
            indc = self._index.tree.query_ball_point((current_line.lastPoint.X,current_line.lastPoint.Y),r=0)
            matches = self._index.point_to_polyline_index[indc]
            matches = matches[matches!=current_index]
            if len(matches)==0:
                break
            match = matches[0]

            next_line = self._lines[match][0]
            rev = False
            a = np.array([[next_line.lastPoint.X,next_line.lastPoint.Y],[next_line.firstPoint.X,next_line.firstPoint.Y]])
            b = np.array([[current_line.lastPoint.X,current_line.lastPoint.Y]])
            if np.sqrt(np.sum(np.square(a[0]-b))) < np.sqrt(np.sum(np.square(a[1]-b))):
                rev = True
            new_line = self.copy_line(next_line,rev)
            newlines.append([new_line,0,match])
            newidx +=1
            self.orig_index_to_sort_index[match] = newidx
            current_line = new_line
            current_index = match
            self.total_length += current_line.length

        self.total_length_inches = self.convert_length_to_inches(self.total_length)
        self._left_point = newlines[0][0].firstPoint
        self.message(self._left_point)
        self._right_point = newlines[-1][0].lastPoint
        self.message(self._right_point)
        return newlines
        
    def copy_line(self,line,reverse=False):
        """Make a copy of the line. Flip it if reverse == True"""
        parts = arcpy.Array()
        for part in line:
            partlst = []
            for pnt in part:
                partlst.append(pnt)
            if reverse:
                partlst = partlst[::-1]
            parts.append(arcpy.Array(partlst))
        return arcpy.Polyline(parts)

    def convert_length_to_inches(self,length):
        "Figure out the length of the line in inches."
        self.message("Changing scale")
        units = self.SPATIAL_REF.linearUnitName
        self.message(units)
        self.unit_label = "Meters"
        if units == "Meter":
            self.unit_label = "Meters"
            length_inches = length*39.3700787
        elif units == "Foot_US":
            self.unit_label = "Feet"
            length_inches = length*12
        else:
            length_inches = length
        self.message(length_inches)
        return length_inches


    def convert_length_to_miles(self,length):
        "Figure out the length of the line in miles."
        #self.message("Units to miles")
        units = self.SPATIAL_REF.linearUnitName
        #self.message(units)
        self.unit_label = "Meters"
        if units == "Meter":
            self.unit_label = "Meters"
            length_miles = length*0.000621371
        elif units == "Foot_US":
            self.unit_label = "Feet"
            length_miles = length*0.000189394
        else:
            length_miles = length
        #self.message(length_miles)
        return length_miles

    def convert_miles_to_length(self,length):
        "Figure out the length of the line in miles."
        #self.message("Miles to Units")
        units = self.SPATIAL_REF.linearUnitName
        #self.message(units)
        self.unit_label = "Meters"
        if units == "Meter":
            self.unit_label = "Meters"
            length_ = length*1609.34
        elif units == "Foot_US":
            self.unit_label = "Feet"
            length_ = length*5280
        else:
            length_ = length
        #self.message(length_)
        return length_

    def create_stations(self,sorted_lines):
        """Create the station lines."""
        stations_distance = self.STATION_DISTANCES

        stations = []
        distance_tracker = 0
        self.message("start loop")
        stations.append([self.tangent_line(sorted_lines[0][0],0),0])
        for l in sorted_lines:
            linelength = l[0].length
            if self.CONVERT_TO_MILES:
                linelength = self.convert_length_to_miles(linelength)

            distance_tracker += linelength

            if distance_tracker > stations_distance:
                td = distance_tracker - stations_distance
                #self.message(td)
                if self.CONVERT_TO_MILES:
                    td = self.convert_miles_to_length(td)
                #self.message(td)
                stations.append([self.tangent_line(l[0],td),stations_distance])
                stations_distance += self.STATION_DISTANCES
        stations.append([self.tangent_line(sorted_lines[-1][0],sorted_lines[-1][0].length),distance_tracker])
        with arcpy.da.InsertCursor(str(self.stations_fc),["SHAPE@","Station"]) as ic:
            for sta in stations:
                ic.insertRow(sta)


    def tangent_line(self,main_line,position):
        """Calculate a line that is perpindicular at the input line."""
        epsilon = 1e-5
        basepoint = main_line.positionAlongLine(position, False)
        before = main_line.positionAlongLine(position - epsilon, False)
        after = main_line.positionAlongLine(position + epsilon, False)
        dX = after.firstPoint.X - before.firstPoint.X
        dY = after.firstPoint.Y - before.firstPoint.Y
        angle =np.arctan2(dX, dY) * 180 / np.pi
        ll = self.TANGENT_LINE_LENGTH
        #if self.CONVERT_TO_MILES:
        #    ll = self.convert_miles_to_length(ll)
        oneside = basepoint.pointFromAngleAndDistance(angle + 90, ll,"PLANAR")
        secondside = basepoint.pointFromAngleAndDistance(angle - 90, ll,"PLANAR")
        return arcpy.Polyline(arcpy.Array([oneside.firstPoint,secondside.firstPoint]))


    def import_stripchart_layout(self):
        """Description: imports a pagx file and then finds it."""

        current_layouts = [l.name for l in self.MAP_PROJ.listLayouts()]
        self.MAP_PROJ.importDocument(str(self.LAYOUT_TEMPLATE))

        for l in self.MAP_PROJ.listLayouts():
            if l.name not in current_layouts:
                self._current_layout = l
                return True
        return False
    
                
    def update_map(self,rotation):
        self._current_map_frame = self._current_layout.listElements("mapframe_element","StripChartMap")[0]
        self._current_map_frame.map = self.BASE_MAP
        
        self._current_map_frame.camera.heading = rotation
        self._current_map_frame.camera.setExtent(self._current_map_frame.getLayerExtent(self.CENTERLINE_SOURCE, False, True))

                
    def create_plot_from_points(self,pntsource,sorted_lines,ax,selection_distance,color=["Black"],markers=["."],patches=None):
        self.message("Select by distance")
        self.message(selection_distance)
        #sel_pnts = arcpy.management.SelectLayerByLocation(str(pntsource), 'WITHIN_A_DISTANCE', self.CENTERLINE_SOURCE, "{} {}".format(wad,self.unit_label.upper())).getOutput(0)
        sel_pnts = arcpy.management.SelectLayerByLocation(str(pntsource), 'WITHIN_A_DISTANCE', self.CENTERLINE_SOURCE,search_distance=selection_distance).getOutput(0)
        self.message(arcpy.GetCount_management(sel_pnts).getOutput(0))
        with arcpy.da.SearchCursor(sel_pnts,["SHAPE@"]) as sc:
            point_chart = {}
            for row in sc:
                if row[0]:
                    point_chart_row = [row[0].firstPoint]
                    npolyindx = self._index.nearest_polyline_single_point((row[0].firstPoint.X,row[0].firstPoint.Y))
                    sortedIndex = self.orig_index_to_sort_index[npolyindx]
                    npoly = sorted_lines[sortedIndex][0]
                    res = npoly.queryPointAndDistance(row[0].firstPoint)
                    point_chart_row.append(res[0])
                    point_chart_row.append(res[1])
                    if self.CONVERT_TO_MILES:
                        miles = self.convert_length_to_miles(res[2])
                        point_chart_row.append(miles)
                    else:
                        point_chart_row.append(res[2])
                    
                    try:
                        point_chart_row.append(color[npolyindx])
                    except:
                        point_chart_row.append(color[0])
                    try:
                        point_chart_row.append(markers[npolyindx])
                    except:
                        point_chart_row.append(markers[0])

                    point_chart[sortedIndex] = point_chart_row

        xdist=0
        xs=[]
        ys=[]
        cs = []
        ms = []
        for i,x in enumerate(sorted_lines):
            cl = x[0]
            if i in point_chart:
                info = point_chart[i]
                xval = xdist
                if self.CONVERT_TO_MILES:
                    xval = self.convert_length_to_miles(xdist)
                xs.append(xval)
                ys.append(info[3])
                cs.append(info[4])
                ms.append(info[5])

            xdist += cl.length
        self.message("point chart")
        self.message(ys)
        self.message(cs)
        for x,y,c,m in zip(xs,ys,cs,ms):
            ax.scatter(x,y,color=c,marker=m)
        if patches:
            ax.legend(handles=patches,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    def create_plot_from_polygons(self,psource,sorted_lines,ax,fname,color=["Black"],patches=None):
        wad = self.TANGENT_LINE_LENGTH*.1
        self.message("Select by distance")
        self.message(wad)
        #sel_pnts = arcpy.management.SelectLayerByLocation(str(pntsource), 'WITHIN_A_DISTANCE', self.CENTERLINE_SOURCE, "{} {}".format(wad,self.unit_label.upper())).getOutput(0)
        sel_polys = arcpy.management.SelectLayerByLocation(str(psource), 'WITHIN_A_DISTANCE', self.CENTERLINE_SOURCE,search_distance=wad).getOutput(0)
        self.message(arcpy.GetCount_management(sel_polys).getOutput(0))
        if arcpy.GetCount_management(sel_polys).getOutput(0) == '0':
            return False
        current_length = 0
        ranges = []


        polygons = []
        with arcpy.da.SearchCursor(sel_polys,["SHAPE@",fname]) as sc:
            for row in sc:
                if row[0]:
                    polygons.append([row[0].buffer(wad),row[1]])
        indices = []
        while len(polygons)>0:
            cp = polygons.pop()
            poly = cp[0]
            lbl = cp[1]
            pieces = []        
            for i,d in enumerate(sorted_lines):
                line = d[0]
                if line.disjoint(poly)==False:
                    if i not in indices:
                        pieces.append(i)
                        indices.append(i)
            if len(pieces)>0:
                pieces.sort()
                self.message(pieces)
                start_length = 0
                for i in range(0,pieces[0]):
                    line = sorted_lines[i][0]
                    start_length += line.length
                end_length = start_length
                for i in pieces:
                    line = sorted_lines[i][0]
                    end_length += line.length

                if self.CONVERT_TO_MILES:
                    start_length = self.convert_length_to_miles(start_length)
                    end_length = self.convert_length_to_miles(end_length)
                ranges.append([start_length,end_length,color,lbl])

        ranges = sorted(ranges,key=lambda x:x[0])
        if len(ranges)!=0:
            self.message(ranges)
            ax.set_ylim(0,1)
            ax.set_yticks([])
            intervals = .5/(len(ranges)*.5)
            y_vals = list(np.arange(.08,.5,intervals))+list(np.arange(.6,.02,-intervals))
            for span,yt in zip(ranges,y_vals):
                spanLabel = span[3]
                spanColor = span[2]
                #spanLabel = spanLabel.replace(" ","\n")
                half = (span[1]-span[0])*.5
                ax.axvspan(span[0], span[1], color=spanColor, alpha=0.4)
                ax.text(half+span[0],yt,spanLabel,fontsize='xx-small',ha="center",va="center")
        else:
            ax.set_yticks([])





    def distance_to(self,pnt1,pnt2):
        pnt1 = np.array([pnt1.X,pnt1.Y])
        pnt2 = np.array([pnt2.X,pnt2.Y])
        delta = pnt1-pnt2
        return np.sqrt(np.sum(delta*delta))



    def create_att_line_plots(self,values,sorted_lines,ax,color="Black",add_yticklabels=False):
        xdist=0
        xs=[]
        ys=[]
        for x in sorted_lines:
            cl = x[0]
            xval = xdist
            xvalO = xdist+cl.length
            if self.CONVERT_TO_MILES:
                xval = self.convert_length_to_miles(xdist)
                xvalO = self.convert_length_to_miles(xvalO)
            xs.append(xval)
            xs.append(xvalO)
            ys.append(values[x[2]])
            ys.append(values[x[2]])
            xdist += cl.length

        ax.plot(xs,ys,color=color)
        if add_yticklabels:
            uys = sorted(list(set(ys)))
            #ax.set_yticklabels(uys)
        return


    def set_font_size(self,ax:plt.axes,fs=9):
        """
        Description:
            Set the font size for the title, axis label (x and y), ticklabels (x and y) to the specified font size
        Args:
            ax (pyplot.axes) plot axes
            fs (int): font size
        Returns:
            None"""
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(fs)

