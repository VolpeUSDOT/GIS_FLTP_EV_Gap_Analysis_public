
from pathlib import Path

import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rcParams

import arcpy
import arcgis
from arcgis.features import GeoAccessor, GeoSeriesAccessor
from arcgis.geometry import Point, Polyline
import math

"""
This script creates custom tables with histograms of network distances. 

Create Date: 2022
Created By: The Volpe Center's Geospatial Team 
"""

def create_visitation_dictionary(visit_table:Path,visit_fc:Path)->dict:
    """
    parameters:
        visit_table (Path): path to geodatabase table of visits
        visit_fc (Path): path to geodatabase feature class of all centroids
    returns:
        dictionary of visitation data and labels by matchid\sites
    """
    visits_dict = {}

    with arcpy.da.SearchCursor(str(visit_table),["matchsites","label","visitation"]) as sc:
        visits_dict = {row[0]:{"label":row[1],"visitation":row[2]} for row in sc}

    with arcpy.da.SearchCursor(str(visit_fc),["matchid","label","visitation"]) as sc:
        for row in sc:
            if row[0] not in visits_dict:
                visits_dict[row[0]] = {"label":row[1],"visitation":row[2]}
            #if row[0] in visits_dict:
            #    visits_dict[row[0]]["postal"] = row[3]

    return visits_dict


def get_max_miles_by_route(od_gdb:Path)->dict:
    """
    Parameters:
        od_gdb:Path - path to od-database
    Returns:
        dictionary of route name keys and their max miles
    """
    max_miles={}
    arcpy.env.workspace = str(od_gdb)
    for fc in arcpy.ListFeatureClasses():
        if "charging_stations" not in fc:
            mm = max([row[0] for row in arcpy.da.SearchCursor(fc,["Total_Miles"])])
            max_miles[fc] = mm
    return max_miles        
    

def get_max_miles_route_by_agency_state(routes_table:Path,max_miles:dict,route_gdb:Path,visits_dict:dict):
    """
    Parameters:
        route_tables:Path - path ot route table
        max_miles:dict - dictionary of route names and max miles
    Returns:
        two dictionaries, route names (key is agency, value is list of org, agency, routename, and matchid) and max route by agency (key is agency, value is max miles)
    """

    route_names_dict = {}
    max_routes = {}
    #states_pth = admin_data["state"][0]
    #states = pd.DataFrame.spatial.from_featureclass(states_pth)
    #sind = states.spatial.sindex()

    #for k in admin_data.keys():
        #if k != "state":
            #print(k)
            #arcpy.AddSpatialIndex_management(str(admin_data[k][0]))
            #admin_data[k].append(pd.DataFrame.spatial.from_featureclass(admin_data[k][0]))
            #print("index")
            #_ = admin_data[k][-1].spatial.sindex()
    print("state search")
    processed = []
    duplicates = []

    with arcpy.da.SearchCursor(str(routes_table),["org","agency","routefcid","matchid","labeltxt","postal"],sql_clause=(None, "ORDER BY org")) as sc:
        for row in sc:
            #agency
            #get visitation information
            if row[2]:
                routes =  pd.DataFrame.spatial.from_featureclass(route_gdb/row[2])
                #routes["sp"] = routes["SHAPE"].apply(lambda x: Point({"x":x.coordinates()[0][0][0],"y":x.coordinates()[0][0][1],"spatialReference" : {"wkid" : 4326}}))
                #routes["SHAPE"] = routes["sp"]
                #routes.spatial.set_geometry("SHAPE")
                #routes = routes.reset_index()
                #stvs = states.spatial.select(routes)["postal"].values
                
                #adminvs = admin_data[row[1]][2].spatial.select(routes)[admin_data[row[1]][1]].values
                #if len(stvs)>0:
                #    st = stvs[0]
                #else:
                #    st = ""
                #if len(adminvs)>0:
                    #adminlbl = adminvs[0]
                #else:
                    #adminlbl = ""
                #del routes
                label = row[4]
                
                try:
                    label = visits_dict[row[3]]["label"]
                except:
                    if row[4]:
                        label = row[4].title()
                    else:
                        label = None


                if label != None:
                    if label not in processed:
                        processed.append(label)
                    else:
                        duplicates.append(label)
                        #label = label + "**"

                    if row[1] in route_names_dict:

                        if row[5] in route_names_dict[row[1]]:
                            route_names_dict[row[1]][row[5]].append([row[0],row[1],row[2],row[3],label])
                        else:
                            route_names_dict[row[1]][row[5]] = [[row[0],row[1],row[2],row[3],label]]
                    else:

                        route_names_dict[row[1]] = {row[5]:[[row[0],row[1],row[2],row[3],label]]}
                    mm = max_miles[row[2]]
                    if row[1] in max_routes:
                        if row[5] in max_routes[row[1]]:
                            if mm > max_routes[row[1]][row[5]]:
                                max_routes[row[1]][row[5]] = int(mm)
                        else:
                            max_routes[row[1]][row[5]] = int(mm)
                    else:
                        max_routes[row[1]] = {row[5]:int(mm)}

    for agency,sd in route_names_dict.items():
        for state,values in sd.items():
            #values = route_names_dict[agency][state]
            try:
                values = sorted(values,key=lambda x: x[4])
                for v in values:
                    if v[4] in duplicates:
                        if v[4]!=None:
                            v[4] = v[4] + "**"
                route_names_dict[agency][state] = values
            except:
                pass
    print("returning")

    return route_names_dict,max_routes


def get_max_miles_route_by_agency(routes_table:Path,max_miles:dict):
    """
    Parameters:
        route_tables:Path - path ot route table
        max_miles:dict - dictionary of route names and max miles
    Returns:
        two dictionaries, route names (key is agency, value is list of org, agency, routename, and matchid) and max route by agency (key is agency, value is max miles)
    """
    route_names_dict = {}
    max_routes = {}
    with arcpy.da.SearchCursor(str(routes_table),["org","agency","routefcid","matchid","label","state"],sql_clause=(None, "ORDER BY label")) as sc:
        for row in sc:
            if row[1] in route_names_dict:
                route_names_dict[row[1]].append([row[0],row[1],row[2],row[3],row[4],row[5]])
            else:
                route_names_dict[row[1]] = [[row[0],row[1],row[2],row[3],row[4],row[5]]]
            mm = max_miles[row[2]]
            if row[1] in max_routes:
                if mm > max_routes[row[1]]:
                    max_routes[row[1]] = int(mm)
            else:
                max_routes[row[1]] = int(mm)
    return route_names_dict,max_routes

def check_width(ax,text,renderer,props={"ha":"center","va":"center", "size":13,"color":"k"},step=0):
    txtbox = ax.text(0.5, 0.5, text,ha=props["ha"], va=props["va"], size=props["size"],color=props["color"])

    tw = txtbox.get_window_extent(renderer=renderer).width
    axw = ax.get_window_extent(renderer=renderer).width
    if tw > axw-.1:
        txtbox.remove()
        t = text.split(" ")
        if step == 3:
            txtbox = ax.text(0.5, 0.5, text,ha=props["ha"], va=props["va"], size=int(props["size"]*.7),color=props["color"])
            return None
        stop = math.floor(len(t)/2)
        newtext = " ".join(t[:stop]) + "\n" + " ".join(t[stop:])

        step+=1
        check_width(ax,newtext,renderer,props,step)
    else:
        return None

def create_page_fig(ncols,nrows,figsize,gridspec,titles):
    
    fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=figsize,gridspec_kw=gridspec,dpi=300)
    fig.subplots_adjust(0.05,0.05,0.95,0.95, wspace=0.00, hspace=.05)
    fig.subplots_adjust()
    r = fig.canvas.get_renderer()
    for row,iaxes in enumerate(axes):
        bc = "k"
        if row == 0:
            spinetop = True
            spinebot = True
        elif row == (len(axes)-1):
            spinetop = False
            spinebot = True
        else:
            spinetop = False
            spinebot = False
        for axi,ax in enumerate(iaxes):
            #if axi!=3:
            ax.tick_params(labelbottom=0, labelleft=0, bottom=0, top=0, left=0, right=0)
            ax.ticklabel_format(useOffset=False, style="plain")
            for d,s in ax.spines.items():
                #print(d)
                #s.set_visible(False)
                if d in ['right','left']:
                    s.set_visible(False)
                if d in ['top']:
                    s.set_color(bc)
                    s.set_visible(spinetop)
                if d in ['bottom']:
                    s.set_color(bc)
                    s.set_visible(spinebot)
        if row == 0:
            for t_i, t in enumerate(titles):
                check_width(iaxes[t_i],t,r)


    return fig,axes

def histplot(ax,df,avg,avgcolor,minstr,maxstr,cnts,binrange,binwidth):
    ax.margins(0.1)
    try:
        sns.histplot(df["Total_Miles"],ax=ax,kde=False,color="#994455",binwidth=binwidth,binrange=binrange)
    except:
        sns.histplot(df["Total_Miles"],ax=ax,kde=False,color="#994455",binwidth=binwidth)
    ylim = ax.get_ylim()
    xlim = ax.get_xlim()
    ax.set_ylim((0,ylim[1]))
    #ax.set_xlim((xlim[0],xlim[1]+1))
    #ax.set_xlim((0,xlim[1]))
    #sns.violinplot(x=sdf_dcfc["Total_Miles"],ax=axes[5],color="#BEE8FF",inner="quartile")
    #sns.boxplot(x=sdf_dcfc["Total_Miles"],ax=axes[5],color="#BEE8FF")
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.axvline(x=avg,color=avgcolor,linewidth="1")
    if REFERENCE_LINE < binrange[1]:
        ax.axvline(x=REFERENCE_LINE,color="k",linestyle="--",linewidth="1")
    ax.axhline(y=0,color="k")
    ax.axvline(x=0,color="k")
    ax.text(.98, .1, binrange[1],transform=ax.transAxes,ha="right", va="center", size=9)
    ax.text(.98, .92, f"Lots: {cnts}*",transform=ax.transAxes,ha="right", va="center", size=9)
    #ax.text(.15, .1, binrange[0],transform=ax.transAxes,ha="left", va="center", size=9)
    ax.text(.05, .1, "0",transform=ax.transAxes,ha="right", va="center", size=9)
    ax.text(.05, .9, int(ylim[1]),transform=ax.transAxes,ha="right", va="center", size=9)
    #ax.text(50, 1, 50,ha="left", va="center", size=9)


def process_routes(od_fgdb:Path,visits_dict:dict,sdf_sites,sdf_altfuel,states:dict,agency:str,mm:int,outputFolder:Path,numCols:int,numRows:int,figSize:tuple,gs:dict,titles:list,adddraft=True,pdfname="table",processWithoutVisits=True):
    pdfPath = outputFolder / f'{pdfname}_{agency}.pdf'
    csvPath = outputFolder / f'{pdfname}_{agency}.csv'
    state_order = ['AK', 'AL', 'AR', 'AZ', 'CA', 'CO', 'CT', 'DC', 'DE', 'FL', 'GA', 'HI', 'IA', 'ID', 'IL', 'IN', 'KS', 'KY',
         'LA', 'MA', 'MD', 'ME', 'MI', 'MN', 'MO', 'MS', 'MT', 'NC', 'ND', 'NE', 'NH', 'NJ', 'NM', 'NV', 'NY',
         'OH', 'OK', 'OR', 'PA', 'RI', 'SC', 'SD', 'TN', 'TX', 'UT', 'VA', 'VT', 'WA', 'WI', 'WV', 'WY']
    print(str(pdfPath))
    page_num = 1
    heading = "{0}\nDistance from Parking Lots to Nearest Stations".format(FULL_NAMES[agency])
    output_table_list = []
    with PdfPages(pdfPath) as pdf:
        row_counter = 1
        fig,axes = create_page_fig(numCols,numRows,figSize,gs,titles)
        r = fig.canvas.get_renderer()
        bbox = axes[0][4].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        width = bbox.width
        print(width)
        plt.text(0.5,0.01,page_num, transform=fig.transFigure, ha="center", va="center", size=10)
        plt.text(0.08,0.035,"*May not include all lots.\n**Visitation for multiple states.", transform=fig.transFigure, ha="left", va="center", size=8)
        plt.text(0.5,0.97,heading, transform=fig.transFigure, ha="center", va="center", size=15)
        if adddraft:
            plt.text(0.5,0.5,"DRAFT", transform=fig.transFigure, ha="center", va="center", size=72,color=(1,0,0,.2),rotation=45)
        for k in state_order:
            if k in list(states.keys()):
                routes = states[k]
                
                for v in routes:
                    print(v)
                    
                    if v[2]:
                        table_list = [v[4],k,]
                        ex_fc = od_fgdb / v[2]
                        visitation_raw = -1
                        if processWithoutVisits == True:
                            visitation_raw = -2
                        site_name = v[4] #v[0].replace("''","'")
                        state_lbl = k
                        try:
                            visitation_raw = visits_dict[v[3]]["visitation"]
                            #site_name = visits_dict[v[3]]["label"]
                            
                        except:
                            pass
                        try:
                            state_lbl = k#visits_dict[v[3]]["postal"]
                        except:
                            pass
                        
                        try:
                            visitation_raw = float(visitation_raw)
                        except:
                            if processWithoutVisits == True:
                                visitation_raw = -2
                            else:
                                visitation_raw = -1
                        

                        if visitation_raw!=-1:
                            table_list.append(visitation_raw)
                            sdf_dcfc = pd.DataFrame.spatial.from_featureclass(ex_fc)
                            cnts = len(sdf_dcfc)
                            sdf_dcfc.spatial.project(sdf_altfuel.spatial.sr)
                            sdf_sites_temp = sdf_sites[(sdf_sites["org"]==v[0].replace("''","'"))&(sdf_sites["agency"]==v[1])].copy().reset_index()

                            selected = sdf_altfuel.spatial.select(sdf_sites_temp)
                            sdf_sites_temp = sdf_sites[(sdf_sites["org"]==v[0].replace("''","'"))&(sdf_sites["agency"]==v[1])].copy().reset_index()
                            
                            selected = sdf_sites_temp.spatial.select(sdf_altfuel)

                            p = len(selected) / len(sdf_sites_temp)*100
                            if p > 98:
                                
                                pie_data = [p,100-p]
                                pie_keys = ["Within","Outside"]
                                pie_colors = ["#007FBF","#C04000"]
                                #axes[row_counter][4].pie(pie_data,textprops={'fontsize': 8},colors=pie_colors)
                                axes[row_counter][5].text(0.5, 0.5, "100%",ha="center", va="center", size=13,color="#007FBF")
                            elif p < 1:
                                
                                pie_data = [p,100-p]
                                pie_keys = ["Within","Outside"]
                                pie_colors = ["#007FBF","#C04000"]
                                #axes[row_counter][4].pie(pie_data,textprops={'fontsize': 8},colors=pie_colors)
                                axes[row_counter][5].text(0.5, 0.5, "0%",ha="center", va="center", size=13,color="#C04000")
                            else:
                                pie_data = [p,100-p]
                                pie_keys = ["Within","Outside"]
                                pie_colors = ["#007FBF","#C04000"]
                                axes[row_counter][5].pie(pie_data,textprops={'fontsize': 8},colors=pie_colors, autopct='%.0f%%')
                            #axes[row_counter][4].text(0.5, 0.5, "{:,.0f}%".format(p),ha="center", va="center", size=10,color="#210203")
                            

                            
                            del sdf_sites_temp
                            average_dcfc_value =sdf_dcfc["Total_Miles"].mean()
                            table_list.append(average_dcfc_value)
                            average_dcfc_value_str = "{:,.1f}".format(sdf_dcfc["Total_Miles"].mean())
                            max_dcfc_value_str = "{:,.0f}".format(sdf_dcfc["Total_Miles"].max())
                            table_list.append(sdf_dcfc["Total_Miles"].max())
                            min_dcfc_value_str = "{:,.0f}".format(sdf_dcfc["Total_Miles"].min())
                            table_list.append(sdf_dcfc["Total_Miles"].min())
                            table_list.append(p)
                            if visitation_raw >0:
                                visitation_value_str = "{:,.0f}".format(visitation_raw)
                            else:
                                visitation_value_str = " "
                                site_name = site_name.replace("**","")
                            #axes[row_counter][0].text(0.5, 0.5, site_name, ha="center", va="center", size=10)
                            check_width(axes[row_counter][0],site_name,r,{"ha":"center","va":"center", "size":10,"color":"k"}) 
                            axes[row_counter][1].text(0.5, 0.5, state_lbl,ha="center", va="center", size=10)
                            axes[row_counter][2].text(0.5, 0.5, visitation_value_str,  ha="center", va="center", size=10)
                            avgcolor = "green"
                            if average_dcfc_value >REFERENCE_LINE:
                                axes[row_counter][3].text(0.5, 0.5, average_dcfc_value_str,ha="center", va="center", size=13,color="white")
                                axes[row_counter][3].set_facecolor("red")
                                avgcolor = "red"
                            else:
                                axes[row_counter][3].text(0.5, 0.5, average_dcfc_value_str,ha="center", va="center", size=13,color="white")
                                axes[row_counter][3].set_facecolor("green")
                                avgcolor = "green"
                            binwidth = .1/width * (mm[k]+1)
                            binrng = (0,mm[k]+1)
                            print(binwidth)
                            print()

                            if mm[k] ==0:
                                #sdf_dcfc["Total_Miles"].max()
                                binrng = (0,sdf_dcfc["Total_Miles"].max()+1)
                                binwidth = .1 / width * (sdf_dcfc["Total_Miles"].max()+1)
                            histplot(axes[row_counter][4],sdf_dcfc,average_dcfc_value,avgcolor,min_dcfc_value_str,max_dcfc_value_str,cnts,binrng,binwidth)
                            if site_name == "Sand Creek Massacre NHS":
                                print("breaking")
                            row_counter +=1
                            if row_counter >= NUM_ROWS:
                                #plt.tight_layout()
                                if fig:
                                    pdf.savefig()  # saves the current figure into a pdf page
                                    plt.close()
                                    #break
                                    row_counter = 1
                                    page_num+=1
                                    #if page_num >5:
                                        #break
                                    fig,axes = create_page_fig(NUM_COLS,NUM_ROWS,FIG_SIZE,GRID_SPEC,COL_TITLES)
                                    plt.text(0.5,0.01,page_num, transform=fig.transFigure, ha="center", va="center", size=10)
                                    plt.text(0.08,0.035,"*May not include all lots.\n**Visitation for multiple states.", transform=fig.transFigure, ha="left", va="center", size=8)
                                    plt.text(0.5,0.97,heading, transform=fig.transFigure, ha="center", va="center", size=15)
                                    if adddraft:
                                        plt.text(0.5,0.5,"DRAFT", transform=fig.transFigure, ha="center", va="center", size=72,color=(1,0,0,.2),rotation=45)
                                
                            output_table_list.append(table_list)
                    # else:
                        
                    #     visitation_raw = -1
                    #     site_name = v[0]
                    #     state_lbl = k
                    #     try:
                    #         visitation_raw = visits_dict[v[4]]["visitation"]
                    #         site_name = visits_dict[v[4]]["label"]
                    #     except:
                    #         pass
                        
                    #     axes[row_counter][0].text(0.5, 0.5, site_name, ha="center", va="center", size=10)
                    #     if visitation_raw >0:
                    #         visitation_value_str = "{:,.0f}".format(visitation_raw)
                    #     else:
                    #         visitation_value_str = " "
                    #     axes[row_counter][1].text(0.5, 0.5, state_lbl,r,{"ha":"center","va":"center", "size":10,"color":"k"})
                    #     axes[row_counter][2].text(0.5, 0.5, " ",ha="center", va="center", size=13,color="white")
                    #     axes[row_counter][2].set_facecolor("red")
                    #     avgcolor = "red"
                        
                        





                    
        try:
            pdf.savefig()
            plt.close()
        except:
            pass
    df = pd.DataFrame(output_table_list,columns=["LABEL","STATE","VISITATION","AVERAGE","MAX","MIN","PERCENTAGE"])
    df.to_csv(csvPath)


OUTPUTS = Path("C:\\Users\\David.Lamb\\OneDrive - DOT OST\\Documents\\GitHub\\GIS_FLTP_EV_Gap_Analysis\\data\\output")
INPUTS = Path("C:\\Users\\David.Lamb\\OneDrive - DOT OST\\Documents\\GitHub\\GIS_FLTP_EV_Gap_Analysis\\data\\input")
GDB_PATH = OUTPUTS / "od_pairs_dcfc_all.gdb"
ALT_FUEL = INPUTS / "ALTFUEL" / "Alternative_Fuel_Corridors.shp"
VISITATION_TABLE = OUTPUTS / "ev_gap_analysis.gdb" / "visitation_table"
ROUTES_TABLE = GDB_PATH  / "org_agency_route"
VISITATION_FC = OUTPUTS / "ev_gap_analysis.gdb" / "flma_points_visits_final"
SITES_FC = OUTPUTS / "ev_gap_analysis.gdb" / "flma_points_prj_charge_stations"

REFERENCE_LINE =   50 #10 #


#ADMIN UNITS
STATE_BOUNDARIES = INPUTS / "Mapping_Info.gdb" / "all_states_natural_earth"
# BLM_UNITS = INPUTS / "VISITATION" / "DOI_BLM" / "admu.gdb" / "blm_natl_admu_field_poly_webpub"
# BLM_LABEL = "ADMU_NAME"
# BOR_UNITS = INPUTS / "VISITATION" / "DOI_BOR" / "Area_Office_Boundaries.shp"
# BOR_LABEL = "AreaOffice"
# FWS_UNITS = INPUTS / "VISITATION" / "DOI_FWS" / "FWSBoundaries.shp"
# FWS_LABEL = "ORGNAME"
# NPS_UNITS = INPUTS / "VISITATION" / "DOI_NPS" / "NPS_-_Land_Resources_Division_Boundary_and_Tract_Data_Service.shp"
# NPS_LABEL = "UNIT_NAME"
# FS_UNITS = INPUTS / "VISITATION" / "USDA_FS" / "AdministrativeForest.shp"
# FS_LABEL = "FORESTNAME"

# ADMIN_UNITS = {"state":[STATE_BOUNDARIES],'doi_blm':[BLM_UNITS,BLM_LABEL], 'doi_bor':[BOR_UNITS,BOR_LABEL], 'doi_fws':[FWS_UNITS,FWS_LABEL],
# 'doi_nps':[NPS_UNITS,NPS_LABEL], 'usda_fs':[FS_UNITS,FS_LABEL], 'usace':[]}

FULL_NAMES = {'doi_blm':'Bureau of Land Management', 'doi_bor':'Bureau of Reclamation', 'doi_fws':'Fish and Wildlife Service',
'doi_nps':'National Park Service', 'usda_fs':'US Forest Service', 'usace':'US Army Corps of Engineers'}

#FIGURE SETTINGS
NUM_ROWS = 12
NUM_COLS = 6
GRID_SPEC = {"width_ratios":[1.8,.8,.8,1,2.8,.8],'height_ratios':[1,.8,.8,.8,.8,.8,.8,.8,.8,.8,.8,.8,.8,.8,.8]}
GRID_SPEC = {"width_ratios":[1.8,.8,.8,1,2.8,.8],'height_ratios':[1,.8,.8,.8,.8,.8,.8,.8,.8,.8,.8,.8]}
COL_TITLES = ["Name","State","Visitation","Average Miles to Existing DCFC","Miles to Existing DCFC","% Near Ready or Pending Corridor"]
FIG_SIZE = (8.5,11)
AGENCIES = ['doi_blm', 'doi_bor', 'doi_fws', 'doi_nps', 'usda_fs', 'usace']
#AGENCIES = ['doi_blm', 'doi_bor']
AGENCIES = ['usace']




if __name__ == "__main__":
    # Set the style properties for the charts
    sns.set() #imports seaborn style for use with plt
    sns.set_style("ticks")
    rcParams['font.family'] = 'serif'
    rcParams['font.serif'] = ['Gill Sans MT']
    print("visitsDict")
    visits_dict = create_visitation_dictionary(VISITATION_TABLE,VISITATION_FC)
    print("maxmiles")
    max_miles = get_max_miles_by_route(GDB_PATH)
    print("routenames")
    print(ROUTES_TABLE)
    print(GDB_PATH)
    route_names_dict,max_routes = get_max_miles_route_by_agency_state(ROUTES_TABLE,max_miles,GDB_PATH,visits_dict)
    
    sdf_sites = pd.DataFrame.spatial.from_featureclass(SITES_FC)
    sdf_altfuel = pd.DataFrame.spatial.from_featureclass(ALT_FUEL)
    evtypes = ['Signage Ready', 'Signage Pending', 'Signage  Pending']
    sdf_altfuel = sdf_altfuel[sdf_altfuel["ELECTRICVE"].isin(evtypes)].copy().reset_index()
    polygons=sdf_altfuel["SHAPE"].apply(lambda x: x.buffer(80467.2))
    sdf_altfuel["SHAPE"] = polygons
    sdf_altfuel.spatial.set_geometry("SHAPE")
    #build fuel index now
    sindex = sdf_altfuel.spatial.sindex()
    for k,routes in route_names_dict.items():
        print(k)
        mm = max_routes[k]
        if k in AGENCIES:
            process_routes(GDB_PATH,visits_dict,sdf_sites,sdf_altfuel,routes,k,mm,OUTPUTS,NUM_COLS,NUM_ROWS,FIG_SIZE,GRID_SPEC,COL_TITLES,pdfname="DCFCtable_all",adddraft=False)
        #without draft
        #process_routes(GDB_PATH,visits_dict,sdf_sites,sdf_altfuel,routes,k,mm,OUTPUTS,NUM_COLS,NUM_ROWS,FIG_SIZE,GRID_SPEC,COL_TITLES,False)
    
##Placeholders for future implementation of data structure


