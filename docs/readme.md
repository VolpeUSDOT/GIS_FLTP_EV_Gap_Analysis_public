
# Notes
1. Run import_data.py 
    - loads data from multiple sources, combines visitation data with the site locations
    - Need api key from https://developer.nrel.gov/signup/ for latest charging stations
    - Change key str variable
    - This will create ev_gap_analysis.gdb
2. Run quick_search.py
    - Finds DCFC stations with 50 miles and Level 2 stations within 10 miles.
3. Run closest_facility_nax.py
    - Finds routes to closest stations. Creates od_pairs_dcfc.gdb (or the name set by the parameters)
4. Run generate_table_plots.py
    - Creates tables with plots for each agency.
5. Run stripcharts.py
    - Modify stripcarts_data.py with the data.
    - Creates the ArcGIS projcets for the strip charts.
    
## Maps

### State level charging stations

After running steps 1 and 2 above, to create the maps use the EV_GAP_Map project. You will need the ExportByAgency.ipynb and ExportByState.ipynb notebooks to export.

Run ExportByState to create the state level maps showing all sites.
Run ExportByAgency to create the state level pdf maps showing individual agency sites.

