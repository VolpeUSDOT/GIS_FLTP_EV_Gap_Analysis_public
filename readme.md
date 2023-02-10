# EV Gap Analysis

## Data Structure

```
-data**
|
|    -input
|    |
|    |---- ALTFUEL -> Alternative Fuel Corridors
|    |    |
|    |    |---- (Alternative_Fuel_Corridors.shp)[https://www.fhwa.dot.gov/environment/alternative_fuel_corridors/]
|    |    |---- Used for mapping
|    |    |
|    |---- DOI_BLM -> BLM Parking locations
|    |    |
|    |    |---- (DOI_BLM\BLM_National_Recreation_Site_pts\recs.gdb)[https://www.landscape.blm.gov/geoportal/catalog/search/browse/browse.page]
|    |    |---- \RECS\recs_pt <- Points for parking
|    |    |
|    |---- DOI_BOR -> FLMA Parking locations
|    |    |
|    |    |---- Jan_2021_FLMA_Submittal_UPDATE2021June22.gdb
|    |    |---- \ BOR_2021 <- FLTP Road network
|    |    |---- \ NPS_PAVED_2021 <- FLTP Road network Supplement for NPS Parking
|    |    |
|    |---- DOI_FWS -> Fish and Wildlife Parking locations
|    |    |
|    |    |---- JFWS_Public_Parking_Lot_All\FWS_Public_Parking_Lot_All.shp <- Used for parking lots
|    |    |---- FWS_NWRS_HQ_VisitorServiceAmenities_V2\FWS_NWRS_HQ_VisitorServiceAmenities_V2.shp <- Needed to match up to the visitation tables
|    |    |
|    |---- DOI_NPS -> National Parks Parking locations
|    |    |
|    |    |---- (NPS_-_Parking_Areas_-_Web_Mercator\NPS_-_Parking_Areas_-_Web_Mercator.shp)[https://public-nps.opendata.arcgis.com/datasets/nps-parking-areas-web-mercator-1/explore?location=38.484649%2C-120.417298%2C10.98]
|    |    |----  FLTP Road network Supplement for NPS Parking (see the DOI_BOR \ Jan_2021_FLMA_Submittal)
|    |    |
|    |---- templates -> files for creating strip charts
|    |    |
|    |    |---- ArcGIS Project, layout pagx and layer lyrx
|    |    |
|    |---- EV_stations -> use the API
|    |    |
|    |    |---- (Alternative Fuel Stations)[https://developer.nrel.gov/docs/transportation/alt-fuel-stations-v1/]
|    |    |
|    |---- USDA_FS -> US Forest Service Parking locations
|    |    |
|    |    |---- FS_Points_of_Interest_mod.xlsx <- Converted to points
|    |    |
|    |---- USACE -> US Army Corps of Engineers Parking locations
|    |    |
|    |    |---- USACE_Recreation_Project_Site_Areas_(Points).shp
|    |    |
|    |---- VISITATION -> Visitation Data Spreadsheets
|    |    |
|    |    |----DOI_BLM
|    |    |    |----(admu.gdb)[https://gbp-blm-egis.hub.arcgis.com/datasets/blm-national-administrative-unit-boundary-polygons-and-office-points-national-geospatial-data-asset-ngda-1/about] -> for the centroids
|    |    |    |----Visitation_Cleaned_bor.xlsx <- From agency, converted from pdf to spreadsheet.
|    |    |
|    |    |----DOI_BOR
|    |    |    |----Area_Office_Boundaries.shp -> for the centroids
|    |    |    |----Visitation_Cleaned_blm.xlsx <- From agency
|    |    |
|    |    |----DOI_FWS
|    |    |    |----(FWSBoundaries.shp)[https://gis-fws.opendata.arcgis.com/datasets/fws::fws-national-realty-tracts/about] -> for the centroids
|    |    |    |----Visitation_Cleaned_fws.xlsx <- From agency
|    |    |
|    |    |----DOI_NPS
|    |    |    |----(NPS_-_Land_Resources_Division_Boundary_and_Tract_Data_Service.shp)[https://public-nps.opendata.arcgis.com/search?collection=Dataset&q=boundaries] -> for the centroids
|    |    |    |----Visitation_Cleaned_nps.xlsx <- From agency
|    |    |
|    |    |----USDA_FS
|    |    |    |----(AdministrativeForest.shp)[https://data-usfs.hub.arcgis.com/datasets/usfs::forest-administrative-boundaries-feature-layer/about] -> for the centroids
|    |    |    |----Visitation_Cleaned_usfs.xlsx <- From agency
|    |---- Mapping_Info.gdb -> Custom, used for mapping
|
|    -output
|    |
|    |---- ev_gap_analysis.gdb
|    |---- od_pairs_dcfc.gdb
```

```
**Names of files may have changed from the source.
```
# Documentation


[Information on the scripts](docs/readme.md)

## Notes

Data files are not provided in this repository. The source of the data is noted above when available.
