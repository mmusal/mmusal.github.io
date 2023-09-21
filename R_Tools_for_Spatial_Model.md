Both R and Python have libraries that can process files referred to as "shape files". Details on the definition and development of the format can be read at [wikipedia] (https://en.wikipedia.org/wiki/Shapefile). An important thing to know is that whereas shape file refers to data that is used to represent geographical features, it needs 3 files to be operational.
# .shp 
where coordinates of geographical features are kept.
# .shx
an indexing file for geometric features.
# .dbf
data that is associated with the geometric features.
The folder we are going to download exists within the US Census Bureau website. From the website we navigate to where the 
[shape files for all the counties in USA are located](https://www2.census.gov/geo/tiger/TIGER2021/COUNTY/)
