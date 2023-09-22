Both R and Python have libraries that can process files referred to as "shape files". Details on the definition and development of the format can be read at [wikipedia] (https://en.wikipedia.org/wiki/Shapefile). An important thing to know is that whereas shape file refers to data that is used to represent geographical features, it needs 3 files to be operational.
### .shp 
where coordinates of geographical features are kept.
### .shx
an indexing file for geometric features.
### .dbf
data that is associated with the geometric features.
The folder we are going to download exists within the US Census Bureau website. From the website we navigate to where the 
[shape files for all the counties in USA are located](https://www2.census.gov/geo/tiger/TIGER2021/COUNTY/)

As you can see from the download, there are additional files present besides the .shp, .shx and .dbf. These are files that are not required for R and Python to create informative maps but we will keep them in our folder. 
In this project we make use of a multitude of spatial tools and functions for statistical analysis but for now we will focus only on the mapping aspects and the simple features, [sf package](https://r-spatial.github.io/sf/). 

```
<!-- Loading the library, if you are using Windows you will have to make sure you have Rtools before being able to use the library-->
library(sf)
<!--getting the shape data and reading it with the read_sf command-->
shape <- read_sf(dsn = "C:/Users/rm84/Desktop/research/HMM/data/tl_2021_us_county.shp")
<!--selecting California via FIPS state code as you can see selection is done via base R via the %in% statement-->
shape=shape[(shape$STATEFP %in% '06'),]
```
