library(rgdal)

# Preliminary in QGIS : clip the "all_roads" shapefile from the "Wallonia" project 

roads = readOGR(dsn=paste("./",sep=""), layer="Clipped")
primary = subset(roads, type=="primary")
writeOGR(primary, dsn="./", layer="Primary_roads", driver="ESRI Shapefile")
