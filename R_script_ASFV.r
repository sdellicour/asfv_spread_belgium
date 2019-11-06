library(doMC)
library(fields)
library(gdistance)
library(geometry)
library(ks)
library(lubridate)
library(maptools)
library(OpenStreetMap)
library(raster)
library(RColorBrewer)
library(rgdal)
library(rgeos)
library(sampSurf)
library(seraphim)
library(sp)
library(spatstat)

# A. Preparing the environmental data
	# A.1. Defining the projections and extents
	# A.2. Loading the rasters and shapefiles
	# A.3. Uniformising the different projections
	# A.4. Generating the overall European map
# B. Analysing the GPS collar data
	# B.1. Loading the GPS collar traces
	# B.2. Visualising the GPS traces on OpenStreetMap
	# B.3. Extracting environmental values and distances
	# B.4. Analyses of the impact of factors on movement velocity
	# B.5. Analysis of the time proportion spent in non-forest areas
	# B.6. Estimating the frequency of crossing motorway events
# C. Analysing the occurrence data
	# C.1. Loading and preparing the reported cases
	# C.2. Subsetting the data based on kernel densities
	# C.3. Estimating and plotting the wavefront velocities
		# C.3.1. Defining the analysis mask
		# C.3.2. Interpolation of first time of invasion
		# C.3.3. Estimating the friction and spread rate
		# C.3.4. Estimating the spread rate
		# C.3.5. Plotting the resulting rasters
	# C.4. Investigating the impact of environmental factors on the dispersal
		# C.4.1. Analysing the impact of factors on the wavefront dispersal velocity
		# C.4.2. Analysis of the impact of barriers on the dispersal frequency
		
showingPlots = FALSE

# A. Preparing the environmental data

	# A.1. Defining the projections and extents

wgs84 = CRS("+init=epsg:4326")
lambert72 = CRS("+init=epsg:31370")
lambert93 = CRS("+init=epsg:2154")
e1_lam93 = extent(760000, 950000, 6910000, 7080000)
e1_temp = as(e1_lam93, "SpatialPolygons")
sp::proj4string(e1_temp) = lambert93
e1_temp = sp::spTransform(e1_temp, wgs84)
e1_wgs84 = extent(e1_temp)
e2_lam93 = extent(865000, 911000, 6935700, 6975000)
e2_temp = as(e2_lam93, "SpatialPolygons")
sp::proj4string(e2_temp) = lambert93
e2_temp = sp::spTransform(e2_temp, wgs84)
e2_wgs84 = extent(e2_temp)

	# A.2. Loading the rasters and shapefiles

borders = shapefile("Raster_shapefiles/Borders_GADM_0_shapefile/Borders_GADM_0.shp")
communes = shapefile("Raster_shapefiles/Belgian_communes_shapefile/Belgian_communes.shp")
motorways = shapefile("Raster_shapefiles/Isolated_motorways_shapefile/Motorways_shapefile.shp")
primaries = shapefile("Raster_shapefiles/Isolated_primary_roads_shp/Primary_roads.shp")
fences_1 = shapefile("Raster_shapefiles/Clotures_20190731_shapefile/Clotures_20190731.shp")
barriers = shapefile("Raster_shapefiles/Clotures_20190731_merged/Barriers_20190731.shp")
clc_raster = raster("Raster_shapefiles/Corine_Land_Cover_rasters/CLC_original_raster.tif")
forest_areas = raster("Raster_shapefiles/Corine_Land_Cover_rasters/Forest_areas.asc")
agricultural_areas = raster("Raster_shapefiles/Corine_Land_Cover_rasters/Agricultural_areas.asc")
artificial_areas = raster("Raster_shapefiles/Corine_Land_Cover_rasters/Artificial_areas.asc")
land_cover = forest_areas; land_cover[agricultural_areas[]==1] = 2; land_cover[artificial_areas[]==1] = 3

	# A.3. Uniformising the different projections

crs(fences_1) = crs(communes)
crs(barriers) = crs(fences_1)
fences_2 = union(fences_1, barriers)
borders = crop(spTransform(borders, lambert93), e1_lam93)
communes = crop(spTransform(communes, lambert93), e1_lam93)
motorways = crop(spTransform(motorways, lambert93), e1_lam93)
primaries = crop(spTransform(primaries, lambert93), e1_lam93)
fences_1 = crop(spTransform(fences_1, lambert93), e1_lam93)
fences_2 = crop(spTransform(fences_2, lambert93), e1_lam93)
fences_2 = union(fences_2, motorways)

crs(forest_areas) = crs(clc_raster)
crs(agricultural_areas) = crs(clc_raster)
crs(artificial_areas) = crs(clc_raster)
crs(land_cover) = crs(clc_raster)
forest_areas = projectRaster(forest_areas, crs=lambert93)
forest_areas = crop(forest_areas, e1_lam93); r = forest_areas
forest_areas[] = (r[]-min(r[],na.rm=T))/(max(r[],na.rm=T)-min(r[],na.rm=T))
forest_areas[] = round(forest_areas[]); names(forest_areas) = "forest_areas"
agricultural_areas = projectRaster(agricultural_areas, crs=lambert93)
agricultural_areas = crop(agricultural_areas, e1_lam93); r = agricultural_areas
agricultural_areas[] = (r[]-min(r[],na.rm=T))/(max(r[],na.rm=T)-min(r[],na.rm=T))
agricultural_areas[] = round(agricultural_areas[]); names(agricultural_areas) = "agricultural_areas"
crs(forest_areas) = crs(lambert93); crs(agricultural_areas) = crs(lambert93)
land_covers_3 = projectRaster(land_cover, crs=lambert93)
land_covers_3 = crop(land_covers_3, e1_lam93); r = land_covers_3
land_covers_3[] = ((r[]-min(r[],na.rm=T))/(max(r[],na.rm=T)-min(r[],na.rm=T)))*3
land_covers_3[] = round(land_covers_3[]); crs(land_covers_3) = crs(lambert93)
study_areas = forest_areas; study_areas[!is.na(study_areas[])] = 1

# B. Analysing the GPS collar data

	# B.1. Loading the GPS collar traces

if (!file.exists("GSP_collar_2.csv"))
	{
		gps = read.csv("GSP_collar_1.csv") # removing outliers:
		gps = gps[which(gps[,"gps_validity_code"]==1),]
		times = matrix(nrow=dim(gps)[1], ncol=1)
		colnames(times) = "time"
		for (i in 1:dim(gps)[1])
			{
				time = as.character(gps[i,"acquisition_time"])
				if ((!is.na(time))&&(time != "\\N"))
					{
						time = unlist(strsplit(time," "))
						day1 = decimal_date(ymd(time[1]))
						hours1 = unlist(strsplit(gsub("\\+00","",time[2]),":"))
						day2 = ((as.numeric(hours1[1])/24)/365) + (((as.numeric(hours1[2])/60)/24)/365) + ((((as.numeric(hours1[3])/60)/60)/24)/365)
						time = day1+day2; times[i,1] = time
					}	else	{
						times[i,1] = NA
					}
			}
		gps = cbind(gps, times)
		gps = gps[!is.na(gps[,"time"]),]
		gps = gps[!is.na(gps[,"longitude"]),]
		gps = gps[!is.na(gps[,"latitude"]),]
		temp1 = gps[,c("longitude","latitude")]
		coordinates(temp1) = ~ longitude + latitude; crs(temp1) = wgs84
		temp2 = spTransform(temp1, lambert93)
		gps$x_lambert93 = temp2@coords[,1]
		gps$y_lambert93 = temp2@coords[,2]
		temp3 = spTransform(temp1, osm())
		gps$x_openStreetMap = temp3@coords[,1]
		gps$y_openStreetMap = temp3@coords[,2]
		write.csv(gps, "GSP_collar_2.csv", row.names=F, quote=F)
	}

	# B.2. Visualising the GPS traces on OpenStreetMap

gps = read.csv("GSP_collar_2.csv"); individuals = list(); individual_names = list()
individual_ids = unique(gps[,"animals_original_id"])
for (i in 1:length(individual_ids))
	{
		lines = which(gps[,"animals_original_id"]==individual_ids[i])
		individual = gps[lines,c("time","longitude","latitude","x_lambert93","y_lambert93","x_openStreetMap","y_openStreetMap")]
		individual = individual[order(individual[,"time"]),]
		individual_name = paste(gps[lines[1],"short_name"],gps[lines[1],"animals_original_id"],sep="_")
		if (individual_ids[i] == "su398bleu2005")
			{
				index = which(individual[,"time"]==2006.41495494039)
				individual = individual[-index,] # to manually remove a remaining outlier
			}
		write.csv(individual, paste0("GSP_collar_data/",individual_name,".csv"), row.names=F, quote=F)
		individuals[[i]] = individual; individual_names[[i]] = individual_name
	}
upperLeft = c(e1_wgs84@ymax,e1_wgs84@xmin); lowerRight=c(e1_wgs84@ymin,e1_wgs84@xmax)
openStreetMap = openmap(upperLeft, lowerRight, type="osm", zoom=9); OSMap = TRUE
for (i in 1:length(individuals))
	{
		if (!file.exists(paste0("GSP_collar_data/",individual_names[[i]],".pdf")))
			{
				if (OSMap == TRUE)
					{
						pdf(paste0("GSP_collar_data/",individual_names[[i]],".pdf"), width=6.0, height=5.5)
						r = raster(openStreetMap); # dev.new(width=6.0, height=5.5)
						par(mar=c(0,0,0,0), oma=c(0,2,0,2), mgp=c(0,0.4,0), lwd=0.2, bty="o")
						plot(openStreetMap)
						for (j in 2:dim(individuals[[i]])[1])
							{
								pts = cbind(individuals[[i]][j-1,c("x_openStreetMap")],individuals[[i]][j-1,c("y_openStreetMap")])
								pts = rbind(pts,cbind(individuals[[i]][j,c("x_openStreetMap")],individuals[[i]][j,c("y_openStreetMap")]))
								lines(pts, lwd=0.3, col="gray10")
							}
					}	else	{
						pdf(paste0(individual_names[[i]],".pdf"), width=6.0, height=5.5)
						r = forest_areas; # dev.new(width=6.0, height=5.5)
						par(mar=c(0,4,0,0), oma=c(0,0,0,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
						plot(forest_areas, col=c(rgb(0,1,0,0),rgb(68/255,165/255,68/255,0.3)), box=F, axes=F, legend=F)
						plot(belgium, add=T, lwd=1.0, border="white", lty=1)
						plot(belgium, add=T, lwd=0.25, border="gray30", lty=1)
						plot(motorways, add=T, lwd=0.75, col="red", lty=1)
						for (j in 2:dim(individuals[[i]])[1])
							{
								pts = cbind(individuals[[i]][j-1,c("x_lambert93")],individuals[[i]][j-1,c("y_lambert93")])
								pts = rbind(pts,cbind(individuals[[i]][j,c("x_lambert93")],individuals[[i]][j,c("y_lambert93")]))
								lines(pts, lwd=0.3, col="gray10")
							}
					}
				rect(xmin(r), ymin(r), xmax(r), ymax(r), xpd=T, lwd=0.2, border="gray30")
				dev.off()
			}
	}
if (!file.exists(paste0("GSP_data_2.pdf")))
	{
		if (OSMap == FALSE) data1 = read.csv("Data_130919.csv", header=T)
		l = length(individuals); p = round(length(individuals)/3)+1; OSMap = FALSE
		cols = colorRampPalette(brewer.pal(11,"Spectral"))(l+p)[c(1:(l/2),((l/2)+p):(l+p))]
		if (OSMap == TRUE)
			{
				pdf(paste0("GSP_data_2.pdf"), width=6.0, height=5.5)
				r = raster(openStreetMap); # dev.new(width=6.0, height=5.5)
				par(mar=c(0,0,0,0), oma=c(0,2,0,2), mgp=c(0,0.4,0), lwd=0.2, bty="o")
				plot(openStreetMap)
				for (i in 1:length(individuals))
					{
						for (j in 2:dim(individuals[[i]])[1])
							{
								pts = cbind(individuals[[i]][j-1,c("x_openStreetMap")],individuals[[i]][j-1,c("y_openStreetMap")])
								pts = rbind(pts,cbind(individuals[[i]][j,c("x_openStreetMap")],individuals[[i]][j,c("y_openStreetMap")]))
								lines(pts, lwd=0.3, col="gray10")
							}
					}
				rect(xmin(r), ymin(r), xmax(r), ymax(r), xpd=T, lwd=0.2, border="gray30")
			}	else 	{
				pdf(paste0(individual_names[[i]],".pdf"), width=6.0, height=5.5)
				r = forest_areas; # dev.new(width=6.0, height=5.5)
				par(mar=c(0,4,0,0), oma=c(0,0,0,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
				plot(land_covers_3, col=c(rgb(0,1,0,0),rgb(68/255,165/255,68/255,0.2),rgb(0,1,0,0),"gray80"), box=F, axes=F, legend=F)
				plot(borders, add=T, lwd=3, col="white", lty=1)
				plot(borders, add=T, lwd=0.5, col="blue")
				plot(motorways, add=T, lwd=0.75, col="red", lty=1)
				for (i in 1:length(individuals))
					{
						GPS_segments = cbind(individuals[[i]][2:dim(individuals[[i]])[1],4:5],individuals[[i]][1:dim(individuals[[i]])[1]-1,4:5])
						segments(GPS_segments[,1],GPS_segments[,2],GPS_segments[,3],GPS_segments[,4], lwd=0.3, col="gray10")
					}
				points(data1[,c("x_transformed","y_transformed")], cex=0.5,  lwd=0.5, pch=3, col="green3")
				segments(e2_lam93@xmin, e2_lam93@ymin, e2_lam93@xmax, e2_lam93@ymin, col="gray30", lwd=0.75, lty=1)
				segments(e2_lam93@xmax, e2_lam93@ymin, e2_lam93@xmax, e2_lam93@ymax, col="gray30", lwd=0.75, lty=1)
				segments(e2_lam93@xmax, e2_lam93@ymax, e2_lam93@xmin, e2_lam93@ymax, col="gray30", lwd=0.75, lty=1)
				segments(e2_lam93@xmin, e2_lam93@ymax, e2_lam93@xmin, e2_lam93@ymin, col="gray30", lwd=0.75, lty=1)
				rect(xmin(r), ymin(r), xmax(r), ymax(r), xpd=T, lwd=0.2, border="gray30")
			}
		dev.off()
	}

	# B.3. Extracting environmental values and distances
	
environmentalExtraction1 = function(traces, rasters, resistances)
	{
		extractions = matrix(nrow=dim(traces)[1], ncol=length(traces))
		for (i in 2:dim(traces)[1])
			{
				point1 = traces[i-1,c("x_lambert93","y_lambert93")]; names(point1) = c("x","y")
				point2 = traces[i,c("x_lambert93","y_lambert93")]; names(point2) = c("x","y")
				points = rbind(point1,point2); linesList = list()
				linesList[[1]] = Lines(list(Line(points)),1)
				spatialLine = SpatialLines(linesList)
				for (j in 1:length(rasters))
					{
						values = raster::extract(rasters[[j]], spatialLine)[[1]]
						if (resistances[j] == FALSE) values = 1/values
						extractions[i,j] = sum(values, na.rm=T)
					}
			}
		colNames = c()
		for (i in 1:length(rasters)) colNames = c(colNames, names(rasters[[i]]))
		colnames(extractions) = colNames
		return(extractions)
	}
environmentalExtraction2 = function(traces, rasters)
	{
		extractions = matrix(nrow=dim(traces)[1], ncol=length(rasters))
		for (i in 1:dim(traces)[1])
			{
				point = traces[i,c("x_lambert93","y_lambert93")]; names(point) = c("x","y")
				for (j in 1:length(rasters))
					{
						extractions[i,j] = raster::extract(rasters[[j]], point)[[1]]
					}
			}
		colNames = c()
		for (i in 1:length(rasters)) colNames = c(colNames, names(rasters[[i]]))
		colnames(extractions) = colNames
		return(extractions)
	}

rasters1 = list(); resistances = c(); c = 0
null_raster = forest_areas
names(null_raster) = "null_raster"
null_raster[!is.na(null_raster[])] = 1
c = c + 1; rasters1[[c]] = null_raster
resistances = c(resistances, TRUE)
kS = c(10, 100, 1000)
for (k in kS)
	{
		c = c + 1
		r = forest_areas
		r[!is.na(r[])] = 1+(k*(r[!is.na(r[])]))
		names(r) = paste0("forest_areas_k",k)
		rasters1[[c]] = r
		resistances = c(resistances, FALSE)
	}
for (k in kS)
	{
		c = c + 1
		r = rasterize(motorways, null_raster, fun=min)
		r[!is.na(r[])] = 1; r[(is.na(r[]))&(!is.na(null_raster[]))] = 0
		r[!is.na(r[])] = 1+(k*(r[!is.na(r[])]))
		names(r) = paste0("motorways_k",k)
		rasters1[[c]] = r
		resistances = c(resistances, TRUE)
	}
rasters2 = list()
null_raster = forest_areas
null_raster[!is.na(null_raster[])] = 1
rasters2[[1]] = forest_areas

registerDoMC(10); buffer = list()
buffer = foreach(i = 1:length(individuals)) %dopar% {
# for (i in 1:length(individuals)) {
		extractions_obs_1 = environmentalExtraction1(individuals[[i]], rasters1, resistances)
		write.csv(extractions_obs_1, paste0("GPS_extractions/Individual_",i,"_extractions_1.csv"), row.names=F, quote=F)
		extractions_obs_2 = environmentalExtraction2(individuals[[i]], rasters2)
		write.csv(extractions_obs_2, paste0("GPS_extractions/Individual_",i,"_extractions_2.csv"), row.names=F, quote=F)
		i
	}

	# B.4. Analyses of the impact of factors on movement velocity

nberOfRandomisations = 100
colNames = c(); rowNames = c()
for (i in 2:length(rasters1))
	{
		colNames = c(colNames, names(rasters1[[i]]))
	}
for (i in 1:length(individuals))
	{
		rowNames = c(rowNames, individual_names[[i]])
	}
Qs = matrix(nrow=length(individuals), ncol=6); colnames(Qs) = colNames; row.names(Qs) = rowNames
pV = matrix(nrow=length(individuals), ncol=6); colnames(pV) = colNames; row.names(Qs) = rowNames
for (i in 1:length(individuals))
	{
		dt = individuals[[i]][2:dim(individuals[[i]])[1],"time"]-individuals[[i]][1:dim(individuals[[i]])[1]-1,"time"]
		extractions_obs_1 = read.csv(paste0("GPS_extractions/Individual_",i,"_extractions_1.csv"), header=T)
		x = extractions_obs_1[2:dim(extractions_obs_1)[1],"null_raster"]; y = dt
		LR_obs = lm(as.formula(paste0("y ~ x"))); R2_nul_obs = summary(LR_obs)$r.squared
		for (j in 2:length(rasters1))
			{
				x = extractions_obs_1[2:dim(extractions_obs_1)[1],j]; c = 0
				LR_obs = lm(as.formula(paste0("y ~ x"))); R2_env_obs = summary(LR_obs)$r.squared
				Q_obs = R2_env_obs - R2_nul_obs
				for (k in 1:nberOfRandomisations)
					{
						x = extractions_obs_1[2:dim(extractions_obs_1)[1],"null_raster"]; y = sample(dt, length(dt), replace=F)
						LR_ran = lm(as.formula(paste0("y ~ x"))); R2_nul_ran = summary(LR_ran)$r.squared
						x = extractions_obs_1[2:dim(extractions_obs_1)[1],j]
						LR_ran = lm(as.formula(paste0("y ~ x"))); R2_env_ran = summary(LR_ran)$r.squared
						Q_ran = R2_env_ran - R2_nul_ran
						if (Q_obs < Q_ran) c = c+1
					}
				Qs[i,j-1] = Q_obs; pV[i,j-1] = c/nberOfRandomisations
			}
	}
write.csv(Qs, "GPS_Q_statistics.csv", quote=F); write.csv(pV, "GPS_Q_pValues.csv", quote=F)

	# B.5. Analysis of the time proportion spent in non-forest areas

n = 0; n_tot = 0; Fs = matrix(nrow=length(individuals), ncol=1)
colnames(Fs) = "%_in_forest_cells"; row.names(Fs) = rowNames
for (i in 1:length(individuals))
	{
		extractions_obs_2 = read.csv(paste0("GPS_extractions/Individual_",i,"_extractions_2.csv"), header=T)
		if (dim(extractions_obs_2)[1] >= 100)
			{
				Fs[i,1] = (sum(extractions_obs_2[,1]==1)/dim(extractions_obs_2)[1])*100
				n = n+sum(extractions_obs_2[,1]==1); n_tot = n_tot+dim(extractions_obs_2)[1]
			}
	}
# print(mean(Fs, na.rm=T)) # 74.448141 %
# print((n/n_tot)*100) # 69.44071 %
if (showingPlots == TRUE)
	{
		dev.new(width=3.5, height=2); par(mgp=c(1,0.35,0), oma=c(0.5,1.3,1,2), mar=c(2.5,2,0.5,0), lwd=0.2)
		colsa = rgb(120,120,120,255,maxColorValue=255); colsb = rgb(120,120,120,100,maxColorValue=255)
		hist(100-Fs[,"%_in_forest_cells"], border=colsa, col=colsb, breaks=50, main=NA, axes=F, ann=F, xlim=c(0,100))
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,-0.03,0), lwd=0.2, tck=-0.025, col="gray30", col.tick="gray30", col.axis="gray30", at=seq(0,100,10))
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.20,0), lwd=0.2, tck=-0.025, col="gray30", col.tick="gray30", col.axis="gray30", at=seq(0,20,5))
		title(xlab="% of GPS records in non-forest areas", cex.lab=0.7, mgp=c(1.0,0,0), col.lab="gray30")
		title(ylab="frequency", cex.lab=0.7, mgp=c(1.2,0,0), col.lab="gray30")
	}

	# B.6. Estimating the frequency of crossing motorway events

Ns = matrix(nrow=length(individuals), ncol=1); colnames(Ns) = "#_crossing_motorway_events"; row.names(Ns) = rowNames
for (i in 1:length(individuals))
	{
		ind = individuals[[i]]; N = 0
		for (j in 2:dim(ind)[1])
			{
				point1 = ind[j-1,c("x_lambert93","y_lambert93")]; names(point1) = c("x","y")
				point2 = ind[j,c("x_lambert93","y_lambert93")]; names(point2) = c("x","y")
				points = rbind(point1,point2); linesList = list()
				linesList[[1]] = Lines(list(Line(points)),1)
				spatialLine = SpatialLines(linesList); crs(spatialLine) = crs(motorways)
				intersections = gIntersection(motorways, spatialLine)
				if (!is.null(intersections)) N = N + 1
			}
		Ns[i,1] = N; print(c(rowNames[i],Ns[i,1]))
	}

# C. Analysing the occurrence data

	# C.1. Loading and preparing the reported cases

data1 = read.csv("Data_130919.csv", header=T); maxDays = max(data1[,"days"])
temp1 = data1; coordinates(temp1) = ~ longitude + latitude; crs(temp1) = wgs84
myOutProj = CRS("+proj=lcc +lat_1=18 +lat_2=24 +lat_0=21 +lon_0=114 +x_0=500000 +y_0=500000 
				 +ellps=WGS72 +towgs84=0,0,1.9,0,0,0.814,-0.38 +units=m +no_defs") # not used
temp2 = spTransform(temp1, lambert93)
data1$x_transformed = temp2@coords[,1]
data1$y_transformed = temp2@coords[,2]

	# C.2. Subsetting the data based on kernel densties

if (!file.exists("Data_KDE_95.csv"))
	{
		template = crop(forest_areas, e2_lam93)
		rast = template; rast[] = NA
		gridSize = c(rast@ncols, rast@nrows)
		xyMin = c(rast@extent@xmin, rast@extent@ymin)
		xyMax = c(rast@extent@xmax, rast@extent@ymax)
		H = Hpi(data1[,c("x_transformed","y_transformed")])
		if (showingPlots == TRUE)
			{
				kde = kde(data1[,c("x_transformed","y_transformed")], H=H, compute.cont=T, gridsize=gridSize, xmin=xyMin, xmax=xyMax)
				plot(data1[,c("x_transformed","y_transformed")], axes=F, ann=F)
				r = raster(kde); c = rasterToContour(r, levels=kde$cont["30%"]); lines(c, col="gray30")
				r = raster(kde); c = rasterToContour(r, levels=kde$cont["10%"]); lines(c, col="red")
				r = raster(kde); c = rasterToContour(r, levels=kde$cont["5%"]); lines(c, col="orange")
				r = raster(kde); c = rasterToContour(r, levels=kde$cont["1%"]); lines(c, col="green3")
			}
		buffer = data1[which(data1[,"days"]==0),]
		data2 = data1[which(data1[,"days"]==0),]
		for (i in 1:maxDays)
			{
				if (showingPlots == FALSE) cat(paste0("day ",i,"\n"))
				indices1 = which(data1[,"days"]==i)
				if (length(indices1) > 0)
					{
						if (dim(data2)[1] >= 3)
							{
								indices2 = c()
								# H = Hpi(buffer[,c("x_transformed","y_transformed")]) # OLD version with one H per time slice
								kde = kde(buffer[,c("x_transformed","y_transformed")], H=H, compute.cont=T, gridsize=gridSize, xmin=xyMin, xmax=xyMax)
								r = raster(kde); contour = rasterToContour(r, levels=kde$cont["5%"]); col_pt = NA
								threshold = contourLevels(kde, 0.05); r[r[]<threshold] = NA; r[!is.na(r[])] = i; crs(r) = crs(rast)
								writeRaster(r, paste0("./KDE_95_contours/KDE_95_contours_day_",i-1,".tif"), overwrite=T)
								writeOGR(contour, dsn="./KDE_95_contours", layer=paste0("KDE_95_contours_day_",i-1), driver="ESRI Shapefile")
								if (showingPlots == TRUE) plot(contour, main=paste0("day ",i), cex.main=0.9, col.main="gray30")
								for (j in 1:length(indices1))
									{
										point_in_polygon = FALSE
										pt.x = data1[indices1[j],"x_transformed"]
										pt.y = data1[indices1[j],"y_transformed"]
										for (k in 1:length(contour@lines[[1]]@Lines))
											{
												pol.x = contour@lines[[1]]@Lines[[k]]@coords[,1]
												pol.y = contour@lines[[1]]@Lines[[k]]@coords[,2]
												if (point.in.polygon(pt.x, pt.y, pol.x, pol.y) != 0)
													{
														point_in_polygon = TRUE
													}
											}
										if (point_in_polygon == FALSE)
											{
												indices2 = c(indices2, indices1[j])
												col_pt = "green3"
											}	else		{
												col_pt = "red"
											}
										if (showingPlots == TRUE)
											{
												points(data1[indices1[j],c("x_transformed","y_transformed")], col=col_pt)
											}
									}
								if (length(indices2) > 0)
									{
										data2 = rbind(data2, data1[indices2,])
									}
							}	else	{
								data2 = rbind(data2, data1[indices1,])
							}
						buffer = rbind(buffer, data1[indices1,])
					}
			}
		data2 = unique(data2)
		index = which(data2[,"days"]==min(data2[,"days"]))
		buffer = data2[index,]
		lines = paste(data2[index,"x_transformed"],data2[index,"y_transformed"],sep="_")
		for (d in 1:max(data2[,"days"]))
			{
				indices1 = which(data2[,"days"]==d)
				if (length(indices1) > 0)
					{
						for (i in 1:length(indices1))
							{
								line = paste(data2[indices1[i],"x_transformed"],data2[indices1[i],"y_transformed"],sep="_")
								if (sum(lines==line) == 0)
									{
										lines = c(lines, line); buffer = rbind(buffer, data2[indices1[i],])
									}
							}
					}
			}
		data2 = buffer
		write.csv(data2, "Data_KDE_95.csv", quote=F, row.names=F)
		for (i in 1:maxDays)
			{
				if (file.exists(paste0("./KDE_95_contours/KDE_95_contours_day_",i,".tif")))
					{
						kde = raster(paste0("./KDE_95_contours/KDE_95_contours_day_",i,".tif"))
						rast = merge(rast, kde)
					}
			}
		writeRaster(rast, "Data_KDE_95.tif")
	}

	# C.3. Estimating and plotting the wavefront velocities

		# C.3.1. Defining the analysis mask

data2 = read.csv("Data_KDE_95.csv", header=T)
template = forest_areas; template[!is.na(template[])] = 0
H = Hpi(data2[,c("x_transformed","y_transformed")])
kde = kde(data2[,c("x_transformed","y_transformed")], H=H, compute.cont=T, gridsize=c(1000,1000))
rast1 = raster(kde); contour = rasterToContour(rast1, levels=kde$cont["5%"])
threshold = min(raster::extract(rast1, data2[,c("x_transformed","y_transformed")]))
rast2 = rast1; rast2[rast1[]<=threshold] = NA # previous version to include all points
p = Polygon(contour@lines[[1]]@Lines[[1]]@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
rast2 = mask(rast1, sps); mask = raster::resample(rast2, template)

		# C.3.2. Interpolation of first time of invasion

tps_model = Tps(x=data2[,c("x_transformed","y_transformed")], Y=data2[,"days"])
tps = interpolate(template, tps_model)
tps_mask = crop(mask(tps, mask), e2_lam93)

		# C.3.3. Estimating the friction and spread rate

			# C.3.3.1. Measuring the local slope/friction using a 3x3 moving windows filter

raster_resolution = mean(c(res(tps_mask)[1],res(tps_mask)[2]))
f = matrix(1/raster_resolution, nrow=3, ncol=3)
f[c(1,3,7,9)] = 1/(sqrt(2)*raster_resolution); f[5] = 0
fun = function(x, ...)
	{ 
		sum(abs(x-x[5])*f)/8
	}
friction = focal(tps_mask, w=matrix(1,nrow=3,ncol=3), fun=fun, pad=T, padValue=NA, na.rm=F)

			# C.3.3.2. Smoothing the resulting friction surface using an average 11x11 cell filter

myAvSize = 11
friction_sd10 = focal(friction, w=matrix(1/(myAvSize^2),nrow=myAvSize,ncol=myAvSize), pad=T, padValue=NA, na.rm=F)

		# C.3.4. Estimating the spread rate

spreadRate_sd10 = ((1/(friction_sd10))/1000)*7
print(mean(spreadRate_sd10[], na.rm=T)) # 0.39 km/week

		# C.3.5. Plotting the resulting rasters

if (showingPlots == TRUE)
	{
		workflow_figure = FALSE
		forest_areas_cropped = crop(forest_areas, e2_lam93)
		land_covers_3_cropped = crop(land_covers_3, e2_lam93)
		borders_cropped = crop(borders, e2_lam93)
		fences_1_cropped = crop(fences_1, e2_lam93)
		motorways_cropped = crop(motorways, e2_lam93)
		primaries_cropped = crop(primaries, e2_lam93)
		successiveKDEs = raster("Data_KDE_95.tif")
		successiveKDEs[1] = max(data1[,"days"]); successiveKDEs[2] = 0
		tps_mask_cropped = mask(tps_mask, belgium_cropped)
		tps_mask_cropped[tps_mask_cropped[]<0] = 0
		spreadRate_sd10_cropped = mask(spreadRate_sd10, belgium_cropped)
		r = spreadRate_sd10_cropped; r[(!is.na(r[]))&(r[]>1)] = 1
		spreadRate_sd10_truncated = r

		dev.new(width=9.3, height=4)
		par(mfrow=c(1,2), mar=c(0.5,0.5,0.5,2.0), oma=c(1.0,1.7,1.0,1.3), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		cols1 = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(161)[21:121])
		cols2 = rev(colorRampPalette(brewer.pal(11,"RdYlGn"))(181)[21:121])
		colsP = cols1[(((data1[,"days"]-min(data1[,"days"]))/(max(data1[,"days"])-min(data1[,"days"])))*100)+1]
		plot(land_covers_3_cropped, col=c(rgb(0,1,0,0),rgb(68/255,165/255,68/255,0.2),rgb(0,1,0,0),"gray80"), box=F, axes=F, legend=F)
		for (i in dim(data1)[1]:1)
			{
				points(data1[i,c("x_transformed","y_transformed")], pch=16, cex=0.7, col=colsP[i])
				points(data1[i,c("x_transformed","y_transformed")], pch=1, cex=0.7, col="gray30", lwd=0.2)
			}
		mtext("1. First invasion times", side=3, line=-2, at=901000, cex=0.8, font=1, col="gray30")
		mtext("(all, in days)", side=3, line=-2.8, at=904500, cex=0.8, font=1, col="gray30")
		legend(891000, 6940800, c("Fences","Motorways"), lwd=c(0.75,0.75), cex=0.75, col=c("gray30","red"), text.col=c("gray30","red"), border=NA, x.intersp=0.7, y.intersp=0.8, bty="n")
		rect(xmin(tps_mask), ymin(tps_mask), xmax(tps_mask), ymax(tps_mask), xpd=T, lwd=0.2, border="gray30")
		legendRast = raster(as.matrix(seq(0,max(data1[,"days"]),1)))
		plot(legendRast, legend.only=T, col=cols1, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.895,0.910,0.06,0.941),
		     alpha=1, legend.args=list(text="", cex=0.5, line=0.5, col="gray30"), axis.args=list(cex.axis=0.6, lwd=0,
		     lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.45,0), at=seq(0,360,30)))
		colsR = cols1[((((1:max(successiveKDEs[],na.rm=T))-min(data1[,"days"]))/(max(data1[,"days"])-min(data1[,"days"])))*100)+1]
		colsP = cols1[(((data2[,"days"]-min(data1[,"days"]))/(max(data1[,"days"])-min(data1[,"days"])))*100)+1]
		plot(successiveKDEs, main="", cex.main=1, cex.axis=0.7, bty="n", box=F, axes=F, legend=F, axis.args=list(cex.axis=0.7), col=colsR, colNA="white")
		plot(land_covers_3_cropped, add=T, col=c(rgb(0,1,0,0),rgb(68/255,165/255,68/255,0.2),rgb(0,1,0,0),"gray80"), box=F, axes=F, legend=F)
		for (i in dim(data2)[1]:1)
			{
				points(data2[i,c("x_transformed","y_transformed")], pch=16, cex=0.7, col=colsP[i])
				points(data2[i,c("x_transformed","y_transformed")], pch=1, cex=0.7, col="gray30", lwd=0.2)
			}
		mtext("2. First invasion times", side=3, line=-2, at=901000, cex=0.8, font=1, col="gray30")
		mtext("(filtered, in days)", side=3, line=-2.8, at=903000, cex=0.8, font=1, col="gray30")
		rect(xmin(spreadRate_sd10), ymin(spreadRate_sd10), xmax(spreadRate_sd10), ymax(spreadRate_sd10), xpd=T, lwd=0.2, border="gray30")
		plot(legendRast, legend.only=T, col=cols1, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.895,0.910,0.06,0.941),
		     alpha=1, legend.args=list(text="", cex=0.5, line=0.5, col="gray30"), axis.args=list(cex.axis=0.6, lwd=0,
		     lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.45,0), at=seq(0,360,30)))
		dev.new(width=9.3, height=4)
		par(mfrow=c(1,2), mar=c(0.5,0.5,0.5,2.0), oma=c(1.0,1.7,1.0,1.3), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		cols1 = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(161)[21:121])
		cols2 = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(161)[41:121])
		colsR = cols1[((((1:max(tps_mask_cropped[],na.rm=T))-min(data1[,"days"]))/(max(data1[,"days"])-min(data1[,"days"])))*100)+1]
		if (workflow_figure != FALSE) colsP = cols1[(((data2[,"days"]-min(data1[,"days"]))/(max(data1[,"days"])-min(data1[,"days"])))*100)+1]
		plot(tps_mask_cropped, main="", cex.main=1, cex.axis=0.7, bty="n", box=F, axes=F, legend=F, axis.args=list(cex.axis=0.7), col=cols1, colNA="white")
		plot(land_covers_3_cropped, add=T, col=c(rgb(0,1,0,0),rgb(68/255,165/255,68/255,0.2),rgb(0,1,0,0),"gray80"), box=F, axes=F, legend=F)
		# points(data2[,c("x_transformed","y_transformed")], pch=3, cex=0.5, lwd=0.4, col="gray30")
		for (i in dim(data2)[1]:1)
			{
				points(data2[i,c("x_transformed","y_transformed")], pch=16, cex=0.7, col=colsP[i])
				points(data2[i,c("x_transformed","y_transformed")], pch=1, cex=0.7, col="gray30", lwd=0.2)
			}
		mtext("3. First invasion times", side=3, line=-2, at=901000, cex=0.8, font=1, col="gray30")
		mtext("(interpolated, in days)", side=3, line=-2.8, at=901000, cex=0.8, font=1, col="gray30")
		rect(xmin(tps_mask), ymin(tps_mask), xmax(tps_mask), ymax(tps_mask), xpd=T, lwd=0.2, border="gray30")
		plot(legendRast, legend.only=T, col=cols1, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.895,0.910,0.06,0.941),
		     alpha=1, legend.args=list(text="", cex=0.5, line=0.5, col="gray30"), axis.args=list(cex.axis=0.6, lwd=0,
		     lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.45,0), at=seq(0,360,30)))
		plot(spreadRate_sd10_truncated, main="", cex.main=1, cex.axis=0.7, bty="n", box=F, axes=F, legend=F, axis.args=list(cex.axis=0.7), col=cols2, colNA="white")
		plot(land_covers_3_cropped, add=T, col=c(rgb(0,1,0,0),rgb(68/255,165/255,68/255,0.2),rgb(0,1,0,0),"gray80"), box=F, axes=F, legend=F)
		points(data2[,c("x_transformed","y_transformed")], pch=3, cex=0.5, lwd=0.4, col="gray30")
		mtext("4. Wavefront velocity", side=3, line=-2, at=901500, cex=0.8, font=1, col="gray30")
		mtext("(km/week)", side=3, line=-2.8, at=905300, cex=0.8, font=1, col="gray30")
		rect(xmin(spreadRate_sd10), ymin(spreadRate_sd10), xmax(spreadRate_sd10), ymax(spreadRate_sd10), xpd=T, lwd=0.2, border="gray30")
		plot(spreadRate_sd10_truncated, legend.only=T, add=T, col=cols2, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.895,0.910,0.06,0.941), 
		     alpha=1, legend.args=list(text="", cex=0.5, line=0.5, col="gray30"), axis.args=list(cex.axis=0.6, lwd=0,
		     lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.45,0), at=seq(0,1,0.1), labels=c(seq(0,0.9,0.1),">=1")))
		     
		dev.new(width=9.3, height=4)
		par(mfrow=c(1,2), mar=c(0.5,0.5,0.5,2.0), oma=c(1.0,1.7,1.0,1.3), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		cols1 = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(161)[21:121])
		cols2 = rev(colorRampPalette(brewer.pal(11,"RdYlGn"))(181)[21:121])
		colsP = cols1[(((data1[,"days"]-min(data1[,"days"]))/(max(data1[,"days"])-min(data1[,"days"])))*100)+1]
		plot(land_covers_3_cropped, col=c(rgb(0,1,0,0),rgb(68/255,165/255,68/255,0.2),rgb(0,1,0,0),"gray80"), box=F, axes=F, legend=F)
		plot(primaries_cropped, add=T, lwd=0.5, col="gray80")
		plot(borders_cropped, add=T, lwd=3, col="white", lty=1)
		plot(borders_cropped, add=T, lwd=0.5, col="blue")
		plot(fences_1_cropped, add=T, lwd=0.75, col="gray30", lty=1)
		plot(motorways_cropped, add=T, lwd=0.75, col="red", lty=1)
		for (i in dim(data1)[1]:1)
			{
				points(data1[i,c("x_transformed","y_transformed")], pch=16, cex=0.7, col=colsP[i])
				points(data1[i,c("x_transformed","y_transformed")], pch=1, cex=0.7, col="gray30", lwd=0.2)
			}
		mtext("1. First invasion times", side=3, line=-2, at=901000, cex=0.8, font=1, col="gray30")
		mtext("(all, in days)", side=3, line=-2.8, at=904500, cex=0.8, font=1, col="gray30")
		legend(891000, 6940800, c("Fences","Motorways"), lwd=c(0.75,0.75), cex=0.75, col=c("gray30","red"), text.col=c("gray30","red"), border=NA, x.intersp=0.7, y.intersp=0.8, bty="n")
		rect(xmin(tps_mask), ymin(tps_mask), xmax(tps_mask), ymax(tps_mask), xpd=T, lwd=0.2, border="gray30")
		legendRast = raster(as.matrix(seq(0,max(data1[,"days"]),1)))
		plot(legendRast, legend.only=T, col=cols1, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.895,0.910,0.06,0.941),
		     alpha=1, legend.args=list(text="", cex=0.5, line=0.5, col="gray30"), axis.args=list(cex.axis=0.6, lwd=0,
		     lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.45,0), at=seq(0,360,30)))
		colsP = cols1[(((data2[,"days"]-min(data2[,"days"]))/(max(data2[,"days"])-min(data2[,"days"])))*100)+1]
		plot(land_covers_3_cropped, col=c(rgb(0,1,0,0),rgb(68/255,165/255,68/255,0.2),rgb(0,1,0,0),"gray80"), box=F, axes=F, legend=F)
		plot(primaries_cropped, add=T, lwd=0.5, col="gray80")
		plot(borders_cropped, add=T, lwd=3, col="white", lty=1)
		plot(borders_cropped, add=T, lwd=0.5, col="blue")
		plot(fences_1_cropped, add=T, lwd=0.75, col="gray30", lty=1)
		plot(motorways_cropped, add=T, lwd=0.75, col="red", lty=1)
		for (i in dim(data2)[1]:1)
			{
				points(data2[i,c("x_transformed","y_transformed")], pch=16, cex=0.7, col=colsP[i])
				points(data2[i,c("x_transformed","y_transformed")], pch=1, cex=0.7, col="gray30", lwd=0.2)
			}
		mtext("2. First invasion times", side=3, line=-2, at=901000, cex=0.8, font=1, col="gray30")
		mtext("(filtered, in days)", side=3, line=-2.8, at=903000, cex=0.8, font=1, col="gray30")
		rect(xmin(spreadRate_sd10), ymin(spreadRate_sd10), xmax(spreadRate_sd10), ymax(spreadRate_sd10), xpd=T, lwd=0.2, border="gray30")
		legendRast = raster(as.matrix(seq(0,max(data2[,"days"]),1)))
		plot(legendRast, legend.only=T, col=cols1, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.895,0.910,0.06,0.941),
		     alpha=1, legend.args=list(text="", cex=0.5, line=0.5, col="gray30"), axis.args=list(cex.axis=0.6, lwd=0,
		     lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.45,0), at=seq(0,360,30)))
		dev.new(width=9.3, height=4)
		par(mfrow=c(1,2), mar=c(0.5,0.5,0.5,2.0), oma=c(1.0,1.7,1.0,1.3), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		cols1 = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(161)[21:121])
		cols2 = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(161)[41:121])
		if (workflow_figure != FALSE) colsP = cols1[(((data2[,"days"]-min(data1[,"days"]))/(max(data1[,"days"])-min(data1[,"days"])))*100)+1]
		plot(tps_mask_cropped, main="", cex.main=1, cex.axis=0.7, bty="n", box=F, axes=F, legend=F, axis.args=list(cex.axis=0.7), col=cols1, colNA="white")
		plot(land_covers_3_cropped, add=T, col=c(rgb(0,1,0,0),rgb(68/255,165/255,68/255,0.2),rgb(0,1,0,0),"gray80"), box=F, axes=F, legend=F)
		plot(primaries_cropped, add=T, lwd=0.5, col="gray80")
		points(data2[,c("x_transformed","y_transformed")], pch=3, cex=0.5, lwd=0.4, col="gray30")
		plot(borders_cropped, add=T, lwd=3, col="white", lty=1)
		plot(borders_cropped, add=T, lwd=0.5, col="blue")
		plot(fences_1_cropped, add=T, lwd=0.75, col="gray30", lty=1)
		plot(motorways_cropped, add=T, lwd=0.75, col="red", lty=1)
		mtext("3. First invasion times", side=3, line=-2, at=901000, cex=0.8, font=1, col="gray30")
		mtext("(interpolated, in days)", side=3, line=-2.8, at=901000, cex=0.8, font=1, col="gray30")
		rect(xmin(tps_mask), ymin(tps_mask), xmax(tps_mask), ymax(tps_mask), xpd=T, lwd=0.2, border="gray30")
		plot(tps_mask_cropped, legend.only=T, add=T, col=cols1, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.895,0.910,0.06,0.941),
		     alpha=1, legend.args=list(text="", cex=0.5, line=0.5, col="gray30"), axis.args=list(cex.axis=0.6, lwd=0,
		     lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.45,0), at=seq(0,360,30)))
		plot(spreadRate_sd10_truncated, main="", cex.main=1, cex.axis=0.7, bty="n", box=F, axes=F, legend=F, axis.args=list(cex.axis=0.7), col=cols2, colNA="white")
		plot(land_covers_3_cropped, add=T, col=c(rgb(0,1,0,0),rgb(68/255,165/255,68/255,0.2),rgb(0,1,0,0),"gray80"), box=F, axes=F, legend=F)
		plot(primaries_cropped, add=T, lwd=0.5, col="gray80")
		points(data2[,c("x_transformed","y_transformed")], pch=3, cex=0.5, lwd=0.4, col="gray30")
		plot(borders_cropped, add=T, lwd=3, col="white", lty=1)
		plot(borders_cropped, add=T, lwd=0.5, col="blue")
		plot(fences_1_cropped, add=T, lwd=0.75, col="gray30", lty=1)
		plot(motorways_cropped, add=T, lwd=0.75, col="red", lty=1)
		mtext("4. Wavefront velocity", side=3, line=-2, at=901500, cex=0.8, font=1, col="gray30")
		mtext("(km/week)", side=3, line=-2.8, at=905300, cex=0.8, font=1, col="gray30")
		rect(xmin(spreadRate_sd10), ymin(spreadRate_sd10), xmax(spreadRate_sd10), ymax(spreadRate_sd10), xpd=T, lwd=0.2, border="gray30")
		plot(spreadRate_sd10_truncated, legend.only=T, add=T, col=cols2, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.895,0.910,0.06,0.941), 
		     alpha=1, legend.args=list(text="", cex=0.5, line=0.5, col="gray30"), axis.args=list(cex.axis=0.6, lwd=0,
		     lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.45,0), at=seq(0,1,0.1), labels=c(seq(0,0.9,0.1),">=1")))
	}
	
	# C.4. Investigating the impact of environmental factors on the dispersal

segments_obs = c(); index0 = which(data2[,"days"]==0)
for (i in 1:dim(data2)[1])
	{
		if (i != index0)
			{
				segment = cbind(data2[index0,"x_transformed"],data2[index0,"y_transformed"])
				segment = cbind(segment, data2[i,"x_transformed"], data2[i,"y_transformed"])
				segment = cbind(segment, data2[index0,"days"], data2[i,"days"])
				segments_obs = rbind(segments_obs, cbind(segment, data2[i,"days"]-data2[index0,"days"]))
			}
	}
colnames(segments_obs) = c("x1","y1","x2","y2","d1","d2","dt")

		# C.4.1. Analysing the impact of factors on the wavefront dispersal velocity

environmentalExtraction4 = function(segments, rasters, resistances, pathModels, nberOfCorses, ID="obs")
	{
		if (!is.list(rasters)) nberOfColumns = 1
		if (is.list(rasters)) nberOfColumns = length(rasters)
		extractions = matrix(nrow=dim(segments)[1], ncol=nberOfColumns)
		fromCoor = matrix(nrow=dim(segments)[1], ncol=2)
		firstCase = t(segments[which(segments[,"d1"]==0)[1],c("x1","y1")])
		for (i in 1:dim(fromCoor)[1]) fromCoor[i,] = firstCase
		toCoor = segments[,c("x2","y2")]; registerDoMC(cores=nberOfCorses)
		for (i in 1:length(pathModels))
			{
				if (pathModels[i] == 1)
					{
						for (j in 1:dim(segments)[1])
							{
								point1 = segments[which(segments[,"d1"]==0)[1],c("x1","y1")]
								point2 = segments[j,c("x2","y2")]; names(point2) = c("x","y")
								points = rbind(point1,point2); linesList = list()
								linesList[[1]] = Lines(list(Line(points)),1)
								spatialLine = SpatialLines(linesList)
								for (k in 1:length(rasters))
									{
										values = raster::extract(rasters[[k]], spatialLine)[[1]]
										if (resistances[k] == FALSE) values = 1/values
										column = ((i-1)*length(pathModels))+k
										extractions[j,column] = sum(values, na.rm=T)
									}
							}
					}
				if (pathModels[i] == 2)
					{
						buffer = list()
						buffer = foreach(j = 1:length(rasters)) %dopar% {
						# for (j in 1:length(rasters)) {
								if (resistances[j] == FALSE)
									{
										trEnvVariable = transition(rasters[[j]], mean, directions=8)
									}	else	{
										trEnvVariable = transition(rasters[[j]], function(x) 1/mean(x), directions=8)
									}
								trEnvVariableCorr = geoCorrection(trEnvVariable, type="c", multpl=F, scl=T)
								envDistances = costDistance(trEnvVariableCorr, fromCoor[1,], toCoor)
								# column = ((i-1)*length(pathModels))+j; extractions[,column] = t(envDistances)
								envDistances
							}
						for (j in 1:length(buffer))
							{
								column = ((i-1)*length(pathModels))+j; extractions[,column] = buffer[[j]]
							}
					}
				if (pathModels[i] == 3)
					{
						source("circuitscapeFct.r")
						dir.create(file.path("Circuitscape_temp"), showWarnings=F)
						wd1 = getwd(); wd2 = paste0(wd1,"/Circuitscape_temp"); setwd(wd2); buffer = list()
						buffer = foreach(j = 1:length(rasters)) %dopar% {
						# for (j in 1:length(rasters)) {
								rasterName = paste0(names(rasters[[j]]),"_",j)
								if (!file.exists(paste0(rasterName,".asc")))
									{
										writeRaster(rasters[[j]], paste0(rasterName,".asc"), overwrite=T, showWarnings=F)
									}
								envDistances = circuitscapeFct(rasters[[j]], rasterName, resistances[j], resistances[j], 
															   fourCells=F, fromCoor, toCoor, OS="Unix", rasterName, ID)
								# column = ((i-1)*length(pathModels))+j; extractions[,column] = envDistances
								envDistances
							}
						for (j in 1:length(rasters))
							{
								rasterName = paste0(names(rasters[[j]]),"_",j)
								folder = paste(rasterName, "_CStemp_", ID, sep="")		
								# if (file.exists(paste0(folder,"/raster_file_temp_resistances.txt")))
									# {
										tab = read.table(paste0(folder,"/raster_file_temp_resistances.txt"), header=F)
										tab = tab[2:dim(tab)[1], 2:dim(tab)[2]]
										mat = matrix(nrow=dim(as.matrix(fromCoor))[1], ncol=1)
										sameCoordinates = FALSE
										if (sum(fromCoor-toCoor) == 0) sameCoordinates = TRUE
										if (sameCoordinates == FALSE)
											{
												for (k in 1:length(fromCoor[,1])) mat[k] = tab[k,(k+length(fromCoor[,1]))]
											}	else	{
												mat = matrix(nrow=dim(as.matrix(fromCoor))[1], ncol=dim(as.matrix(fromCoor))[1])
												mat = tab[1:dim(fromCoor)[1],1:dim(fromCoor)[1]]
											}
										column = ((i-1)*length(pathModels))+j; extractions[,column] = mat
									# }	else	{
										# column = ((i-1)*length(pathModels))+j; extractions[,column] = buffer[[j]]
									# }
							}
						setwd(wd1)
					}
			}
		colNames = c()
		for (i in 1:length(pathModels))
			{
				if (pathModels[i] == 1) pathModel = "SL"
				if (pathModels[i] == 2) pathModel = "LC"
				if (pathModels[i] == 3) pathModel = "CS"
				for (j in 1:length(rasters)) colNames = c(colNames, paste0(pathModel,"_",names(rasters[[j]])))
			}
		colnames(extractions) = colNames
		return(extractions)
	}

rasters4 = list(); resistances = c(); c = 0
null_raster = forest_areas
names(null_raster) = "null_raster"
null_raster[!is.na(null_raster[])] = 1
c = c + 1; rasters4[[c]] = null_raster
resistances = c(resistances, TRUE)
kS = c(10, 100, 1000)
for (k in kS)
	{
		c = c + 1
		r = forest_areas
		r[!is.na(r[])] = 1+(k*(r[!is.na(r[])]))
		names(r) = paste0("forest_areas_k",k)
		rasters4[[c]] = r
		resistances = c(resistances, FALSE)
	}
r1 = rasterize(motorways, null_raster, fun=min)
r1[!is.na(r1[])] = 1
for (k in kS)
	{
		c = c + 1
		r2 = forest_areas
		r2[!is.na(r[])] = 1+(k*(r2[!is.na(r[])]))
		r2[r1[]==1] = 1/(k*1000)
		names(r2) = paste0("forest_motorways_k",k)
		rasters4[[c]] = r2
		resistances = c(resistances, FALSE)
	}
for (k in kS)
	{
		c = c + 1
		r = rasterize(fences_2, null_raster, fun=min)
		r[!is.na(r[])] = 1; r[(is.na(r[]))&(!is.na(null_raster[]))] = 0
		r[!is.na(r[])] = 1+(k*(r[!is.na(r[])]))
		names(r) = paste0("barriers_k",k)
		rasters4[[c]] = r
		resistances = c(resistances, TRUE)
	}

pathModels = c(2,3); nberOfCorses = 10
extractions_obs_4 = environmentalExtraction4(segments_obs, rasters4, resistances, pathModels, nberOfCorses)
write.csv(extractions_obs_4, "Extractions_LC.csv", row.names=F, quote=F)

pathModel = "LC"; pathModel = "CS"
distancesToFirstCase = matrix(nrow=dim(segments_obs)[1], ncol=1)
firstCasePosition = segments_obs[which(segments_obs[,"d1"]==0)[1],c("x1","y1")]
x1 = firstCasePosition["x1"]; y1 = firstCasePosition["y1"]
for (i in 1:dim(distancesToFirstCase)[1])
	{
		dx = segments_obs[i,"x2"]-firstCasePosition["x1"]
		dy = segments_obs[i,"y2"]-firstCasePosition["y1"]
		distancesToFirstCase[i,1] = sqrt((dx^2)+(dy^2))
	}
extractions_obs_4 = read.csv(paste0("Extractions_",pathModel,".csv"), header=T)
velocities = distancesToFirstCase/segments_obs[,"d2"] # m/day
residualsLR = matrix(nrow=dim(segments_obs)[1], ncol=dim(extractions_obs_4)[2])
colnames(residualsLR) = colnames(extractions_obs_4)
for (i in 2:dim(extractions_obs_4)[2])
	{
		y = segments_obs[,"d2"]; x = extractions_obs_4[,1]
		LR = lm(as.formula(paste0("y ~ x")))
		R2_nul = summary(LR)$r.squared
		y = segments_obs[,"d2"]; x = extractions_obs_4[,i]
		LR = lm(as.formula(paste0("y ~ x")))
		R2_env = summary(LR)$r.squared
		cat(paste0("Q statistic for ",colnames(extractions_obs_4)[i],":\t"))
		cat(round(R2_env-R2_nul,3)); cat("\n")
		y = extractions_obs_4[,i]; x = extractions_obs_4[,1]
		LR = lm(as.formula(paste0("y ~ x")))
		residualsLR[,i] = LR$residuals
	}
if (showingPlots == TRUE)
	{
		dev.new(width=9, height=4); par(mfrow=c(1,2), mgp=c(1,0.35,0), oma=c(0.3,0.5,2.0,2.5), mar=c(3.3,3.3,0,0.5))
		plot(extractions_obs_4[,1], segments_obs[,"d2"], pch=1, cex=0.75, lwd=0.5, axes=F, ann=F, main=NULL, frame=F, col="gray30")
		x = extractions_obs_4[,1]; y = segments_obs[,"d2"]; LR_obs = lm(as.formula(paste0("y ~ x"))); abline(LR_obs, lwd=1.0, col="red")
		axis(side=1, lwd.tick=0.25, cex.axis=0.6, mgp=c(0,0.05,0), lwd=0.25, tck=-0.020, col="gray30", col.tick="gray30", col.axis="gray30") # at=seq(0,300,50))
		axis(side=2, lwd.tick=0.25, cex.axis=0.6, mgp=c(0,0.25,0), lwd=0.25, tck=-0.015, col="gray30", col.tick="gray30", col.axis="gray30", at=seq(0,350,50))
		title(xlab="environmental distances computed on the null raster", cex.lab=0.65, mgp=c(1.0,0,0), col.lab="gray30")
		title(ylab="duration since the first reported case (in days)", cex.lab=0.65, mgp=c(1.5,0,0), col.lab="gray30")
		R2_env_obs = round(summary(LR_obs)$r.squared,5); f = summary(LR_obs)$fstatistic
		p = pf(f[1],f[2],f[3], lower.tail=F); attributes(p) = NULL; pV_env_obs = round(p,3)
		mtext(paste0	("R2 = ",R2_env_obs," (p-value = ",pV_env_obs,")"), side=3, line=-2, col="red", cex=0.7)
		plot(extractions_obs_4[,2], segments_obs[,"d2"], pch=1, cex=0.75, lwd=0.5, axes=F, ann=F, main=NULL, frame=F, col="gray30")
		x = extractions_obs_4[,2]; y = segments_obs[,"d2"]; LR_obs = lm(as.formula(paste0("y ~ x"))); abline(LR_obs, lwd=1.0, col="red")
		axis(side=1, lwd.tick=0.25, cex.axis=0.6, mgp=c(0,0.05,0), lwd=0.25, tck=-0.020, col="gray30", col.lab="gray30", col.axis="gray30") # at=seq(0,60,10))
		axis(side=2, lwd.tick=0.25, cex.axis=0.6, mgp=c(0,0.25,0), lwd=0.25, tck=-0.015, col="gray30", col.lab="gray30", col.axis="gray30", at=seq(0,350,50))
		title(xlab="environmental distances computed on the forest coverage raster (k=10)", cex.lab=0.65, mgp=c(1.0,0,0), col.lab="gray30")
		title(ylab="duration since the first reported case (in days)", cex.lab=0.65, mgp=c(1.5,0,0), col.lab="gray30")
		R2_env_obs = round(summary(LR_obs)$r.squared,5); f = summary(LR_obs)$fstatistic
		p = pf(f[1],f[2],f[3], lower.tail=F); attributes(p) = NULL; pV_env_obs = round(p,3)
		mtext(paste0("R2 = ",R2_env_obs," (p-value = ",pV_env_obs,")"), side=3, line=-2, col="red", cex=0.7)
	}
y = segments_obs[,"d2"]; x1 = extractions_obs_4[,paste0(pathModel,"_null_raster")]
x2 = extractions_obs_4[,paste0(pathModel,"_forest_areas_k10")]
x3 = extractions_obs_4[,paste0(pathModel,"_barriers_k100")]
glm1 = glm(as.formula("y ~ x1 + x2 + x3")); summary(glm1)
glm2 = glm(as.formula("y ~ x2 + x3")); summary(glm2)
glm3 = glm(as.formula("y ~ x2")); summary(glm3)
glm4 = glm(as.formula("y ~ x3")); summary(glm4)
glm5 = glm(as.formula("y ~ x1")); summary(glm5)
library(lmtest); lrtest(glm2, glm3); lrtest(glm4, glm5)

distancesToFirstCase = matrix(nrow=dim(segments_obs)[1], ncol=1)
firstCasePosition = segments_obs[which(segments_obs[,"d1"]==0)[1],c("x1","y1")]
x1 = firstCasePosition["x1"]; y1 = firstCasePosition["y1"]
nberOfNonForestOccurrences = sum(raster::extract(forest_areas,segments_obs[,c("x2","y2")])!=1)
nberOfRandomisations = 1000; # print((nberOfNonForestOccurrences/dim(segments_obs)[1])*100) # 24.59%
for (i in 1:nberOfRandomisations)
	{
		rotations = segments_obs; # print(c(i))
		nonForestRotations = sample(1:dim(rotations)[1],nberOfNonForestOccurrences,replace=F)
		for (j in 1:dim(rotations)[1])
			{
				if (j%in%nonForestRotations)
					{
						nonForestRotation = FALSE; counter = 0
						while (nonForestRotation == FALSE)
							{
								x2 = segments_obs[j,"x2"]; y2 = segments_obs[j,"y2"]
								angle = (2*pi)*runif(1); s = sin(angle); c = cos(angle)
								x = x2-x1; y = y2-y1
								x_new = (x*c)-(y*s); y_new = (x*s)+(y*c)
								x_new = x_new+x1; y_new = y_new+y1
								if (raster::extract(forest_areas,cbind(x_new,y_new)) != 1)
									{
										nonForestRotation = TRUE
										rotations[j,"x2"] = x_new; rotations[j,"y2"] = y_new
									}	else	{
										counter = counter+1
										if (counter == 100)
											{
												nonForestRotation = TRUE; print(c(i,counter))
											}
									}
							}
					}	else	{
						nonForestRotation = TRUE; counter = 0
						while (nonForestRotation == TRUE)
							{
								x2 = segments_obs[j,"x2"]; y2 = segments_obs[j,"y2"]
								angle = (2*pi)*runif(1); s = sin(angle); c = cos(angle)
								x = x2-x1; y = y2-y1
								x_new = (x*c)-(y*s); y_new = (x*s)+(y*c)
								x_new = x_new+x1; y_new = y_new+y1
								if (raster::extract(forest_areas,cbind(x_new,y_new)) == 1)
									{
										nonForestRotation = FALSE
										rotations[j,"x2"] = x_new; rotations[j,"y2"] = y_new
									}	else	{
										counter = counter+1
										if (counter == 100)
											{
												nonForestRotation = FALSE; print(c(i,j,counter))
											}
									}
							}
					}
			}					
		write.csv(rotations, paste0("Rotated_segments/Rotations_",i,".csv"), row.names=F, quote=F)
	}

nberOfRandomisations = 200; rasters4_selected = list(); resistances = rep(NA,3)
rasters4_selected[[1]] = rasters4[[1]]; resistances[1] = TRUE # LC_null_raster
rasters4_selected[[2]] = rasters4[[5]]; resistances[2] = FALSE # LC_forest_motorways_k10
rasters4_selected[[3]] = rasters4[[9]]; resistances[3] = TRUE # LC_barriers_k100
registerDoMC(10); buffer = list()
buffer = foreach(i = 1:nberOfRandomisations) %dopar% {
# for (i in 1:nberOfRandomisations) {
		pathModels = c(2); nberOfCorses = 1; segments_ran = as.matrix(read.csv(paste0("Rotated_segments/Rotations_",i,".csv"), head=T))
		extractions_ran_4 = environmentalExtraction4(segments_ran, rasters4_selected, resistances, pathModels, nberOfCorses, ID=i)
		write.csv(extractions_ran_4, paste0("Rotated_segments/Extractions_",i,"_LC.csv"), row.names=F, quote=F)
		i
	}
nberOfRandomisations = 200
rasters4_selected = list(); resistances = rep(NA,3)
rasters4_selected[[1]] = rasters4[[1]]; resistances[1] = TRUE # CS_null_raster
rasters4_selected[[2]] = rasters4[[5]]; resistances[2] = FALSE # CS_forest_motorways_k10
rasters4_selected[[3]] = rasters4[[8]]; resistances[3] = TRUE # CS_barriers_k10
buffer = foreach(i = 1:nberOfRandomisations) %dopar% {
# for (i in 1:nberOfRandomisations) {
		pathModels = c(3); nberOfCorses = 1; segments_ran = as.matrix(read.csv(paste0("Rotated_segments/Rotations_",i,".csv"), head=T))
		extractions_ran_4 = environmentalExtraction4(segments_ran, rasters4_selected, resistances, pathModels, nberOfCorses, ID=i)
		write.csv(extractions_ran_4, paste0("Rotated_segments/Extractions_",i,"_CS.csv"), row.names=F, quote=F)
		i
	}
for (i in 1:nberOfRandomisations)
	{
		file.rename(paste0("Rotated_segments/Extractions_",i,"_LC.csv"), paste0("Rotated_segments/Extractions_",i,"_LC1.csv"))
		file.rename(paste0("Rotated_segments/Extractions_",i,"_LC_t.csv"), paste0("Rotated_segments/Extractions_",i,"_LC2.csv"))
	}
for (i in 1:nberOfRandomisations)
	{
		tab1 = read.csv(paste0("Rotated_segments/Extractions_",i,"_LC1.csv"), head=T)
		tab2 = read.csv(paste0("Rotated_segments/Extractions_",i,"_LC2.csv"), head=T)
		tab = cbind(tab1[,1:2],tab2[,1],tab1[,3])
		colnames(tab) = c(colnames(tab1)[1:2],colnames(tab2)[1],colnames(tab1)[3])
		write.csv(tab, paste0("Rotated_segments/Extractions_",i,"_LC.csv"), row.names=F, quote=F)
	}

if (showingPlots == TRUE)
	{
		forest_areas_cropped = crop(forest_areas, e2_lam93)
		land_covers_3_cropped = crop(land_covers_3, e2_lam93)
		belgium_cropped = crop(belgium, e2_lam93)
		borders_cropped = crop(borders, e2_lam93)
		fences_1_cropped = crop(fences_1, e2_lam93)
		motorways_cropped = crop(motorways, e2_lam93)
		primaries_cropped = crop(primaries, e2_lam93)
		dev.new(width=9.3, height=4)
		par(mfrow=c(1,2), mar=c(0.5,0.5,0.5,2.0), oma=c(1.0,1.7,1.0,1.3), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		cols1 = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(161)[21:121])
		cols2 = rev(colorRampPalette(brewer.pal(11,"RdYlGn"))(181)[21:121])
		plot(land_covers_3_cropped, col=c(rgb(0,1,0,0),rgb(0,1,0,0),rgb(0,1,0,0),"gray80"), box=F, axes=F, legend=F)
		plot(primaries_cropped, add=T, lwd=0.5, col="gray80")
		# plot(borders_cropped, add=T, lwd=3, col="white", lty=1)
		# plot(borders_cropped, add=T, lwd=0.5, col="blue")
		plot(fences_1_cropped, add=T, lwd=0.75, col="gray30", lty=1)
		plot(motorways_cropped, add=T, lwd=0.75, col="red", lty=1)
		firstCasePosition = segments_obs[which(segments_obs[,"d1"]==0)[1],c("x1","y1")]
		x1 = firstCasePosition["x1"]; y1 = firstCasePosition["y1"]
		colsP = cols1[(((segments_ran[,"d2"]-min(data1[,"days"]))/(max(data1[,"days"])-min(data1[,"days"])))*100)+1]
		for (i in dim(segments_obs)[1]:1)
			{
				x2 = segments_obs[i,"x2"]; y2 = segments_obs[i,"y2"]
				lines(cbind(x1,x2), cbind(y1,y2), col="gray30", lwd=0.1)
			}
		for (i in dim(segments_obs)[1]:1)
			{
				points(t(segments_obs[i,c("x2","y2")]), pch=16, cex=0.7, col=colsP[i])
				points(t(segments_obs[i,c("x2","y2")]), pch=1, cex=0.7, col="gray30", lwd=0.2)
			}
		legend(891000, 6940800, c("Fences","Motorways"), lwd=c(0.75,0.75), cex=0.75, col=c("gray30","red"), text.col=c("gray30","red"), border=NA, x.intersp=0.7, y.intersp=0.8, bty="n")
		rect(xmin(tps_mask), ymin(tps_mask), xmax(tps_mask), ymax(tps_mask), xpd=T, lwd=0.2, border="gray30")
		legendRast = raster(as.matrix(seq(0,max(data1[,"days"]),1)))
		plot(legendRast, legend.only=T, col=cols1, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.895,0.910,0.06,0.941),
		     alpha=1, legend.args=list(text="", cex=0.5, line=0.5, col="gray30"), axis.args=list(cex.axis=0.6, lwd=0,
		     lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.45,0), at=seq(0,360,30)))
		plot(land_covers_3_cropped, col=c(rgb(0,1,0,0),rgb(68/255,165/255,68/255,0.2),rgb(0,1,0,0),rgb(0,1,0,0)), box=F, axes=F, legend=F)
		# plot(primaries_cropped, add=T, lwd=0.5, col="gray80")
		# plot(borders_cropped, add=T, lwd=3, col="white", lty=1)
		# plot(borders_cropped, add=T, lwd=0.5, col="blue")
		# plot(fences_1_cropped, add=T, lwd=0.75, col="gray30", lty=1)
		plot(motorways_cropped, add=T, lwd=0.75, col="red", lty=1)
		segments_ran = as.matrix(read.csv(paste0("Rotated_segments/Rotations_1.csv"), head=T)) # not used
		colsP = cols1[(((segments_obs[,"d2"]-min(data1[,"days"]))/(max(data1[,"days"])-min(data1[,"days"])))*100)+1]
		for (i in dim(segments_obs)[1]:1)
			{
				x2 = segments_obs[i,"x2"]; y2 = segments_obs[i,"y2"]
				lines(cbind(x1,x2), cbind(y1,y2), col="gray30", lwd=0.1)
			}
		for (i in dim(segments_obs)[1]:1)
			{
				points(t(segments_obs[i,c("x2","y2")]), pch=16, cex=0.7, col=colsP[i])
				points(t(segments_obs[i,c("x2","y2")]), pch=1, cex=0.7, col="gray30", lwd=0.2)
			}
		rect(xmin(tps_mask), ymin(tps_mask), xmax(tps_mask), ymax(tps_mask), xpd=T, lwd=0.2, border="gray30")
		legendRast = raster(as.matrix(seq(0,max(data1[,"days"]),1)))
		plot(legendRast, legend.only=T, col=cols1, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.895,0.910,0.06,0.941),
		     alpha=1, legend.args=list(text="", cex=0.5, line=0.5, col="gray30"), axis.args=list(cex.axis=0.6, lwd=0,
		     lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.45,0), at=seq(0,360,30)))
	}

pathModels = c(2,3); Qs_obs = list(); Qs_ran = list(); n = 0
for (i in 2:length(rasters4_selected))
	{
		for (j in 1:length(pathModels))
			{
				if (pathModels[j] == 2) { pathModel = "LC" }
				if (pathModels[j] == 3) { pathModel = "CS" }
				envVariableName = names(rasters4_selected[[i]])
				if ((i==3)&(j==1)) envVariableName = gsub("k10","k100",envVariableName)
				extractions_obs_4 = read.csv(paste0("Extractions_",pathModel,".csv"), header=T)
				rasterNull = paste0(pathModel,"_null_raster")
				rasterName = paste0(pathModel,"_",envVariableName)
				y = segments_obs[,"d2"]; x = extractions_obs_4[,rasterNull]
				LR = lm(as.formula(paste0("y ~ x")))
				R2_obs_nul = summary(LR)$r.squared
				y = segments_obs[,"d2"]; x = extractions_obs_4[,rasterName]
				LR = lm(as.formula(paste0("y ~ x")))
				R2_obs_env = summary(LR)$r.squared
				Q_obs = R2_obs_env-R2_obs_nul; c = 0
				Q_rans = rep(NA, nberOfRandomisations)
				for (k in 1:nberOfRandomisations)
					{
						segments_ran = read.csv(paste0("Rotated_segments/Rotations_",k,".csv"), header=T)
						extractions_ran_4 = read.csv(paste0("Rotated_segments/Extractions_",k,"_",pathModel,".csv"))
						y = segments_ran[,"d2"]; x = extractions_ran_4[,rasterNull]
						LR = lm(as.formula(paste0("y ~ x")))
						R2_ran_nul = summary(LR)$r.squared
						y = segments_ran[,"d2"]; x = extractions_ran_4[,rasterName]
						LR = lm(as.formula(paste0("y ~ x")))
						R2_ran_env = summary(LR)$r.squared
						Q_ran = R2_ran_env-R2_ran_nul; Q_rans[k] = Q_ran
						if (Q_obs < Q_ran) c = c+1
					}
				pValue = round(c/nberOfRandomisations,3)
				cat(pathModel," - ",rasterName,", Qobs = ",round(Q_obs,3),", p-value = ",pValue,"\n", sep="")
				n = n+1; Qs_obs[[n]] = Q_obs; Qs_ran[[n]] = Q_rans
			}
	}
if (showingPlots == TRUE)
	{
		dev.new(width=8, height=3); par(mfrow=c(2,2), mgp=c(1,0.35,0), oma=c(0.5,1.3,1,2), mar=c(2.5,2,0.5,0), lwd=0.2)
		cols0a = rgb(120,120,120,255,maxColorValue=255); cols0b = rgb(120,120,120,100,maxColorValue=255) # grey
		cols1a = rgb(4,153,14,255,maxColorValue=255); cols1b = rgb(4,153,14,100,maxColorValue=255) # greens
		cols2a = rgb(204,0,0,255,maxColorValue=255); cols2b = rgb(204,0,0,100,maxColorValue=255) # red
		plot(density(Qs_ran[[1]]), lwd=0.7, col=cols0a, axes=F, ann=F, xlim=c(-0.2,0.2), ylim=c(0,35))
		polygon(density(Qs_ran[[1]]), col=cols0b, border=NA); abline(v=Qs_obs[[1]], lwd=1.2, col="gray30", lty=2)
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,-0.03,0), lwd=0.2, tck=-0.035, col="gray30", col.tick="gray30", col.axis="gray30", at=seq(-0.25,0.25,0.05))
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.20,0), lwd=0.2, tck=-0.035, col="gray30", col.tick="gray30", col.axis="gray30", at=seq(0,70,10))
		title(ylab="density", cex.lab=0.7, mgp=c(1.2,0,0), col.lab="gray30")		
		plot(density(Qs_ran[[3]]), lwd=0.7, col=cols0a, axes=F, ann=F, xlim=c(-0.2,0.2), ylim=c(0,35))
		polygon(density(Qs_ran[[3]]), col=cols0b, border=NA); abline(v=Qs_obs[[3]], lwd=1.2, col="gray30", lty=2)
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,-0.03,0), lwd=0.2, tck=-0.035, col="gray30", col.tick="gray30", col.axis="gray30", at=seq(-0.25,0.25,0.05))
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.20,0), lwd=0.2, tck=-0.035, col="gray30", col.tick="gray30", col.axis="gray30", at=seq(0,70,10))
		plot(density(Qs_ran[[2]]), lwd=0.7, col=cols0a, axes=F, ann=F, xlim=c(-0.2,0.2), ylim=c(0,35))
		polygon(density(Qs_ran[[2]]), col=cols0b, border=NA); abline(v=Qs_obs[[2]], lwd=1.2, col="gray30", lty=2)
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,-0.03,0), lwd=0.2, tck=-0.035, col="gray30", col.tick="gray30", col.axis="gray30", at=seq(-0.25,0.25,0.05))
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.20,0), lwd=0.2, tck=-0.035, col="gray30", col.tick="gray30", col.axis="gray30", at=seq(0,70,10))
		title(xlab=expression(italic(Q) == {R^{2}}[env] - {R^{2}}[null]), cex.lab=0.7, mgp=c(1.0,0,0), col.lab="gray30")
		title(ylab="density", cex.lab=0.7, mgp=c(1.2,0,0), col.lab="gray30")
		plot(density(Qs_ran[[4]]), lwd=0.7, col=cols0a, axes=F, ann=F, xlim=c(-0.2,0.2), ylim=c(0,35))
		polygon(density(Qs_ran[[4]]), col=cols0b, border=NA); abline(v=Qs_obs[[4]], lwd=1.2, col="gray30", lty=2)
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,-0.03,0), lwd=0.2, tck=-0.035, col="gray30", col.tick="gray30", col.axis="gray30", at=seq(-0.25,0.25,0.05))
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.20,0), lwd=0.2, tck=-0.035, col="gray30", col.tick="gray30", col.axis="gray30", at=seq(0,70,10))
		title(xlab=expression(italic(Q) == {R^{2}}[env] - {R^{2}}[null]), cex.lab=0.7, mgp=c(1.0,0,0), col.lab="gray30")
	}

if (showingPlots == TRUE)
	{
		forest_areas_cropped = crop(forest_areas, e2_lam93)
		land_covers_3_cropped = crop(land_covers_3, e2_lam93)
		belgium_cropped = crop(belgium, e2_lam93)
		borders_cropped = crop(borders, e2_lam93)
		fences_1_cropped = crop(fences_1, e2_lam93)
		motorways_cropped = crop(motorways, e2_lam93)
		primaries_cropped = crop(primaries, e2_lam93)
		currentMap1 = raster("Circuitscape_maps/CS_current_maps/Segments_null_raster.asc")
		currentMap1 = crop(currentMap1, e2_lam93); currentMap1[] = log10(currentMap1[]+1); currentMap1[currentMap1[]>1] = 1
		currentMap2 = raster("Circuitscape_maps/CS_current_maps/Segments_forest_motor.asc")
		currentMap2 = crop(currentMap2, e2_lam93); currentMap2[] = log10(currentMap2[]+1); currentMap2[currentMap2[]>1] = 1
		dev.new(width=9.3, height=4); par(mfrow=c(1,2), mar=c(0.5,0.5,0.5,2.0), oma=c(1.0,1.7,1.0,1.3), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		cols = c(rep("#FFFFFF",7),colorRampPalette(brewer.pal(9,"YlOrBr"))(94))
		plot(currentMap1, main="", cex.main=1, cex.axis=0.7, bty="n", box=F, axes=F, legend=F, axis.args=list(cex.axis=0.7), col=cols, colNA="white")
		# plot(land_covers_3_cropped, add=T, col=c(rgb(0,1,0,0),rgb(68/255,165/255,68/255,0.0),rgb(0,1,0,0),"gray80"), box=F, axes=F, legend=F)
		# plot(primaries_cropped, add=T, lwd=0.5, col="gray80")
		points(data2[,c("x_transformed","y_transformed")], pch=3, cex=0.5, lwd=0.4, col="gray30")
		# plot(borders_cropped, add=T, lwd=3, col="white", lty=1)
		# plot(borders_cropped, add=T, lwd=0.5, col="blue")
		# plot(fences_1_cropped, add=T, lwd=0.75, col="gray30", lty=1)
		# plot(motorways_cropped, add=T, lwd=0.75, col="red", lty=1)
		mtext("1. Circuitscape current map", side=3, line=-2, at=899000, cex=0.8, font=1, col="gray30")
		mtext("   (based on the null raster)", side=3, line=-2.8, at=899000, cex=0.8, font=1, col="gray30")
		rect(xmin(tps_mask), ymin(tps_mask), xmax(tps_mask), ymax(tps_mask), xpd=T, lwd=0.2, border="gray30")
		plot(currentMap1, legend.only=T, add=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.895,0.910,0.06,0.941),
		     alpha=1, legend.args=list(text="", cex=0.5, line=0.5, col="gray30"), axis.args=list(cex.axis=0.6, lwd=0,
		     lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.45,0)))
		plot(currentMap2, main="", cex.main=1, cex.axis=0.7, bty="n", box=F, axes=F, legend=F, axis.args=list(cex.axis=0.7), col=cols, colNA="white")
		# plot(land_covers_3_cropped, add=T, col=c(rgb(0,1,0,0),rgb(68/255,165/255,68/255,0.0),rgb(0,1,0,0),"gray80"), box=F, axes=F, legend=F)
		# plot(primaries_cropped, add=T, lwd=0.5, col="gray80")
		points(data2[,c("x_transformed","y_transformed")], pch=3, cex=0.5, lwd=0.4, col="gray30")
		# plot(borders_cropped, add=T, lwd=3, col="white", lty=1)
		# plot(borders_cropped, add=T, lwd=0.5, col="blue")
		# plot(fences_1_cropped, add=T, lwd=0.75, col="gray30", lty=1)
		# plot(motorways_cropped, add=T, lwd=0.75, col="red", lty=1)
		mtext("2. Circuitscape current map", side=3, line=-2, at=899000, cex=0.8, font=1, col="gray30")
		mtext("(based on the forest areas)", side=3, line=-2.8, at=899000, cex=0.8, font=1, col="gray30")
		rect(xmin(tps_mask), ymin(tps_mask), xmax(tps_mask), ymax(tps_mask), xpd=T, lwd=0.2, border="gray30")
		plot(currentMap2, legend.only=T, add=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.895,0.910,0.06,0.941),
		     alpha=1, legend.args=list(text="", cex=0.5, line=0.5, col="gray30"), axis.args=list(cex.axis=0.6, lwd=0,
		     lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.45,0)))
	}

		# C.4.2. Analysis of the impact of barriers on the dispersal frequency
	
			# C.4.2.1. Determining the number of edges crossing a barrier

countingCrossingBarrierEvents2 = function(segments, shapefiles)
	{
		I = which(segments[,"d1"]==0)[1]
		crossingBarrierEvents = matrix(nrow=dim(segments)[1], ncol=length(shapefiles))
		for (i in 1:dim(segments)[1])
			{
				point1 = segments[I,c("x1","y1")]; names(point1) = c("x","y")
				point2 = segments[i,c("x2","y2")]; names(point2) = c("x","y")
				points = rbind(point1,point2); linesList = list()
				linesList[[1]] = Lines(list(Line(points)),1)
				spatialLine = SpatialLines(linesList)
				crs(spatialLine) = crs(shapefiles[[1]])
				for (j in 1:length(shapefiles))
					{
						crossingBarrierEvent = 0
						intersections = gIntersection(shapefiles[[j]], spatialLine)
						if (!is.null(intersections)) crossingBarrierEvent = 1
						crossingBarrierEvents[i,j] = crossingBarrierEvent
					}
			}
		return(crossingBarrierEvents)
	}

shapefiles = list(fences_2)
crossingBarriers_obs = countingCrossingBarrierEvents2(segments_obs, shapefiles)
write.csv(crossingBarriers_obs, "Crossing_barriers2.csv", row.names=F, quote=F)

registerDoMC(10); buffer = list()
nberOfRandomisations = 200
buffer = foreach(i = 1:nberOfRandomisations) %dopar% {
# for (i in 1:nberOfRandomisations) {
		segments_ran = read.csv(paste0("Rotated_segments/Rotations_",i,".csv"), header=T)
		crossingBarriers_ran = countingCrossingBarrierEvents2(segments_ran, shapefiles)
		fileName = paste0("Rotated_segments/Crossing_barriers_",i,".csv")
		write.csv(crossingBarriers_ran, fileName, row.names=F, quote=F)
		i
	}

			# C.4.2.2. Statistical test to assess the impact of barriers

fences_1_cropped = crop(fences_1, e2_lam93); motorways_cropped = crop(motorways, e2_lam93)
i = 1; randomisation = read.csv(paste0("Rotated_segments/Rotations_",i,".csv"), header=T)
plot(randomisation[,c("x2","y2")], pch=3, lwd=0.75, cex=0.75, col="gray30", ann=F, axes=F)
plot(fences_1_cropped, add=T, lwd=0.75, col="gray30", lty=1); plot(motorways_cropped, add=T, lwd=0.75, col="red", lty=1)

nberOfRandomisations = 200; Ns_ran = list()
crossingBarriers_obs = read.csv("Crossing_barriers2.csv", head=T)
counters = rep(0, dim(crossingBarriers_obs)[2])
observed_Ns = matrix(nrow=1, ncol=dim(crossingBarriers_obs)[2])
randomised_Ns = matrix(nrow=nberOfRandomisations, ncol=dim(crossingBarriers_obs)[2])
colnames(observed_Ns) = colnames(crossingBarriers_obs)
colnames(randomised_Ns) = colnames(crossingBarriers_obs)
for (i in 1:nberOfRandomisations)
	{
		crossingBarriers_ran = read.csv(paste0("Rotated_segments/Crossing_barriers_",i,".csv"), header=T)
		for (j in 1:dim(crossingBarriers_obs)[2])
			{
				N_obs = sum(crossingBarriers_obs[,j]); observed_Ns[1,j] = N_obs
				N_ran = sum(crossingBarriers_ran[,j]); randomised_Ns[i,j] = N_ran
				if (N_ran <= N_obs) counters[j] = counters[j] + 1
			}
	}
pValue = t(counters/nberOfRandomisations); cat(pValue) # = 0
if (showingPlots == TRUE)
	{
		dev.new(width=4, height=2); par(mfrow=c(1,1), mgp=c(1,0.35,0), oma=c(0.5,1.3,1,2), mar=c(2.5,2,0.5,0), lwd=0.2)
		cols0a = rgb(120,120,120,255,maxColorValue=255); cols0b = rgb(120,120,120,100,maxColorValue=255) # grey
		cols1a = rgb(4,153,14,255,maxColorValue=255); cols1b = rgb(4,153,14,100,maxColorValue=255) # greens
		cols2a = rgb(204,0,0,255,maxColorValue=255); cols2b = rgb(204,0,0,100,maxColorValue=255) # red
		plot(density(randomised_Ns), lwd=0.7, col=cols0a, axes=F, ann=F, xlim=c(0,200))
		polygon(density(randomised_Ns), col=cols0b, border=NA); abline(v=observed_Ns, lwd=1.2, col="gray30", lty=2)
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,-0.03,0), lwd=0.2, tck=-0.035, col="gray30", col.tick="gray30", col.axis="gray30", at=seq(0,200,40))
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.20,0), lwd=0.2, tck=-0.035, col="gray30", col.tick="gray30", col.axis="gray30")
		title(ylab="density", cex.lab=0.7, mgp=c(1.2,0,0), col.lab="gray30")		
	}

