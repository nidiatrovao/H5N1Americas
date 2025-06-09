# TO DO:
	# - estimating and comparing the dispersal statistics estimated for the different genotypes (diffusions velocity and IBD)

library(diagram)
library(ks)
library(lubridate)
library(MetBrewer)
library(seraphim)

genotypes = c("B32","B36","C21","D11")
nberOfExtractionFiles = 1000

	# Note: for B3.6, "Ross's_goose" has to changed to "Rosss_goose" in both the ".csv" and ".trees" files

# 1. Preparing the XML and KML files for the continuous/discrete phylogeographic analyses

	# 1.1. Selecting 1,000 empirical trees

tip_labels = list()
for (h in 1:length(genotypes))
	{
		if (file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_skygrid.trees")))
			{
				if (!file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_empirical.trees")))
					{
						all_trees = scan(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_skygrid.trees"), what="", sep="\n", quiet=T)
						index = which(grepl("\t\t;",all_trees)); index = index[length(index)]; txt = all_trees[1:index]; ratio = 0.1
						if (genotypes[h] == "B36") ratio = 0.2
						indices = which(grepl("tree STATE_",all_trees)); burnIn = round(ratio*(length(indices)))+1
						selected_trees = all_trees[sample(indices[(burnIn+1):length(indices)], nberOfExtractionFiles, replace=F)]
						write(c(txt,selected_trees,"End;"), paste0("Genotype_",genotypes[h],"/",genotypes[h],"_empirical.trees"))
					}
				tip_labels[[h]] = read.nexus(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_empirical.trees"))$tip.label[[1]]
			}	else	{
				if (file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_emp_NT.trees")))
					{
						tip_labels[[h]] = read.nexus(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_emp_NT.trees"))$tip.label[[1]]
					}				
			}
	}

	# 1.2. Generating the metadata file

countries_admin1 = c("Canada","USA","Mexico","Colombia","Honduras","Costa Rica")
if (!file.exists("GADMs_1.rds"))
	{
		gadms1 = list()
		gadms1[[1]] = shapefile("All_shapefiles/GADM_files/GADM_CAN_1.shp")
		gadms1[[2]] = shapefile("All_shapefiles/GADM_files/GADM_USA_1.shp")
		gadms1[[3]] = shapefile("All_shapefiles/GADM_files/GADM_MEX_1.shp")
		gadms1[[4]] = shapefile("All_shapefiles/GADM_files/GADM_COL_1.shp")
		gadms1[[5]] = shapefile("All_shapefiles/GADM_files/GADM_HND_1.shp")
		gadms1[[6]] = shapefile("All_shapefiles/GADM_files/GADM_CRI_1.shp")
		gadms1[[4]]@data[which(gadms1[[4]]@data[,"NAME_1"]=="Bolívar"),"NAME_1"] = "Bolivar"
		gadms1[[4]]@data[which(gadms1[[4]]@data[,"NAME_1"]=="Córdoba"),"NAME_1"] = "Cordoba"
		gadms1[[5]]@data[which(gadms1[[5]]@data[,"NAME_1"]=="Atlántida"),"NAME_1"] = "Atlantida"
		saveRDS(gadms1, "GADMs_1.rds")
	}	else	{
		gadms1 = readRDS("GADMs_1.rds")
	}
if (!file.exists("GADMs_1.rds"))
	{
		gadms2_USA = shapefile("All_shapefiles/GADM_files/GADM_USA_2.shp")
		saveRDS(gadms2_USA, "GADMs_2.rds")
	}	else	{
		gadms2_USA = readRDS("GADMs_2.rds")
	}
for (h in 1:length(genotypes))
	{
		if (genotypes[h] == "B32")
			{
				tab1 = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.csv"), head=T)
				tab2 = tab1[,c("trait","country","province.dept.region","county.district","location_used","latitude","longitude")]
				colnames(tab2) = c("trait","country","admin1","admin2","precision","latitude","longitude")
				tab2[which(tab2[,"admin1"]=="NewYork"),"admin1"] = "New York"
				tab2[which(tab2[,"admin2"]=="Caroll"),"admin2"] = "Carroll"
				tab2[which(tab2[,"admin2"]=="Kingsbury\xa0"),"admin2"] = "Kingsbury"
				tab2[which(tab2[,"admin2"]=="LaSalle"),"admin2"] = "La Salle"
				tab2[which(tab2[,"admin2"]=="Mcintosh"),"admin2"] = "McIntosh"
				tab2[which(tab2[,"admin2"]=="Nome Census Area"),"admin2"] = "Nome"
				tab2[which(tab2[,"admin2"]=="St. Louis"),"admin2"] = "Saint Louis"
				tab2[which(tab2[,"admin2"]=="Valdez-Cordoba"),"admin2"] = "Valdez-Cordova"
				tab2[,c("latitude","longitude")] = NA
				for (i in 1:dim(tab2)[1])
					{
						if (grepl("precise",tab2[i,"precision"])|grepl("http",tab2[i,"precision"])|grepl("city",tab2[i,"precision"])
							|grepl("beach",tab2[i,"precision"])|grepl("district",tab2[i,"precision"]))
							{
								tab2[i,c("latitude","longitude")] = tab1[i,c("latitude","longitude")]
							}
						if (grepl("county",tab2[i,"precision"]))
							{
								maxArea = 0; index3 = 0
								index2 = which((gadms2_USA@data[,"NAME_1"]==gsub("_"," ",tab2[i,"admin1"]))&(gadms2_USA@data[,"NAME_2"]==tab2[i,"admin2"]))
								for (j in 1:length(gadms2_USA@polygons[[index2]]@Polygons))
									{
										if (maxArea < gadms2_USA@polygons[[index2]]@Polygons[[j]]@area)
											{
												maxArea = gadms2_USA@polygons[[index2]]@Polygons[[j]]@area; index3 = j
											}
									}
								pol = gadms2_USA@polygons[[index2]]@Polygons[[index3]]
								if (!file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons/",gsub(" ","",tab2[i,"admin1"]),"_",gsub(" ","",tab2[i,"admin2"]),".kml")))
									{
										sink(file=paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons/",gsub(" ","",tab2[i,"admin1"]),"_",gsub(" ","",tab2[i,"admin2"]),".kml"))
										cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"); cat("\n")
										cat("<kml xmlns=\"http://earth.google.com/kml/2.2\">"); cat("\n")
										cat(paste("\t<polygon id=\"",gsub(" ","",tab2[i,"admin1"]),"_",gsub(" ","",tab2[i,"admin2"]),"\" samplingProbability=\"",1,"\">",sep="")); cat("\n")
										cat("\t\t<coordinates>"); cat("\n")
										for (j in 1:dim(pol@coords)[1])
											{
												cat(paste("\t\t\t",pol@coords[j,2],",",pol@coords[j,1],",0",sep="")); cat("\n")
											}
										cat("\t\t</coordinates>"); cat("\n")
										cat("\t</polygon>"); cat("\n")
										cat("</kml>"); cat("\n")
										sink(NULL)				
									}
								p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
								pol = sps; proj4string(pol) = gadms2_USA@proj4string
								random_point = spsample(pol, 1, type="random")@coords
								tab2[i,"latitude"] = random_point[1,2]; tab2[i,"longitude"] = random_point[1,1]
							}		
						if (grepl("centroid",tab2[i,"precision"]))
							{
								index1 = which(countries_admin1==unlist(strsplit(tab2[i,"country"],"-"))[1]); maxArea = 0; index3 = 0
								index2 = which(grepl(paste0("\\.",unlist(strsplit(tab2[i,"country"],"-"))[2]),gadms1[[index1]]@data[,"HASC_1"]))
								for (j in 1:length(gadms1[[index1]]@polygons[[index2]]@Polygons))
									{
										if (maxArea < gadms1[[index1]]@polygons[[index2]]@Polygons[[j]]@area)
											{
												maxArea = gadms1[[index1]]@polygons[[index2]]@Polygons[[j]]@area; index3 = j
											}
									}
								pol = gadms1[[index1]]@polygons[[index2]]@Polygons[[index3]]
								country = unlist(strsplit(tab2[i,"country"],"-"))[1]
								admin1 = unlist(strsplit(tab2[i,"country"],"-"))[2]
								if (!file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons/",country,"_",admin1,".kml")))
									{
										sink(file=paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons/",country,"_",admin1,".kml"))
										cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"); cat("\n")
										cat("<kml xmlns=\"http://earth.google.com/kml/2.2\">"); cat("\n")
										cat(paste("\t<polygon id=\"",country,"_",admin1,"\" samplingProbability=\"",1,"\">",sep="")); cat("\n")
										cat("\t\t<coordinates>"); cat("\n")
										for (j in 1:dim(pol@coords)[1])
											{
												cat(paste("\t\t\t",pol@coords[j,2],",",pol@coords[j,1],",0",sep="")); cat("\n")
											}
										cat("\t\t</coordinates>"); cat("\n")
										cat("\t</polygon>"); cat("\n")
										cat("</kml>"); cat("\n")
										sink(NULL)				
									}
								p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
								pol = sps; proj4string(pol) = gadms1[[index1]]@proj4string
								random_point = spsample(pol, 1, type="random")@coords
								tab2[i,"latitude"] = random_point[1,2]; tab2[i,"longitude"] = random_point[1,1]
							}
						if (tab2[i,"precision"] == "")
							{
								index1 = which(countries_admin1==tab2[i,"country"]); maxArea = 0; index3 = 0
								index2 = which(grepl(tab2[i,"admin1"],gadms1[[index1]]@data[,"NAME_1"]))
								for (j in 1:length(gadms1[[index1]]@polygons[[index2]]@Polygons))
									{
										if (maxArea < gadms1[[index1]]@polygons[[index2]]@Polygons[[j]]@area)
											{
												maxArea = gadms1[[index1]]@polygons[[index2]]@Polygons[[j]]@area; index3 = j
											}
									}
								pol = gadms1[[index1]]@polygons[[index2]]@Polygons[[index3]]
								if (!file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons/",tab2[i,"country"],"_",tab2[i,"admin1"],".kml")))
									{
										sink(file=paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons/",tab2[i,"country"],"_",tab2[i,"admin1"],".kml"))
										cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"); cat("\n")
										cat("<kml xmlns=\"http://earth.google.com/kml/2.2\">"); cat("\n")
										cat(paste("\t<polygon id=\"",tab2[i,"country"],"_",tab2[i,"admin1"],"\" samplingProbability=\"",1,"\">",sep="")); cat("\n")
										cat("\t\t<coordinates>"); cat("\n")
										for (j in 1:dim(pol@coords)[1])
											{
												cat(paste("\t\t\t",pol@coords[j,2],",",pol@coords[j,1],",0",sep="")); cat("\n")
											}
										cat("\t\t</coordinates>"); cat("\n")
										cat("\t</polygon>"); cat("\n")
										cat("</kml>"); cat("\n")
										sink(NULL)				
									}
								p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
								pol = sps; proj4string(pol) = gadms1[[index1]]@proj4string
								random_point = spsample(pol, 1, type="random")@coords
								tab2[i,"latitude"] = random_point[1,2]; tab2[i,"longitude"] = random_point[1,1]
							}
						if (tab2[i,"admin1"] == "") tab2[i,"admin1"] = NA
						if (tab2[i,"admin2"] == "") tab2[i,"admin2"] = NA
						if (tab2[i,"precision"] == "") tab2[i,"precision"] = NA
					}
				indices = which(tab2[,"trait"]%in%tip_labels[[h]]); tab2 = tab2[indices,]
				if (!file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.txt")))
					{
						write.table(tab2, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.txt"), row.names=F, quote=F, sep="\t")
					}
			}
		if (genotypes[h] == "B36")
			{
				tab1 = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.csv"), head=T)
				tab2 = tab1[,c("label","country","location","county","latitude","longitude")]
				colnames(tab2) = c("trait","country","admin1","admin2","latitude","longitude")
				tab2[which(tab2[,"admin2"]=="DeSoto"),"admin2"] = "Desoto"
				tab2[which(tab2[,"admin2"]=="Golden-Valley"),"admin2"] = "Golden Valley"
				tab2[which(tab2[,"admin2"]=="LaMoure"),"admin2"] = "Lamoure"
				tab2[which(tab2[,"admin2"]=="Lewis"),"admin2"] = "Lewis and Clark"
				tab2[which(tab2[,"admin2"]=="Livingston Parish"),"admin2"] = "Livingston"
				tab2[which(tab2[,"admin2"]=="Noble"),"admin2"] = "Nobles"
				tab2[which(tab2[,"admin2"]=="Norfolk City"),"admin2"] = "Norfolk"
				tab2[which(tab2[,"admin2"]=="Otter"),"admin2"] = "Otter Tail"
				tab2[which(tab2[,"admin2"]=="St. Louis"),"admin2"] = "Saint Louis"
				tab2[,c("latitude","longitude")] = NA
				for (i in 1:dim(tab2)[1])
					{
						if (tab2[i,"admin2"] != "")
							{
								maxArea = 0; index3 = 0
								index2 = which(gadms2_USA@data[,"NAME_1"]==(gsub("_"," ",tab2[i,"admin1"]))&(gadms2_USA@data[,"NAME_2"]==tab2[i,"admin2"]))
								for (j in 1:length(gadms2_USA@polygons[[index2]]@Polygons))
									{
										if (maxArea < gadms2_USA@polygons[[index2]]@Polygons[[j]]@area)
											{
												maxArea = gadms2_USA@polygons[[index2]]@Polygons[[j]]@area; index3 = j
											}
									}
								pol = gadms2_USA@polygons[[index2]]@Polygons[[index3]]
								if (!file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons/",gsub(" ","",tab2[i,"admin1"]),"_",gsub(" ","",tab2[i,"admin2"]),".kml")))
									{
										sink(file=paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons/",gsub(" ","",tab2[i,"admin1"]),"_",gsub(" ","",tab2[i,"admin2"]),".kml"))
										cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"); cat("\n")
										cat("<kml xmlns=\"http://earth.google.com/kml/2.2\">"); cat("\n")
										cat(paste("\t<polygon id=\"",gsub(" ","",tab2[i,"admin1"]),"_",gsub(" ","",tab2[i,"admin2"]),"\" samplingProbability=\"",1,"\">",sep="")); cat("\n")
										cat("\t\t<coordinates>"); cat("\n")
										for (j in 1:dim(pol@coords)[1])
											{
												cat(paste("\t\t\t",pol@coords[j,2],",",pol@coords[j,1],",0",sep="")); cat("\n")
											}
										cat("\t\t</coordinates>"); cat("\n")
										cat("\t</polygon>"); cat("\n")
										cat("</kml>"); cat("\n")
										sink(NULL)				
									}
								p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
								pol = sps; proj4string(pol) = gadms2_USA@proj4string
								random_point = spsample(pol, 1, type="random")@coords
								tab2[i,"latitude"] = random_point[1,2]; tab2[i,"longitude"] = random_point[1,1]
							}
						if (tab2[i,"admin2"] == "")
							{
								if (grepl("-",tab2[i,"country"]))
									{
										index1 = which(countries_admin1==unlist(strsplit(tab2[i,"country"],"-"))[1])
										index2 = which(grepl(paste0("\\.",unlist(strsplit(tab2[i,"country"],"-"))[2]),gadms1[[index1]]@data[,"HASC_1"]))
										if (length(index2) == 0) index2 = which(gsub("CA.","",gadms1[[index1]]@data[,"HASC_1"],"\\.") == tab2[i,"admin1"])
									}	else	{
										index1 = which(countries_admin1==tab2[i,"country"])
										index2 = which(gsub("CA.","",gadms1[[index1]]@data[,"HASC_1"])==tab2[i,"admin1"])
										if (length(index2) == 0)
											{
												index2 = which(gadms1[[index1]]@data[,"NAME_1"]==tab2[i,"admin1"])
											}
									}
								maxArea = 0; index3 = 0
								for (j in 1:length(gadms1[[index1]]@polygons[[index2]]@Polygons))
									{
										if (maxArea < gadms1[[index1]]@polygons[[index2]]@Polygons[[j]]@area)
											{
												maxArea = gadms1[[index1]]@polygons[[index2]]@Polygons[[j]]@area; index3 = j
											}
									}
								pol = gadms1[[index1]]@polygons[[index2]]@Polygons[[index3]]
								if (grepl("-",tab2[i,"country"]))
									{
										country = unlist(strsplit(tab2[i,"country"],"-"))[1]; admin1 = unlist(strsplit(tab2[i,"country"],"-"))[2]
									}	else	{
										country = tab2[i,"country"]; admin1 = tab2[i,"admin1"]
									}
								if (!file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons/",country,"_",admin1,".kml")))
									{
										sink(file=paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons/",country,"_",admin1,".kml"))
										cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"); cat("\n")
										cat("<kml xmlns=\"http://earth.google.com/kml/2.2\">"); cat("\n")
										cat(paste("\t<polygon id=\"",country,"_",admin1,"\" samplingProbability=\"",1,"\">",sep="")); cat("\n")
										cat("\t\t<coordinates>"); cat("\n")
										for (j in 1:dim(pol@coords)[1])
											{
												cat(paste("\t\t\t",pol@coords[j,2],",",pol@coords[j,1],",0",sep="")); cat("\n")
											}
										cat("\t\t</coordinates>"); cat("\n")
										cat("\t</polygon>"); cat("\n")
										cat("</kml>"); cat("\n")
										sink(NULL)				
									}
								p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
								pol = sps; proj4string(pol) = gadms1[[index1]]@proj4string
								random_point = spsample(pol, 1, type="random")@coords
								tab2[i,"latitude"] = random_point[1,2]; tab2[i,"longitude"] = random_point[1,1]
							}
						if (tab2[i,"admin1"] == "") tab2[i,"admin1"] = NA
						if (tab2[i,"admin2"] == "") tab2[i,"admin2"] = NA
					}
				indices = which(tab2[,"trait"]%in%tip_labels[[h]]); tab2 = tab2[indices,]
				if (!file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.txt")))
					{
						write.table(tab2, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.txt"), row.names=F, quote=F, sep="\t")
					}
			}
		if (genotypes[h] == "C21")
			{
				tab1 = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.csv"), head=T)
				tab2 = tab1[,c("trait","state","county","latitude","longitude")]
				colnames(tab2) = c("trait","admin1","admin2","latitude","longitude")
				tab2[which(tab2[,"admin2"]=="Plum_Island"),"admin2"] = "Suffolk"
				tab2[which(tab2[,"admin2"]=="Stafford"),"admin2"] = "Strafford"
				tab2[,c("latitude","longitude")] = NA
				for (i in 1:dim(tab2)[1])
					{
						if (tab2[i,"admin2"] != "")
							{
								maxArea = 0; index3 = 0; admin1 = unlist(strsplit(tab2[i,"admin1"],"-"))[2]
								index2 = which((grepl(paste0("US.",admin1),gadms2_USA@data[,"HASC_2"]))&(gadms2_USA@data[,"NAME_2"]==tab2[i,"admin2"]))
								for (j in 1:length(gadms2_USA@polygons[[index2]]@Polygons))
									{
										if (maxArea < gadms2_USA@polygons[[index2]]@Polygons[[j]]@area)
											{
												maxArea = gadms2_USA@polygons[[index2]]@Polygons[[j]]@area; index3 = j
											}
									}
								pol = gadms2_USA@polygons[[index2]]@Polygons[[index3]]
								if (!file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons/",gsub(" ","",tab2[i,"admin1"]),"_",gsub(" ","",tab2[i,"admin2"]),".kml")))
									{
										sink(file=paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons/",gsub(" ","",tab2[i,"admin1"]),"_",gsub(" ","",tab2[i,"admin2"]),".kml"))
										cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"); cat("\n")
										cat("<kml xmlns=\"http://earth.google.com/kml/2.2\">"); cat("\n")
										cat(paste("\t<polygon id=\"",gsub(" ","",tab2[i,"admin1"]),"_",gsub(" ","",tab2[i,"admin2"]),"\" samplingProbability=\"",1,"\">",sep="")); cat("\n")
										cat("\t\t<coordinates>"); cat("\n")
										for (j in 1:dim(pol@coords)[1])
											{
												cat(paste("\t\t\t",pol@coords[j,2],",",pol@coords[j,1],",0",sep="")); cat("\n")
											}
										cat("\t\t</coordinates>"); cat("\n")
										cat("\t</polygon>"); cat("\n")
										cat("</kml>"); cat("\n")
										sink(NULL)				
									}
								p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
								pol = sps; proj4string(pol) = gadms2_USA@proj4string
								random_point = spsample(pol, 1, type="random")@coords
								tab2[i,"latitude"] = random_point[1,2]; tab2[i,"longitude"] = random_point[1,1]
							}
					}
				indices = which(tab2[,"trait"]%in%tip_labels[[h]]); tab2 = tab2[indices,]
				if (!file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.txt")))
					{
						write.table(tab2, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.txt"), row.names=F, quote=F, sep="\t")
					}
			}
		if (genotypes[h] == "D11")
			{
				tab1 = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.csv"), head=T)
				tab2 = tab1[,c("trait","location","latitude","longitude")]; tab2[,c("latitude","longitude")] = NA
				indices = which(tab2[,"trait"]%in%tip_labels[[h]]); tab2 = tab2[indices,]; colnames(tab2)[2] = "admin1"
				for (i in 1:dim(tab2)[1])
					{
						index1 = which(countries_admin1==unlist(strsplit(tab2[i,"admin1"],"-"))[1]); maxArea = 0; index3 = 0
						index2 = which(grepl(paste0("\\.",unlist(strsplit(tab2[i,"admin1"],"-"))[2]),gadms1[[index1]]@data[,"HASC_1"]))
						for (j in 1:length(gadms1[[index1]]@polygons[[index2]]@Polygons))
							{
								if (maxArea < gadms1[[index1]]@polygons[[index2]]@Polygons[[j]]@area)
									{
										maxArea = gadms1[[index1]]@polygons[[index2]]@Polygons[[j]]@area; index3 = j
									}
							}
						pol = gadms1[[index1]]@polygons[[index2]]@Polygons[[index3]]
						if (!file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons/",tab2[i,"admin1"],".kml")))
							{
								sink(file=paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons/",tab2[i,"admin1"],".kml"))
								cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"); cat("\n")
								cat("<kml xmlns=\"http://earth.google.com/kml/2.2\">"); cat("\n")
								cat(paste("\t<polygon id=\"",tab2[i,"admin1"],"\" samplingProbability=\"",1,"\">",sep="")); cat("\n")
								cat("\t\t<coordinates>"); cat("\n")
								for (j in 1:dim(pol@coords)[1])
									{
										cat(paste("\t\t\t",pol@coords[j,2],",",pol@coords[j,1],",0",sep="")); cat("\n")
									}
								cat("\t\t</coordinates>"); cat("\n")
								cat("\t</polygon>"); cat("\n")
								cat("</kml>"); cat("\n")
								sink(NULL)				
							}
						p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
						pol = sps; proj4string(pol) = gadms1[[index1]]@proj4string
						random_point = spsample(pol, 1, type="random")@coords
						tab2[i,"latitude"] = random_point[1,2]; tab2[i,"longitude"] = random_point[1,1]
					}
				if (!file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.txt")))
					{
						write.table(tab2, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.txt"), row.names=F, quote=F, sep="\t")
					}
			}
	}

	# 1.3. Editing the XML file with the polygons

		# 1.3.1. First approach: from an existing XML file
	
for (h in 1:length(genotypes))
	{
		if (file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_emp.xml")))
			{
				directory = paste0(genotypes[h],"_all_polygons")
				tab = read.table(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.txt"), head=T, sep="\t")
				xml = scan(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_emp.xml"), what="", sep="\n", quiet=T, blank.lines.skip=F)
				xml = gsub(paste0(genotypes[h],"_RRW_emp"), paste0(genotypes[h],"_RRW_pols"), xml)
				sink(file=paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_pols1.xml"))
				for (i in 1:length(xml))
					{
						cat(xml[i]); cat("\n")
						if (xml[i] == "\t<!-- INSERT leafTraitParameter -->")
							{
								cat("\n")
								for (j in 1:dim(tab)[1])
									{
										cat(paste("\t<leafTraitParameter id=\"",tab[j,"trait"],".trait\" taxon=\"",tab[j,"trait"],"\">",sep="")); cat("\n")
										cat(paste("\t\t<treeModel idref=\"treeModel\"/>",sep="")); cat("\n")
										cat(paste("\t\t<parameter idref=\"leaf.location\"/>",sep="")); cat("\n")
										cat(paste("\t</leafTraitParameter>",sep="")); cat("\n")
									}
								cat("\n")
								for (j in 1:dim(tab)[1])
									{
										cat(paste("\t<flatGeoSpatialPrior id=\"",tab[j,"trait"],"_polygon\" taxon=\"",tab[j,"trait"],"\" kmlFileName=\"",directory,"/",tab[j,"admin1"],".kml\" inside=\"true\" union=\"true\" cache=\"true\">",sep="")); cat("\n")
										cat(paste("\t\t<data>",sep="")); cat("\n")
										cat(paste("\t\t\t<parameter idref=\"",tab[j,"trait"],".trait\"/>",sep="")); cat("\n")
										cat(paste("\t\t</data>",sep="")); cat("\n")
										cat(paste("\t</flatGeoSpatialPrior>",sep="")); cat("\n")
									}
								cat("\n")		
							}
						if (xml[i]=="\t\t<!-- INSERT uniformGeoSpatialOperator -->")
							{
								cat("\n")
								for (j in 1:dim(tab)[1])
									{
										cat(paste("\t\t<uniformGeoSpatialOperator weight=\"0.01\">",sep="")); cat("\n")
										cat(paste("\t\t\t<parameter idref=\"",tab[j,"trait"],".trait\"/>",sep="")); cat("\n")
										cat(paste("\t\t\t<flatGeoSpatialPrior idref=\"",tab[j,"trait"],"_polygon\"/>",sep="")); cat("\n")
										cat(paste("\t\t</uniformGeoSpatialOperator>",sep="")); cat("\n")
									}
								cat("\n")
							}
						if (xml[i]=="\t\t\t\t<!-- INSERT geoDistributionCollection -->")
							{
								cat("\n")
								cat("\t\t\t\t<geoDistributionCollection id=\"allGeoDistributions\">"); cat("\n")
								for (j in 1:dim(tab)[1])
									{
										cat(paste("\t\t\t\t<flatGeoSpatialPrior idref=\"",tab[j,"trait"],"_polygon\"/>",sep="")); cat("\n")
									}
								cat("\t\t\t\t</geoDistributionCollection>"); cat("\n")
								cat("\n")
							}					
					}
				sink(NULL)
			}
	}

		# 1.3.2. Second approach: from a template XML file

for (h in 1:length(genotypes))
	{
		directory = paste0(genotypes[h],"_all_polygons")
		fasta = scan(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_skygrid.fasta"), what="", sep="\n", quiet=T)
		tab = read.table(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.txt"), head=T, sep="\t")
		xml = scan(paste0("Template.xml"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		xml = gsub(paste0("GENOTYPE"), paste0(genotypes[h]), xml)
		sink(file=paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_pols.xml"))
		for (i in 1:length(xml))
			{
				cat(xml[i]); cat("\n")
				if ((xml[i] == "\t\t</traitParameter>")&(genotypes[h] == "B32"))
					{
						cat("\t\t<jitter window=\"0.05 0.05\" duplicatesOnly=\"true\">"); cat("\n")
						cat("\t\t\t<parameter idref=\"leaf.location\"/>"); cat("\n")
						cat("\t\t</jitter>"); cat("\n")
					}
				if (xml[i] == "	<taxa id=\"taxa\">")
					{
						for (j in 1:dim(tab)[1])
							{
								year = decimal_date(ymd(unlist(strsplit(tab[j,"trait"],"\\|"))[length(unlist(strsplit(tab[j,"trait"],"\\|")))-2]))
								cat(paste("\t\t<taxon id=\"",tab[j,"trait"],"\">",sep="")); cat("\n")
								cat(paste("\t\t\t<date value=\"",year,"\" direction=\"forwards\" units=\"years\"/>",sep="")); cat("\n")
								cat(paste("\t\t\t<attr name=\"latitude\">",sep="")); cat("\n")
								cat(paste("\t\t\t\t",tab[j,"latitude"],sep="")); cat("\n")
								cat(paste("\t\t\t</attr>",sep="")); cat("\n")
								cat(paste("\t\t\t<attr name=\"longitude\">",sep="")); cat("\n")
								cat(paste("\t\t\t\t",tab[j,"longitude"],sep="")); cat("\n")
								cat(paste("\t\t\t</attr>",sep="")); cat("\n")
								cat(paste("\t\t\t",sep="")); cat("\n")
								cat(paste("\t\t\t<!-- START Multivariate diffusion model                                      -->",sep="")); cat("\n")
								cat(paste("\t\t\t<attr name=\"location\">",sep="")); cat("\n")
								cat(paste("\t\t\t\t",tab[j,"latitude"]," ",tab[j,"longitude"],sep="")); cat("\n")
								cat(paste("\t\t\t</attr>",sep="")); cat("\n")
								cat(paste("\t\t\t",sep="")); cat("\n")
								cat(paste("\t\t\t<!-- END Multivariate diffusion model                                        -->",sep="")); cat("\n")
								cat(paste("\t\t\t",sep="")); cat("\n")
								cat(paste("\t\t</taxon>",sep="")); cat("\n")								
							}
					}
				if (xml[i] == "	<alignment id=\"alignment\" dataType=\"nucleotide\">")
					{
						for (j in 1:dim(tab)[1])
							{
								sequence = fasta[which(fasta==paste0(">",tab[j,"trait"]))+1]
								cat(paste("\t\t<sequence>",sep="")); cat("\n")
								cat(paste("\t\t\t<taxon idref=\"",tab[j,"trait"],"\"/>",sep="")); cat("\n")
								cat(paste("\t\t\t",sequence,sep="")); cat("\n")
								cat(paste("\t\t</sequence>",sep="")); cat("\n")								
							}
					}				
				if (xml[i] == "\t<!-- INSERT leafTraitParameter -->")
					{
						cat("\n")
						for (j in 1:dim(tab)[1])
							{
								usePolygon = TRUE
								if (genotypes[h] == "B32")
									{
										if (grepl("precise",tab[j,"precision"])|grepl("http",tab[j,"precision"])|grepl("city",tab[j,"precision"])
											|grepl("beach",tab[j,"precision"])|grepl("district",tab[j,"precision"]))
											{
												usePolygon = FALSE
											}
									}
								if (usePolygon)
									{
										cat(paste("\t<leafTraitParameter id=\"",tab[j,"trait"],".trait\" taxon=\"",tab[j,"trait"],"\">",sep="")); cat("\n")
										cat(paste("\t\t<treeModel idref=\"treeModel\"/>",sep="")); cat("\n")
										cat(paste("\t\t<parameter idref=\"leaf.location\"/>",sep="")); cat("\n")
										cat(paste("\t</leafTraitParameter>",sep="")); cat("\n")
									}
							}
						cat("\n")
						for (j in 1:dim(tab)[1])
							{
								usePolygon = TRUE
								if (genotypes[h] == "B32")
									{
										if (grepl("precise",tab[j,"precision"])|grepl("http",tab[j,"precision"])|grepl("city",tab[j,"precision"])
											|grepl("beach",tab[j,"precision"])|grepl("district",tab[j,"precision"]))
											{
												usePolygon = FALSE
											}
									}
								if (usePolygon)
									{
										if (genotypes[h] != "D11")
											{
												if (is.na(tab[j,"admin1"])) fileName = tab[j,"country"]
												if ((!is.na(tab[j,"admin1"]))&(is.na(tab[j,"admin2"]))) fileName = paste0(unlist(strsplit(tab[j,"country"],"-"))[1],"_",gsub(" ","",tab[j,"admin1"]))
												if (!is.na(tab[j,"admin2"])) fileName = paste0(gsub(" ","",tab[j,"admin1"]),"_",gsub(" ","",tab[j,"admin2"]))
											}
										if (genotypes[h] == "D11")
											{
												fileName = tab[j,"admin1"]
											}
										cat(paste("\t<flatGeoSpatialPrior id=\"",tab[j,"trait"],"_polygon\" taxon=\"",tab[j,"trait"],"\" kmlFileName=\"",directory,"/",fileName,".kml\" inside=\"true\" union=\"true\" cache=\"true\">",sep="")); cat("\n")
										cat(paste("\t\t<data>",sep="")); cat("\n")
										cat(paste("\t\t\t<parameter idref=\"",tab[j,"trait"],".trait\"/>",sep="")); cat("\n")
										cat(paste("\t\t</data>",sep="")); cat("\n")
										cat(paste("\t</flatGeoSpatialPrior>",sep="")); cat("\n")
									}
							}
						cat("\n")		
					}
				if (xml[i]=="\t\t<!-- INSERT uniformGeoSpatialOperator -->")
					{
						cat("\n")
						for (j in 1:dim(tab)[1])
							{
								usePolygon = TRUE
								if (genotypes[h] == "B32")
									{
										if (grepl("precise",tab[j,"precision"])|grepl("http",tab[j,"precision"])|grepl("city",tab[j,"precision"])
											|grepl("beach",tab[j,"precision"])|grepl("district",tab[j,"precision"]))
											{
												usePolygon = FALSE
											}
									}
								if (usePolygon)
									{
										cat(paste("\t\t<uniformGeoSpatialOperator weight=\"0.01\">",sep="")); cat("\n")
										cat(paste("\t\t\t<parameter idref=\"",tab[j,"trait"],".trait\"/>",sep="")); cat("\n")
										cat(paste("\t\t\t<flatGeoSpatialPrior idref=\"",tab[j,"trait"],"_polygon\"/>",sep="")); cat("\n")
										cat(paste("\t\t</uniformGeoSpatialOperator>",sep="")); cat("\n")
									}
							}
						cat("\n")
					}
				if (xml[i]=="\t\t\t\t<!-- INSERT geoDistributionCollection -->")
					{
						cat("\n")
						cat("\t\t\t\t<geoDistributionCollection id=\"allGeoDistributions\">"); cat("\n")
						for (j in 1:dim(tab)[1])
							{
								usePolygon = TRUE
								if (genotypes[h] == "B32")
									{
										if (grepl("precise",tab[j,"precision"])|grepl("http",tab[j,"precision"])|grepl("city",tab[j,"precision"])
											|grepl("beach",tab[j,"precision"])|grepl("district",tab[j,"precision"]))
											{
												usePolygon = FALSE
											}
									}
								if (usePolygon)
									{
										cat(paste("\t\t\t\t<flatGeoSpatialPrior idref=\"",tab[j,"trait"],"_polygon\"/>",sep="")); cat("\n")
									}
							}
						cat("\t\t\t\t</geoDistributionCollection>"); cat("\n")
						cat("\n")
					}					
			}
		sink(NULL)
	}

# 2. Extracting the spatio-temporal information embedded in postrtior trees

for (h in 1:length(genotypes))
	{
		trees = scan(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_pols.trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		indices3 = which(grepl("tree STATE",trees)); index3 = indices3[length(indices3)]
		burnIn = (0.1*(length(indices3)-1))+1
		if (length(indices3) == 1101) burnIn = 101
		index1 = which(trees=="\t\t;")[length(which(trees=="\t\t;"))]; index2 = index1 + burnIn + 1
		nberOfTreesToSample = nberOfExtractionFiles; interval = floor((index3-(index1+burnIn))/nberOfTreesToSample)
		indices = seq(index3-((nberOfTreesToSample-1)*interval),index3,interval)
		selected_trees = c(trees[c(1:index1,indices)],"End;")
		write(selected_trees, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000.trees"))	
		file1 = paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000.trees")
		file2 = paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000.tree")
		trees = readAnnotatedNexus(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000.trees"))
		system(paste0("BEAST_1_10_5/bin/treeannotator -burninTrees 0 -heights keep -hpd2D 0.8 ",file1," ",file2), ignore.stdout=F, ignore.stderr=F)
		tab = read.table(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.txt"), head=T, sep="\t"); collection_dates = rep(NA, dim(tab)[1])
		for (i in 1:dim(tab)[1]) collection_dates[i] = unlist(strsplit(tab[i,"trait"],"\\|"))[length(unlist(strsplit(tab[i,"trait"],"\\|")))-2]
		mostRecentSamplingDatum = max(decimal_date(ymd(collection_dates)))
		mcc_tre = readAnnotatedNexus(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000.tree"))
		mcc = mccTreeExtractions(mcc_tre, mostRecentSamplingDatum)
		write.csv(mcc, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000.csv"), row.names=F, quote=F)
		dir.create(file.path(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000_ext")), showWarnings=F)
		trees = readAnnotatedNexus(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000.trees"))
		for (i in 1:length(trees))
			{
				tab = postTreeExtractions(trees[[i]], mostRecentSamplingDatum)
				write.csv(tab, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000_ext/TreeExtractions_",i,".csv"), row.names=F, quote=F)
			}
	}
for (h in 1:length(genotypes))
	{
		if (h == 1)
			{
				mccs = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000.csv"), head=T, sep=",")
				mccs = cbind(mccs, rep(genotypes[h],dim(mccs)[1])); colnames(mccs)[length(colnames(mccs))] = "genotype"
			}	else	{
				mcc = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000.csv"), head=T, sep=",")
				mcc = cbind(mcc, rep(genotypes[h],dim(mcc)[1])); colnames(mcc)[length(colnames(mcc))] = "genotype"
				mccs = rbind(mccs, mcc)
			}
	}
mccs = mccs[order(mccs[,"endYear"],mccs[,"startYear"],decreasing=F),]
write.csv(mccs, "ALL_RRW.csv", row.names=F, quote=F)
dir.create(file.path(paste0("Genotype_ALL")), showWarnings=F)
for (i in 1:nberOfExtractionFiles)
	{
		for (h in 1:length(genotypes))
			{
				if (h == 1)
					{
						tabs = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000_ext/TreeExtractions_",i,".csv"), head=T, sep=",")
						tabs = cbind(tabs, rep(genotypes[h],dim(tabs)[1])); colnames(tabs)[length(colnames(tabs))] = "genotype"
					}	else	{
						tab = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000_ext/TreeExtractions_",i,".csv"), head=T, sep=",")
						tab = cbind(tab, rep(genotypes[h],dim(tab)[1])); colnames(tab)[length(colnames(tab))] = "genotype"
						tabs = rbind(tabs, tab)
					}
				}
		tabs = tabs[order(tabs[,"endYear"],tabs[,"startYear"],decreasing=F),]
		write.csv(tabs, paste0("Genotype_ALL/TreeExtractions_",i,".csv"), row.names=F, quote=F)
	}

# 3. Estimating lineage dispersal statistics from the continuous phylogeographic analyses

timeSlices = 100; onlyTipBranches = F; showingPlots = F; nberOfCores = 10; slidingWindow = 1/24
for (h in 1:length(genotypes))
	{
		localTreesDirectory = paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000_ext")
		outputName = paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_statistics/",genotypes[h])
		dir.create(file.path(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_statistics/")), showWarnings=F)
		spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)
	 }

# 4. Visualising the continuous phylogeographic reconstructions based on the polygons

e1 = extent(-175, -25, -80, 77); e2 = extent(-100, -60, 23, 50); e3 = extent(-165, -49, 13, 77)
countries_1 = crop(shapefile("All_shapefiles/World_borders/World_countries_shapefile.shp"), e1)
borders_1 = crop(shapefile("All_shapefiles/Natural_Earth/International_borders/Only_international_borders.shp"), e1)
coasts_1 = crop(shapefile("All_shapefiles/Natural_Earth/Coast_lines_borders/Only_coast_lines_borders.shp"), e1)
countries_2 = crop(shapefile("All_shapefiles/World_borders/World_countries_shapefile.shp"), e2)
borders_2 = crop(shapefile("All_shapefiles/Natural_Earth/International_borders/Only_international_borders.shp"), e2)
coasts_2 = crop(shapefile("All_shapefiles/Natural_Earth/Coast_lines_borders/Only_coast_lines_borders.shp"), e2)
countries_3 = crop(shapefile("All_shapefiles/World_borders/World_countries_shapefile.shp"), e3)
borders_3 = crop(shapefile("All_shapefiles/Natural_Earth/International_borders/Only_international_borders.shp"), e3)
coasts_3 = crop(shapefile("All_shapefiles/Natural_Earth/Coast_lines_borders/Only_coast_lines_borders.shp"), e3)
colour_scale = met.brewer(name="Hiroshige", n=111, type="continuous")[1:101]

	# 4.1. Individual figures maps per each genotype

for (h in 1:length(genotypes))
	{
		mcc = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000.csv"), head=T, sep=",")
		minYear = min(mcc[,"startYear"]); maxYear = max(mcc[,"endYear"])
		endYears_indices = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
		endYears_colours = colour_scale[endYears_indices]
		localTreesDirectory = paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000_ext")
		prob = 0.80; precision = 1/12; startDatum = minYear
		polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))
		polygons_colours = rep(NA, length(polygons))
		for (i in 1:length(polygons))
			{
				date = as.numeric(names(polygons[[i]]))
				polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
				polygons_colours[i] = paste0(colour_scale[polygon_index],"20")
			}
		cutOffs = c(maxYear); g = 1; croppingPolygons = FALSE; plottingAllNodes = FALSE
		if (genotypes[h] == "B32")
			{
				pdf(paste0(genotypes[h],"_RRW.pdf"), width=5.5, height=5.7)
				par(oma=c(0,0,0,0), mar=c(1,0.4,1,0.4), mgp=c(0,0.1,0), lwd=0.3, bty="o")
				plot(countries_1, col="gray90", border=NA, ann=F, axes=F)
				plot(borders_1, col="white", lwd=0.3, add=T)
				plot(coasts_1, col="gray70", lwd=0.5, add=T)
				rast = raster(matrix(nrow=1, ncol=2))
				rast[1] = minYear; rast[2] = maxYear
				if (g == length(cutOffs))
					{
						plot(rast, legend.only=T, add=T, col=colour_scale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.15,0.50,0.280,0.287),
							 legend.args=list(text="", cex=0.7, col="gray30"), horizontal=T,
							 axis.args=list(cex.axis=0.5, lwd=0, lwd.tick=0.3, col.tick="gray30", tck=-1, col="gray30", col.axis="gray30", line=0, mgp=c(0,-0.15,0)))
		 			}
				for (i in 1:length(polygons))
					{
						if (as.numeric(names(polygons[[i]])) <= cutOffs[g]) plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
					}
				for (i in 1:dim(mcc)[1])
					{
						if (mcc[i,"endYear"] <= cutOffs[g])
							{
								curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
						  		  		    arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
						  	}
					}
				for (i in dim(mcc)[1]:1)
					{
						if ((mcc[i,"startYear"] <= cutOffs[g])&(!mcc[i,"node1"]%in%mcc[,"node2"]))
							{
								startYears_index = (((mcc[i,"startYear"]-minYear)/(maxYear-minYear))*100)+1
								points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=colour_scale[startYears_index], cex=0.4)
								points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", lwd=0.2, cex=0.4)
							}
						if (mcc[i,"endYear"] <= cutOffs[g])
							{
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.4)
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", lwd=0.2, cex=0.4)
							}
					}
				rect(e1@xmin, e1@ymin, e1@xmax, e1@ymax, lwd=0.3, border="gray30")
				dev.off()
			}
		if (genotypes[h] == "C21")
			{
				pdf(paste0(genotypes[h],"_RRW.pdf"), width=5.5, height=4.7)
				par(oma=c(0,0,0,0), mar=c(1,0.4,1,0.4), mgp=c(0,0.1,0), lwd=0.3, bty="o")
				plot(countries_2, col="gray90", border=NA, ann=F, axes=F)
				plot(borders_2, col="white", lwd=0.4, add=T)
				plot(coasts_2, col="gray70", lwd=0.6, add=T)
				rast = raster(matrix(nrow=1, ncol=2))
				rast[1] = minYear; rast[2] = maxYear
				if (g == length(cutOffs))
					{
						plot(rast, legend.only=T, add=T, col=colour_scale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.57,0.88,0.290,0.298),
							 legend.args=list(text="", cex=0.7, col="gray30"), horizontal=T,
							 axis.args=list(cex.axis=0.5, lwd=0, lwd.tick=0.3, col.tick="gray30", tck=-0.9, col="gray30", col.axis="gray30", line=0, mgp=c(0,-0.2,0)))
		 			}
				for (i in 1:length(polygons))
					{
						if (as.numeric(names(polygons[[i]])) <= cutOffs[g]) plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
					}
				for (i in 1:dim(mcc)[1])
					{
						if (mcc[i,"endYear"] <= cutOffs[g])
							{
								curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
						  		  		    arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
						  	}
					}
				for (i in dim(mcc)[1]:1)
					{
						if ((mcc[i,"startYear"] <= cutOffs[g])&(!mcc[i,"node1"]%in%mcc[,"node2"]))
							{
								startYears_index = (((mcc[i,"startYear"]-minYear)/(maxYear-minYear))*100)+1
								points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=colour_scale[startYears_index], cex=0.55)
								points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", lwd=0.2, cex=0.55)
							}
						if (mcc[i,"endYear"] <= cutOffs[g])
							{
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.55)
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", lwd=0.2, cex=0.55)
							}
					}
				rect(e2@xmin, e2@ymin, e2@xmax, e2@ymax, lwd=0.3, border="gray30")
				dev.off()
			}
		if ((genotypes[h] == "B36")|(genotypes[h] == "D11"))
			{	
				pdf(paste0(genotypes[h],"_RRW.pdf"), width=5.3, height=4.2)
				par(oma=c(0,0,0,0), mar=c(0.0,0.4,0.2,0.4), mgp=c(0,0.1,0), lwd=0.3, bty="o")
				plot(countries_3, col="gray90", border=NA, ann=F, axes=F)
				plot(borders_3, col="white", lwd=0.40, add=T)
				plot(coasts_3, col="gray70", lwd=0.40, add=T)
				rast = raster(matrix(nrow=1, ncol=2))
				rast[1] = minYear; rast[2] = maxYear
				if (g == length(cutOffs))
					{
						plot(rast, legend.only=T, add=T, col=colour_scale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.08,0.40,0.108,0.117),
							 legend.args=list(text="", cex=0.7, col="gray30"), horizontal=T,
							 axis.args=list(cex.axis=0.45, lwd=0, lwd.tick=0.3, col.tick="gray30", tck=-1.0, col="gray30", col.axis="gray30", line=0, mgp=c(0,-0.20,0)))
		 			}
				for (i in 1:length(polygons))
					{
						if (as.numeric(names(polygons[[i]])) <= cutOffs[g]) plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
					}
				for (i in 1:dim(mcc)[1])
					{
						if (mcc[i,"endYear"] <= cutOffs[g])
							{
								curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
						  		  		    arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
						  	}
					}
				for (i in dim(mcc)[1]:1)
					{
						if ((mcc[i,"startYear"] <= cutOffs[g])&(!mcc[i,"node1"]%in%mcc[,"node2"]))
							{
								startYears_index = (((mcc[i,"startYear"]-minYear)/(maxYear-minYear))*100)+1
								points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=colour_scale[startYears_index], cex=0.4)
								points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", lwd=0.2, cex=0.4)
							}
						if (mcc[i,"endYear"] <= cutOffs[g])
							{
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.4)
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", lwd=0.2, cex=0.4)
							}
					}
				rect(e3@xmin, e3@ymin, e3@xmax, e3@ymax, lwd=0.3, border="gray30")
				dev.off()
			}
	}

	# 4.2. Composite figure with one map/genotype

minYear = 9999; maxYear = -9999
for (h in 1:length(genotypes))
	{
		mcc = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000.csv"), head=T, sep=",")
		if (minYear > min(mcc[,"startYear"])) minYear = min(mcc[,"startYear"])
		if (maxYear < max(mcc[,"endYear"])) maxYear = max(mcc[,"endYear"])	
	}
endYears_colours_list = list(); polygons_list = list(); polygons_colours_list = list()
for (h in 1:length(genotypes))
	{	
		mcc = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000.csv"), head=T, sep=",")
		endYears_indices = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
		endYears_colours_list[[h]] = colour_scale[endYears_indices]
		localTreesDirectory = paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000_ext")
		prob = 0.80; precision = 1/12; startDatum = minYear
		polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))
		polygons_colours = rep(NA, length(polygons))
		for (i in 1:length(polygons))
			{
				date = as.numeric(names(polygons[[i]]))
				polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
				polygons_colours[i] = paste0(colour_scale[polygon_index],"20")
			}
		polygons_list[[h]] = polygons; polygons_colours_list[[h]] = polygons_colours
	}
rast = raster(matrix(nrow=1, ncol=2)); rast[1] = minYear; rast[2] = maxYear
for (h in 1:length(genotypes))
	{
		mcc = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000.csv"), head=T, sep=",")
		pdf(paste0("All1_RRW_",genotypes[h],".pdf"), width=5.5, height=5.7)
		par(oma=c(0,0,0,0), mar=c(1,0.4,1,0.4), mgp=c(0,0.1,0), lwd=0.3, bty="o")
		plot(countries_1, col="gray90", border=NA, ann=F, axes=F)
		plot(borders_1, col="white", lwd=0.3, add=T)
		plot(coasts_1, col="gray70", lwd=0.5, add=T)
		plot(rast, legend.only=T, add=T, col=colour_scale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.15,0.50,0.280,0.287),
			 legend.args=list(text="", cex=0.7, col="gray30"), horizontal=T,
			 axis.args=list(cex.axis=0.5, lwd=0, lwd.tick=0.3, col.tick="gray30", tck=-1, col="gray30", col.axis="gray30", line=0, mgp=c(0,-0.15,0)))
		for (i in 1:length(polygons))
			{
				plot(polygons_list[[h]][[i]], axes=F, col=polygons_colours_list[[h]][i], add=T, border=NA)
			}
		for (i in 1:dim(mcc)[1])
			{
				curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
							arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
			}
		for (i in dim(mcc)[1]:1)
			{
				if (!mcc[i,"node1"]%in%mcc[,"node2"])
					{
						startYears_index = (((mcc[i,"startYear"]-minYear)/(maxYear-minYear))*100)+1
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=colour_scale[startYears_index], cex=0.4)
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", lwd=0.2, cex=0.4)
					}
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours_list[[h]][i], cex=0.4)
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", lwd=0.2, cex=0.4)
			}
		rect(e1@xmin, e1@ymin, e1@xmax, e1@ymax, lwd=0.3, border="gray30")
		dev.off()
	}
polygon_list = list(); prob = 0.80
for (h in 1:length(genotypes))
	{
		coords = c()
		localTreesDirectory = paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000_ext")
		for (i in 1:nberOfExtractionFiles)
			{
				tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), head=T)
				coords = rbind(coords, tab[which(!tab["node1"]%in%tab[,"node2"])[1],c("endLon","endLat")])
			}
		kde = kde(coords, H=Hpi(coords), compute.cont=T, gridsize=c(1000,1000))
		contourLevel = contourLevels(kde, prob=(1-prob)); polygons = list()
		contourLines = contourLines(kde$eval.points[[1]], kde$eval.points[[2]], kde$estimate, level=contourLevel)
		for (i in 1:length(contourLines)) polygons[[i]] = Polygon(cbind(contourLines[[i]]$x,contourLines[[i]]$y))
		ps = Polygons(polygons,1); contourPolygons = SpatialPolygons(list(ps))
		spdf = SpatialPolygonsDataFrame(contourPolygons, data.frame(ID=1:length(contourPolygons)))
		names(spdf) = genotypes[h]; polygon_list[[h]] = spdf
	}
cols1 = list(); cols2 = list()
cols1[[1]] = rgb(222,67,39,220,maxColorValue=255); cols2[[1]] = rgb(222,67,39,100,maxColorValue=255) # red
cols1[[2]] = rgb(250,165,33,220,maxColorValue=255); cols2[[2]] = rgb(250,165,33,100,maxColorValue=255) # orange
cols1[[3]] = rgb(106,61,154,255,maxColorValue=255); cols2[[3]] = rgb(106,61,154,100,maxColorValue=255) # purple
cols1[[4]] = rgb(70,118,187,220,maxColorValue=255); cols2[[4]] = rgb(70,118,187,100,maxColorValue=255) # blue
pdf(paste0("Ancestors_NEW.pdf"), width=5.5, height=5.7)
par(oma=c(0,0,0,0), mar=c(1,0.4,1,0.4), mgp=c(0,0.1,0), lwd=0.3, bty="o")
plot(countries_1, col="gray90", border=NA, ann=F, axes=F)
plot(borders_1, col="white", lwd=0.3, add=T)
plot(coasts_1, col="gray70", lwd=0.5, add=T)
for (h in 1:length(genotypes))
	{
		plot(polygon_list[[h]], axes=F, col=cols2[[h]], border=NA, add=T)
	}
for (h in 1:length(genotypes))
	{
		plot(polygon_list[[h]], axes=F, col=NA, border=cols1[[h]], lwd=0.5, add=T)
	}
rect(e1@xmin, e1@ymin, e1@xmax, e1@ymax, lwd=0.3, border="gray30")
dev.off()

	# 4.3. Unique map with all the dispersal histories

mcc = read.csv("ALL_RRW.csv", head=T); localTreesDirectory = "Genotype_ALL"
minYear = min(mcc[,"startYear"]); maxYear = max(mcc[,"endYear"])
endYears_indices = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
endYears_colours = colour_scale[endYears_indices]
prob = 0.80; precision = 1/12; startDatum = minYear
polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))
polygons_colours = rep(NA, length(polygons))
for (i in 1:length(polygons))
	{
		date = as.numeric(names(polygons[[i]]))
		polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
		polygons_colours[i] = paste0(colour_scale[polygon_index],"20")
	}
rast = raster(matrix(nrow=1, ncol=2)); rast[1] = minYear; rast[2] = maxYear
pdf(paste0("ALL_RRW_NEW.pdf"), width=11.0, height=11.4); cutOffs = c(2022.5,2023,2023.5,2025)
par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(1,0.4,1,0.4), mgp=c(0,0.1,0), lwd=0.3, bty="o")
for (g in 1:length(cutOffs)) # date_decimal(c(2022.5,2023,2023.5,2025))
	{
		plot(countries_1, col="gray90", border=NA, ann=F, axes=F)
		plot(borders_1, col="white", lwd=0.3, add=T)
		plot(coasts_1, col="gray70", lwd=0.5, add=T)
		plot(rast, legend.only=T, add=T, col=colour_scale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.15,0.50,0.280,0.287),
			 legend.args=list(text="", cex=0.7, col="gray30"), horizontal=T,
			 axis.args=list(cex.axis=0.5, lwd=0, lwd.tick=0.3, col.tick="gray30", tck=-1, col="gray30", col.axis="gray30", line=0, mgp=c(0,-0.15,0)))
		for (i in 1:length(polygons))
			{
				if (as.numeric(names(polygons[[i]])) <= cutOffs[g]) plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
			}
		for (i in 1:dim(mcc)[1])
			{
				if (mcc[i,"endYear"] <= cutOffs[g])
					{
						curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
									arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
					}
			}
		for (i in 1:dim(mcc)[1])
			{
				if ((mcc[i,"startYear"] <= cutOffs[g])&(!mcc[i,"node1"]%in%mcc[,"node2"]))
					{
						startYears_index = (((mcc[i,"startYear"]-minYear)/(maxYear-minYear))*100)+1
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=colour_scale[startYears_index], cex=0.4)
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", lwd=0.2, cex=0.4)
					}
				if (mcc[i,"endYear"] <= cutOffs[g])
					{
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.4)
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", lwd=0.2, cex=0.4)
					}
			}
		rect(e1@xmin, e1@ymin, e1@xmax, e1@ymax, lwd=0.3, border="gray30")
	}
dev.off()

system(paste0("magick -units PixelsPerInch -density 1000 All1_RRW.pdf -background white -alpha remove -flatten All1_RRW.png"))
system(paste0("magick -units PixelsPerInch -density 1000 All2_RRW.pdf -background white -alpha remove -flatten All2_RRW.png"))

