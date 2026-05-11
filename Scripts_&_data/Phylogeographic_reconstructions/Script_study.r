library(diagram)
library(ks)
library(lubridate)
library(MetBrewer)
library(seraphim)
library(tibble)
library(tidytree)
library(treeio)

genotypes = c("A1","A3")
genotypes = c("AC","MM") 
genotypes = c("B32","B36","C21","D11")
nberOfExtractionFiles = 1000

	# Notes:
		# - some extra spaces had to be removed next to the country names in the metadata files exported in ".csv" files
		# - for A1, B3.2 and B3.6, "Ross's_goose" has to changed to "Ross_s_goose" in both the ".csv" and ".trees" files
		# - "AC" and "MM" refer to the Argentinian avian and marine mammal clades of the B3.2 genotype, respectively

# 1. Subsampling the sequences by phylogenetic clustering and then by month and state

	# 1.1. Subsampling by phylogenetic clustering

for (h in 1:length(genotypes))
	{
		if (!file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_2.fas")))
			{
				# tree = readAnnotatedNexus(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_2.tre")) # does not work...
				tree = treeio::read.beast(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_2.tre"))
				tab1 = as.data.frame(as_tibble(tree)); tab1 = tab1[-which(tab1[,"parent"]==tab1[,"node"]),]
				colnames(tab1)[which(colnames(tab1)=="label.x")] = "label"; colnames(tab1)[which(colnames(tab1)=="label.y")] = "bs"
				threshold = 0.7; branches = as.phylo(tree)$edge; bootstraps = rep(NA, dim(branches)[1])
				for (i in 1:length(bootstraps)) bootstraps[i]= tab1[which((tab1[,"parent"]==branches[i,1])&(tab1[,"node"]==branches[i,2])),"bs"]
				tab2 = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_1.csv"), head=T, sep=";")
				tab2[which(tab2[,"admin2"]==""),"admin2"] = "Unknown"; tab2$location = rep(NA, dim(tab2)[1])
				for (i in 1:dim(tab2)[1])
					{
						tab2[i,"location"] = paste0(tab2[i,"country"],"_",tab2[i,"admin1"],"_",tab2[i,"admin2"])
					}
				tree = as.phylo(tree); subTrees = subtrees(tree, wait=F); tipLabels = tree$tip.label
				considering_bootstrap_supports = TRUE # considering_bootstrap_supports = FALSE
				c1 = 0; c2 = 0; clusters1 = list(); clusters2 = list(); bootstraps1 = list(); bootstraps2 = list()
				for (i in 2:length(subTrees)) # the 1st subtree is the entire one
					{
						subTree = subTrees[i][[1]]; subEdges = tree$edge[which(tree$edge[,1]%in%subTree$node.label),]
						if (considering_bootstrap_supports)
							{
								root = subEdges[which(!subEdges[,1]%in%subEdges[,2])[1],1]
								bs1 = tab1[which(tab1[,"node"]==root),"bs"]; bs2 = bootstraps[which(tree$edge[,2]==root)]; 
								if ((!is.na(bs1))&&(!is.na(bs1))&&(bs1!=bs2)) print(cbind(bs1,bs2))
								bootstrap = bs1; labels = subTree$tip.label; states = rep(NA, length(labels))
								for (j in 1:length(labels)) states[j] = tab2[which(tab2[,"label"]==labels[j]),"location"]
								if ((!is.na(bootstrap))&&(bootstrap >= threshold)&&(length(unique(states)) == 1))
									{
										c1 = c1+1; clusters1[[c1]] = labels; bootstraps1[[c1]] = bootstrap
									}
							}	else	{
								labels = subTree$tip.label; states = rep(NA, length(labels))
								for (j in 1:length(labels)) states[j] = tab2[which(tab2[,"label"]==labels[j]),"location"]
								if (length(unique(states)) == 1)
									{
										c1 = c1+1; clusters1[[c1]] = labels
									}
							}
					}
				for (i in 1:length(clusters1))
					{
						nested = FALSE
						for (j in 1:length(clusters1))
							{
								if (i != j)
									{
										if (sum(clusters1[[i]]%in%clusters1[[j]]) == length(clusters1[[i]])) nested = TRUE
									}
							}
						if (nested == FALSE)
							{
								c2 = c2 + 1; clusters2[[c2]] = clusters1[[i]]
								if (considering_bootstrap_supports) bootstraps2[[c2]] = bootstraps1[[i]]
							}	
					}
				sequencesToRemove2 = c()
				for (i in 1:length(clusters2))
					{
						sequencesToRemove2 = c(sequencesToRemove2, sample(clusters2[[i]],length(clusters2[[i]])-1,replace=F))
					}
				tab2 = tab2[which(!tab2[,"label"]%in%sequencesToRemove2),]
				tab2 = tab2[which((tab2[,"admin1"]!="Unknown")&(!is.na(tab2[,"admin1"]))),]
				fasta1 = scan(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_1.fas"), what="", sep="\n", quiet=T)
				sequences1 = fasta1[which(!grepl(">",fasta1))]; ids1 = fasta1[which(grepl(">",fasta1))]
				sequences2 = sequences1[!gsub(">","",ids1)%in%sequencesToRemove2]
				ids2 = ids1[!gsub(">","",ids1)%in%sequencesToRemove2]; fasta2 = c()
				for (i in 1:length(ids2)) fasta2 = c(fasta2, ids2[i], sequences2[i])
				if (considering_bootstrap_supports)
					{
						write.csv(tab2, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_2.csv"), row.names=F, quote=F)
						write(fasta2, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_2.fas"))
					}	else	{
						write.csv(tab2, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_2_without_bootstraps.csv"), row.names=F, quote=F)
						write(fasta2, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_2_without_bootstraps.fas"))
					}
			}
		if ((genotypes[h] == "D11")&(!file.exists(paste0("Genotype_",genotypes[h],"/D11_DTA_samples/",genotypes[h],"_alignment_2.fas"))))
			{
				tree = treeio::read.beast(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_2.tre"))
				tab1 = as.data.frame(as_tibble(tree)); tab1 = tab1[-which(tab1[,"parent"]==tab1[,"node"]),]
				colnames(tab1)[which(colnames(tab1)=="label.x")] = "label"; colnames(tab1)[which(colnames(tab1)=="label.y")] = "bs"
				threshold = 0.7; branches = as.phylo(tree)$edge; bootstraps = rep(NA, dim(branches)[1])
				for (i in 1:length(bootstraps)) bootstraps[i]= tab1[which((tab1[,"parent"]==branches[i,1])&(tab1[,"node"]==branches[i,2])),"bs"]
				tab2 = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_1.csv"), head=T, sep=";")
				tab2$location = tab2$admin1; tree = as.phylo(tree); subTrees = subtrees(tree, wait=F); tipLabels = tree$tip.label
				considering_bootstrap_supports = TRUE # considering_bootstrap_supports = FALSE
				c1 = 0; c2 = 0; clusters1 = list(); clusters2 = list(); bootstraps1 = list(); bootstraps2 = list()
				for (i in 2:length(subTrees)) # the 1st subtree is the entire one
					{
						subTree = subTrees[i][[1]]; subEdges = tree$edge[which(tree$edge[,1]%in%subTree$node.label),]
						if (considering_bootstrap_supports)
							{
								root = subEdges[which(!subEdges[,1]%in%subEdges[,2])[1],1]
								bs1 = tab1[which(tab1[,"node"]==root),"bs"]; bs2 = bootstraps[which(tree$edge[,2]==root)]; 
								if ((!is.na(bs1))&&(!is.na(bs1))&&(bs1!=bs2)) print(cbind(bs1,bs2))
								bootstrap = bs1; labels = subTree$tip.label; states = rep(NA, length(labels))
								for (j in 1:length(labels)) states[j] = tab2[which(tab2[,"label"]==labels[j]),"location"]
								if ((!is.na(bootstrap))&&(bootstrap >= threshold)&&(length(unique(states)) == 1))
									{
										c1 = c1+1; clusters1[[c1]] = labels; bootstraps1[[c1]] = bootstrap
									}
							}	else	{
								labels = subTree$tip.label; states = rep(NA, length(labels))
								for (j in 1:length(labels)) states[j] = tab2[which(tab2[,"label"]==labels[j]),"location"]
								if (length(unique(states)) == 1)
									{
										c1 = c1+1; clusters1[[c1]] = labels
									}
							}
					}
				for (i in 1:length(clusters1))
					{
						nested = FALSE
						for (j in 1:length(clusters1))
							{
								if (i != j)
									{
										if (sum(clusters1[[i]]%in%clusters1[[j]]) == length(clusters1[[i]])) nested = TRUE
									}
							}
						if (nested == FALSE)
							{
								c2 = c2 + 1; clusters2[[c2]] = clusters1[[i]]
								if (considering_bootstrap_supports) bootstraps2[[c2]] = bootstraps1[[i]]
							}	
					}
				sequencesToRemove2 = c()
				for (i in 1:length(clusters2))
					{
						sequencesToRemove2 = c(sequencesToRemove2, sample(clusters2[[i]],length(clusters2[[i]])-1,replace=F))
					}
				tab2 = tab2[which(!tab2[,"label"]%in%sequencesToRemove2),]
				tab2 = tab2[which((tab2[,"admin1"]!="")&(tab2[,"admin1"]!="Unknown")&(!is.na(tab2[,"admin1"])),] # (tab2[,"admin1"]!="") was missing...!
				fasta1 = scan(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_1.fas"), what="", sep="\n", quiet=T)
				sequences1 = fasta1[which(!grepl(">",fasta1))]; ids1 = fasta1[which(grepl(">",fasta1))]
				sequences2 = sequences1[!gsub(">","",ids1)%in%sequencesToRemove2]
				ids2 = ids1[!gsub(">","",ids1)%in%sequencesToRemove2]; fasta2 = c()
				for (i in 1:length(ids2)) fasta2 = c(fasta2, ids2[i], sequences2[i])
				if (considering_bootstrap_supports)
					{
						write.csv(tab2, paste0("Genotype_",genotypes[h],"/D11_DTA_samples/",genotypes[h],"_alignment_2.csv"), row.names=F, quote=F)
						write(fasta2, paste0("Genotype_",genotypes[h],"/D11_DTA_samples/",genotypes[h],"_alignment_2.fas"))
					}	else	{
						write.csv(tab2, paste0("Genotype_",genotypes[h],"/D11_DTA_samples/",genotypes[h],"_alignment_2_without_bootstraps.csv"), row.names=F, quote=F)
						write(fasta2, paste0("Genotype_",genotypes[h],"/D11_DTA_samples/",genotypes[h],"_alignment_2_without_bootstraps.fas"))
					}
			}
	}

	# 1.2. Subsampling by month and by admin-1

n = 2 # sampling two sequences (when available) per combination of month and admin-1 location
for (h in c(1,4)) # this second subsampling step was only conducted for B3.2 and D1.1
	{
		if (!file.exist(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_3.fas")))
			{
				fasta2 = scan(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_2.fas"), what="", sep="\n", quiet=T)
				tab2 = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_2.csv"), head=T, sep=",")
				locations = paste0(tab2[,"country"],"_",tab2[,"admin1"]); dates = tab2[,"collection_date"]; months = rep(NA, length(dates))
				for (i in 1:length(months)) months[i] = paste(unlist(strsplit(dates[i],"-"))[1],unlist(strsplit(dates[i],"-"))[2],sep="-")
				different_locations = unique(locations); different_locations = different_locations[order(different_locations)]
				different_months = unique(months); different_months = different_months[order(different_months)]
				categories = paste0(months,"$",locations); different_categories = unique(categories); selected_sequences = c()
				for (i in 1:length(different_categories))
					{
						month = unlist(strsplit(different_categories[i],"\\$"))[1]; location = unlist(strsplit(different_categories[i],"\\$"))[2]
						indices = which((grepl(month,tab2[,"collection_date"]))&(paste0(tab2[,"country"],"_",tab2[,"admin1"])==location))
						sub = tab2[indices,]; sub[which(sub[,"admin2"]==""),"admin2"] = "Unknown"
						if (dim(sub)[1] <= n)
							{
								# print(c(i,"a"))
								selected_sequences = c(selected_sequences,sub[,"label"])					
							}	else	{
								n2 = n; sub2 = sub
								if (sum(sub[,"location_used"]!="") > 0)
									{
										n1 = sum(sub[,"location_used"]!=""); sub1 = sub[which(sub[,"location_used"]!=""),]
										n2 = n-n1; sub2 = sub[which(sub[,"location_used"]==""),]
										if (n1 <= n)
											{
												# print(c(i,"b"))
												selected_sequences = c(selected_sequences,sub1[sample(1:dim(sub1)[1],n1),"label"])
											}
										if (n1 > n)
											{
												coordinates = paste0(sub1[,"longitude"],"_",sub1[,"latitude"])
												unique_coordinates = unique(coordinates)
												if (length(unique_coordinates) < n)
													{
														# print(c(i,"c"))
														selected_sequences = c(selected_sequences, sub1[sample(1:dim(sub1)[1],n,replace=F),"label"])
													}	else	{
														# print(c(i,"d"))
														selected_coordinates = unique_coordinates[sample(1:length(unique_coordinates), n, replace=F)]
														for (j in 1:n)
															{
																coordinates_indices = which(coordinates==selected_coordinates[j])
																coordinates_index_selected = coordinates_indices[sample(1:length(coordinates_indices),1)]
																selected_sequences = c(selected_sequences, sub1[coordinates_index_selected,"label"])
															}
													}
											}
									}
								if (n2 > 0)
									{
										if (sum(sub2[,"admin2"]!="Unknown") == 0)
											{
												# print(c(i,"e"))
												selected_sequences = c(selected_sequences, sub2[sample(1:dim(sub2)[1],n2,replace=F),"label"])
											}
										if ((sum(sub2[,"admin2"]!="Unknown") > 0)&(sum(sub2[,"admin2"]!="Unknown") < n2))
											{
												# print(c(i,"f"))
												selected_sequences = c(selected_sequences, sub2[which(sub2[,"admin2"]!="Unknown"),"label"])
												possible_indices = which(sub2[,"admin2"]=="Unknown"); nn = n2-sum(sub2[,"admin2"]!="Unknown")
												remaining_indices = possible_indices[sample(1:length(possible_indices),nn,replace=F)]										
												selected_sequences = c(selected_sequences, sub2[remaining_indices,"label"])
											}
										if (sum(sub2[,"admin2"]!="Unknown") >= n2)
											{
												# print(c(i,"g"))
												possible_indices = which(sub2[,"admin2"]!="Unknown")								
												selected_indices = possible_indices[sample(1:length(possible_indices),n2,replace=F)]
												selected_sequences = c(selected_sequences, sub2[selected_indices,"label"])
											}
									}
							}
					}
				cat(genotypes[h],": ",length(selected_sequences)," selected sequences\n",sep="")
				tab3 = tab2[which(tab2[,"label"]%in%selected_sequences),]
				sequences1 = fasta2[which(!grepl(">",fasta2))]; ids1 = fasta2[which(grepl(">",fasta2))]
				sequences2 = sequences1[gsub(">","",ids1)%in%selected_sequences]
				ids2 = ids1[gsub(">","",ids1)%in%selected_sequences]; fasta3 = c()
				for (i in 1:length(ids2)) fasta3 = c(fasta3, ids2[i], sequences2[i])
				write.csv(tab3, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_3.csv"), row.names=F, quote=F)
				write(fasta3, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_3.fas"))
			}
		if ((genotypes[h] == "D11")&(!file.exist(paste0("Genotype_",genotypes[h],"/D11_DTA_samples/",genotypes[h],"_alignment_3.fas"))))
			{
				fasta2 = scan(paste0("Genotype_",genotypes[h],"/D11_DTA_samples/",genotypes[h],"_alignment_2.fas"), what="", sep="\n", quiet=T)
				tab2 = read.csv(paste0("Genotype_",genotypes[h],"/D11_DTA_samples/",genotypes[h],"_alignment_2.csv"), head=T, sep=",")
				locations = paste0(tab2[,"country"],"_",tab2[,"admin1"]); dates = tab2[,"collection_date"]; months = rep(NA, length(dates))
				for (i in 1:length(months)) months[i] = paste(unlist(strsplit(dates[i],"-"))[1],unlist(strsplit(dates[i],"-"))[2],sep="-")
				different_locations = unique(locations); different_locations = different_locations[order(different_locations)]
				different_months = unique(months); different_months = different_months[order(different_months)]
				categories = paste0(months,"$",locations); different_categories = unique(categories)
				selected_sequences = c(); all_indices = c(); all_selected_indices = c()
				for (i in 1:length(different_categories))
					{
						month = unlist(strsplit(different_categories[i],"\\$"))[1]; location = unlist(strsplit(different_categories[i],"\\$"))[2]
						indices = which((grepl(month,tab2[,"collection_date"]))&(paste0(tab2[,"country"],"_",tab2[,"admin1"])==location))
						all_indices = c(all_indices, indices); admins2 = tab2[indices,"admin2"]
						if (length(indices) == 1)
							{
								selected_indices = indices
							}	else	{
								selected_indices = indices[sample(c(1:length(indices)), size=n, replace=F)]
							}
						selected_sequences = c(selected_sequences, tab2[selected_indices,"label"])
						all_selected_indices = c(all_selected_indices, selected_indices)
					}
				cat(genotypes[h],": ",length(selected_sequences)," selected sequences\n",sep="")
				tab3 = tab2[which(tab2[,"label"]%in%selected_sequences),]
				sequences1 = fasta2[which(!grepl(">",fasta2))]; ids1 = fasta2[which(grepl(">",fasta2))]
				sequences2 = sequences1[gsub(">","",ids1)%in%selected_sequences]
				ids2 = ids1[gsub(">","",ids1)%in%selected_sequences]; fasta3 = c()
				for (i in 1:length(ids2)) fasta3 = c(fasta3, ids2[i], sequences2[i])
				write.csv(tab3, paste0("Genotype_",genotypes[h],"/D11_DTA_samples/",genotypes[h],"_alignment_3.csv"), row.names=F, quote=F)
				write(fasta3, paste0("Genotype_",genotypes[h],"/D11_DTA_samples/",genotypes[h],"_alignment_3.fas"))
			}
	}

	# 1.3. Post hoc correction in the specific case of B3.2 (due to error)

h = 1
tab2 = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_2.csv"), head=T, sep=",")
tab3 = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_3.csv"), head=T, sep=",")
fasta2 = scan(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_2.fas"), what="", sep="\n", quiet=T)
fasta3 = scan(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_3.fas"), what="", sep="\n", quiet=T)
sequencesToRemove = tab2[which((tab2[,"admin1"]=="")|(tab2[,"admin1"]=="Unknown")|(is.na(tab2[,"admin1"]))|(tab2[,"admin1"]=="Gough_Island")),"label"]
if (length(sequencesToRemove) > 0)
	{
		tab2 = tab2[which(!tab2[,"label"]%in%sequencesToRemove),]; tab3 = tab3[which(!tab3[,"label"]%in%sequencesToRemove),]
		sequences2 = fasta2[which(!grepl(">",fasta2))]; ids2 = fasta2[which(grepl(">",fasta2))]
		sequences3 = fasta3[which(!grepl(">",fasta3))]; ids3 = fasta3[which(grepl(">",fasta3))]
		sequences2 = sequences2[!gsub(">","",ids2)%in%sequencesToRemove]
		sequences3 = sequences3[!gsub(">","",ids3)%in%sequencesToRemove]
		ids2 = ids2[!gsub(">","",ids2)%in%sequencesToRemove]; fasta2 = c()
		ids3 = ids3[!gsub(">","",ids3)%in%sequencesToRemove]; fasta3 = c()
		for (i in 1:length(ids2)) fasta2 = c(fasta2, ids2[i], sequences2[i])
		for (i in 1:length(ids3)) fasta3 = c(fasta3, ids3[i], sequences3[i])
		write.csv(tab2, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_2.csv"), row.names=F, quote=F)
		write(fasta2, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_2.fas"))
		write.csv(tab3, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_3.csv"), row.names=F, quote=F)
		write(fasta3, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_3.fas"))
	}

	# 1.4. Subsampling the extented metadata file for each genotype accordingly

indices = c(3,2,2,3)
for (h in 1:length(genotypes))
	{
		fileName = paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_",indices[h],".fas")
		file.copy(fileName, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_skygrid.fasta"), overwrite=T)
		tab2a = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_",indices[h],".csv"), head=T, sep=",")
		if ("location_used"%in%colnames(tab2a))
			{
				tab2b = tab2a; colnames(tab2b)[1] = "trait"; tab2b = tab2b[c("trait","country","admin1","admin2","collection_date","latitude","longitude","location_used")]
			}	else	{
				tab2b = tab2a; colnames(tab2b)[1] = "trait"; tab2b = tab2b[c("trait","country","admin1","admin2","collection_date")]
			}
		write.csv(tab2b, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.csv"), row.names=F, quote=F)				
	}
	
	# 1.5. Retrieving the difference between the maximum and minimum collection dates

indices = c(3,2,2,3)
for (h in 1:length(genotypes))
	{
		tab = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_subsampling/",genotypes[h],"_alignment_",indices[h],".csv"), head=T, sep=",")
		dates = decimal_date(ymd(tab[,"collection_date"])); cat(genotypes[h],": ",max(dates)-min(dates),"\n",sep="")
	}

# 2. Preparing the XML and KML files for the continuous/discrete phylogeographic analyses

	# 2.1. Selecting 1,000 empirical trees

for (h in 1:length(genotypes))
	{
		if (file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_skygrid.trees")))
			{
				if (!file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_empirical.trees")))
					{
						all_trees = scan(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_skygrid.trees"), what="", sep="\n", quiet=T)
						index = which(grepl("\t\t;",all_trees)); index = index[length(index)]; txt = all_trees[1:index]; ratio = 0.1
						for (i in 1:length(txt)) txt = gsub("'","",gsub("\"","",txt))
						indices = which(grepl("tree STATE_",all_trees)); burnIn = round(ratio*(length(indices)))+1
						selected_trees = all_trees[sample(indices[(burnIn+1):length(indices)], nberOfExtractionFiles, replace=F)]
						write(c(txt,selected_trees,"End;"), paste0("Genotype_",genotypes[h],"/",genotypes[h],"_empirical.trees"))
					}
			}
	}

	# 2.2. Re-generating metadata file

tip_labels = list()
for (h in 1:length(genotypes))
	{
		if (!genotypes[h]%in%c("AC","MM"))
			{
				tip_labels[[h]] = read.nexus(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_empirical.trees"))$tip.label[[1]]
			}	else	{
				tip_labels[[h]] = read.nexus(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_run.trees"))$tip.label[[1]]
			}
		tip_labels[[h]] = gsub("'","",tip_labels[[h]])
		if (!file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_skygrid.fasta")))
			{
				txt2 = c()
				if (!genotypes[h]%in%c("AC","MM"))
					{
						txt1 = scan(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_skygrid.xml"), what="", sep="\n", quiet=T, blank.lines.skip=F)
					}	else	{
						txt1 = scan(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_run.xml"), what="", sep="\n", quiet=T, blank.lines.skip=F)
					}
				for (i in 1:length(tip_labels[[h]]))
					{
						index = which(gsub("'","",txt1) == paste0("\t\t\t<taxon idref=\"",tip_labels[[h]][i],"\"/>"))
						txt2 = c(txt2, paste0(">",tip_labels[[h]][i]), gsub("\t\t\t","",txt1[index+1]))
					}
				write(txt2, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_skygrid.fasta"))
			}
		if ((!genotypes[h]%in%c("AC","MM"))&(!file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.csv"))))
			{
				tab = matrix(nrow=length(tip_labels[[h]]), ncol=3); tab[,1] = tip_labels[[h]]; colnames(tab) = c("trait","country","admin1")
				txt1 = scan(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_skygrid.xml"), what="", sep="\n", quiet=T, blank.lines.skip=F)
				for (i in 1:length(tip_labels[[h]]))
					{
						index = which(gsub("'","",txt1) == paste0("\t\t<taxon id=\"",tip_labels[[h]][i],"\">"))
						txt2 = unlist(strsplit(gsub("\t\t\t\t","",txt1[index+3]),"-"))
						tab[i,"country"] = txt2[1]; tab[i,"admin1"] = txt2[2]
					}
				write.csv(tab, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.csv"), row.names=F, quote=F)	
			}
		if ((genotypes[h]%in%c("AC","MM"))&(!file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.txt"))))
			{
				tab = matrix(nrow=length(tip_labels[[h]]), ncol=4); tab[,1] = tip_labels[[h]]; colnames(tab) = c("trait","collection_date","latitude","longitude")
				for (i in 1:length(tip_labels[[h]]))
					{
						tab[i,"trait"] = tip_labels[[h]][i]; txts = unlist(strsplit(tip_labels[[h]][i],"\\|"))
						if (!grepl("ghost",tip_labels[[h]][i]))
							{
								tab[i,"collection_date"] = txts[length(txts)-3]
								if (txts[1] == "EPI_ISL_17468386") tab[i,"collection_date"] = "2023-03-24"
								if (txts[1] == "EPI_ISL_17805999") tab[i,"collection_date"] = "2023-02-08"
								if (txts[1] == "EPI_ISL_18690677") tab[i,"collection_date"] = "2023-03-18"
								tab[i,"latitude"] = txts[length(txts)-1]; tab[i,"longitude"] = txts[length(txts)]
							}	else	{
								xml = scan(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_run.xml"), what="", sep="\n", quiet=T, blank.lines.skip=F)
								index = which(xml==paste0("		<taxon id=\"",tip_labels[[h]][i],"\">"))
								temp = date_decimal(as.numeric(unlist(strsplit(xml[index+1],"\""))[2]))
								tab[i,"collection_date"] = unlist(strsplit(as.character(temp)," "))[1]
								tab[i,"latitude"] = as.numeric(gsub("\t","",xml[index+3]))
								tab[i,"longitude"] = as.numeric(gsub("\t","",xml[index+6]))
							}
					}
				write.table(tab, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.txt"), row.names=F, quote=F, sep="\t")
			}
	}

	# 2.3. Generating the metadata file

countries_admins0 = c("Canada","USA","Mexico","Colombia","Honduras","Costa Rica","Argentina","Bolivia",
					  "Chile","Brazil","Ecuador","Guatemala","Peru","Uruguay","Panama","Falkland Islands")
if (!file.exists("GADMs_0.rds"))
	{
		gadms0 = list()
		gadms0[[1]] = shapefile("All_shapefiles/GADM_files/GADM_CAN_0.shp")
		gadms0[[2]] = shapefile("All_shapefiles/GADM_files/GADM_USA_0.shp")
		gadms0[[3]] = shapefile("All_shapefiles/GADM_files/GADM_MEX_0.shp")
		gadms0[[4]] = shapefile("All_shapefiles/GADM_files/GADM_COL_0.shp")
		gadms0[[5]] = shapefile("All_shapefiles/GADM_files/GADM_HND_0.shp")
		gadms0[[6]] = shapefile("All_shapefiles/GADM_files/GADM_CRI_0.shp")
		gadms0[[7]] = shapefile("All_shapefiles/GADM_files/GADM_ARG_0.shp")
		gadms0[[8]] = shapefile("All_shapefiles/GADM_files/GADM_BOL_0.shp")
		gadms0[[9]] = shapefile("All_shapefiles/GADM_files/GADM_CHL_0.shp")
		gadms0[[10]] = shapefile("All_shapefiles/GADM_files/GADM_BRA_0.shp")
		gadms0[[11]] = shapefile("All_shapefiles/GADM_files/GADM_ECU_0.shp")
		gadms0[[12]] = shapefile("All_shapefiles/GADM_files/GADM_GTM_0.shp")
		gadms0[[13]] = shapefile("All_shapefiles/GADM_files/GADM_PER_0.shp")
		gadms0[[14]] = shapefile("All_shapefiles/GADM_files/GADM_URY_0.shp")
		gadms0[[15]] = shapefile("All_shapefiles/GADM_files/GADM_PAN_0.shp")
		gadms0[[16]] = shapefile("All_shapefiles/GADM_files/GADM_FLK_0.shp")
		saveRDS(gadms0, "GADMs_0.rds")
	}	else	{
		gadms0 = readRDS("GADMs_0.rds")
	}
countries_admins1 = c("Canada","USA","Mexico","Colombia","Honduras","Costa Rica","Argentina",
					  "Bolivia","Chile","Brazil","Ecuador","Guatemala","Peru","Uruguay","Cayman Islands")
if (!file.exists("GADMs_1.rds"))
	{
		gadms1 = list()
		gadms1[[1]] = shapefile("All_shapefiles/GADM_files/GADM_CAN_1.shp")
		gadms1[[2]] = shapefile("All_shapefiles/GADM_files/GADM_USA_1.shp")
		gadms1[[3]] = shapefile("All_shapefiles/GADM_files/GADM_MEX_1.shp")
		gadms1[[4]] = shapefile("All_shapefiles/GADM_files/GADM_COL_1.shp")
		gadms1[[5]] = shapefile("All_shapefiles/GADM_files/GADM_HND_1.shp")
		gadms1[[6]] = shapefile("All_shapefiles/GADM_files/GADM_CRI_1.shp")
		gadms1[[7]] = shapefile("All_shapefiles/GADM_files/GADM_ARG_1.shp")
		gadms1[[8]] = shapefile("All_shapefiles/GADM_files/GADM_BOL_1.shp")
		gadms1[[9]] = shapefile("All_shapefiles/GADM_files/GADM_CHL_1.shp")
		gadms1[[10]] = shapefile("All_shapefiles/GADM_files/GADM_BRA_1.shp")
		gadms1[[11]] = shapefile("All_shapefiles/GADM_files/GADM_ECU_1.shp")
		gadms1[[12]] = shapefile("All_shapefiles/GADM_files/GADM_GTM_1.shp")
		gadms1[[13]] = shapefile("All_shapefiles/GADM_files/GADM_PER_1.shp")
		gadms1[[14]] = shapefile("All_shapefiles/GADM_files/GADM_URY_1.shp")
		gadms1[[15]] = shapefile("All_shapefiles/GADM_files/GADM_CYM_1.shp")
		gadms1[[4]]@data[which(gadms1[[4]]@data[,"NAME_1"]=="Bolívar"),"NAME_1"] = "Bolivar"
		gadms1[[4]]@data[which(gadms1[[4]]@data[,"NAME_1"]=="Córdoba"),"NAME_1"] = "Cordoba"
		gadms1[[5]]@data[which(gadms1[[5]]@data[,"NAME_1"]=="Atlántida"),"NAME_1"] = "Atlantida"
		gadms1[[8]]@data[which(gadms1[[8]]@data[,"NAME_1"]=="Potosí"),"NAME_1"] = "Potosi"
		gadms1[[9]]@data[which(gadms1[[9]]@data[,"NAME_1"]=="Valparaíso"),"NAME_1"] = "Valparaiso"
		saveRDS(gadms1, "GADMs_1.rds")
	}	else	{
		gadms1 = readRDS("GADMs_1.rds")
	}
countries_admins2 = c("USA","Costa Rica","Uruguay")
if (!file.exists("GADMs_2.rds"))
	{
		gadms2 = list()
		gadms2[[1]] = shapefile("All_shapefiles/GADM_files/GADM_USA_2.shp")
		gadms2[[2]] = shapefile("All_shapefiles/GADM_files/GADM_CRI_2.shp")
		gadms2[[3]] = shapefile("All_shapefiles/GADM_files/GADM_URY_2.shp")
		saveRDS(gadms2, "GADMs_2.rds")
	}	else	{
		gadms2 = readRDS("GADMs_2.rds")
	}
for (h in 1:length(genotypes))
	{
		dir.create(file.path(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons")), showWarnings=F)
		if (genotypes[h] == "B32")
			{
				tab1 = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.csv"), head=T); tab2 = tab1
				colnames(tab2)[which(colnames(tab2)=="location_used")] = "precision"
				tab2[which(tab2[,"admin1"]=="NewYork"),"admin1"] = "New York"
				gadms1[[9]]@data[which(gadms1[[9]]@data[,"NAME_1"]=="Magallanes y Antártica Chilena"),"NAME_1"] = "Magallanes y Antartica Chilena"
				gadms1[[9]]@data[which(gadms1[[9]]@data[,"NAME_1"]=="Ñuble"),"NAME_1"] = "Nuble"
				gadms1[[9]]@data[which(gadms1[[9]]@data[,"NAME_1"]=="Bío-Bío"),"NAME_1"] = "Bio-Bio"
				gadms1[[9]]@data[which(gadms1[[9]]@data[,"NAME_1"]=="Tarapacá"),"NAME_1"] = "Tarapaca"	
				gadms2[[2]]@data[which(gadms2[[2]]@data[,"NAME_1"]=="Limón"),"NAME_1"] = "Limon"
				gadms2[[2]]@data[which(gadms2[[2]]@data[,"NAME_2"]=="Limón"),"NAME_2"] = "Limon"
				gadms2[[2]]@data[which(gadms2[[2]]@data[,"NAME_2"]=="Pococí"),"NAME_2"] = "Pococi"	
				tab2[which(tab2[,"admin1"]=="Magallanes y Antartica Chilena"),"admin1"] = "Magallanes y Antartica Chilena"
				tab2[which(tab2[,"admin1"]=="Region_del_Biobio"),"admin1"] = "Bio-Bio"
				tab2[which(tab2[,"admin1"]=="Region_Metropolitana_de_Santiago"),"admin1"] = "Santiago Metropolitan"
				tab2[which(tab2[,"admin1"]=="Tarapaca"),"admin1"] = "Tarapaca"
				tab2[which(tab2[,"admin1"]=="Departamento_del_Magdalena"),"admin1"] = "Magdalena"			
				tab2[which(tab2[,"admin2"]=="Caroll"),"admin2"] = "Carroll"
				tab2[which(tab2[,"admin2"]=="Kingsbury\xa0"),"admin2"] = "Kingsbury"
				tab2[which(tab2[,"admin2"]=="LaSalle"),"admin2"] = "La Salle"
				tab2[which(tab2[,"admin2"]=="Mcintosh"),"admin2"] = "McIntosh"
				tab2[which(tab2[,"admin2"]=="Nome Census Area"),"admin2"] = "Nome"
				tab2[which(tab2[,"admin2"]=="St_Johns"),"admin2"] = "Saint Johns"
				tab2[which(tab2[,"admin2"]=="St. Louis"),"admin2"] = "Saint Louis"
				tab2[which(tab2[,"admin2"]=="Valdez-Cordoba"),"admin2"] = "Valdez-Cordova"
				tab2[which(tab2[,"admin2"]=="Matanuska-Susitna Borough"),"admin2"] = "Matanuska-Susitna"
				tab2[which(tab2[,"admin2"]=="Mahnhomen"),"admin2"] = "Mahnomen"
				tab2[which(tab2[,"admin2"]=="Juan Mart\xedn de Pueyrred\xf3n"),"admin2"] = "Juan Martin de Pueyrredon"
				tab2[which(tab2[,"precision"]=="island"),"precision"] = NA # no coordinate provided
				tab2[which(tab2[,"admin2"]=="Metan"),"precision"] = "precise (Google_Earth)"
				tab2[which(tab2[,"admin2"]=="Metan"),c("latitude","longitude")] = c(-25.4960,-64.9730)
				tab2[which(tab2[,"admin1"]=="Beak_Island"),"precision"] = "precise (Google_Earth)"
				tab2[which(tab2[,"admin1"]=="Beak_Island"),c("latitude")] = c(-63.6139)
				tab2[which(tab2[,"admin1"]=="Beak_Island"),c("longitude")] = c(-57.3232)
				tab2[which(tab2[,"admin1"]=="Devil_Island"),"precision"] = "precise (Google_Earth)"
				tab2[which(tab2[,"admin1"]=="Devil_Island"),c("latitude","longitude")] = c(-63.8000,-57.2801)
				tab2[which(tab2[,"admin1"]=="Gough_Island"),"precision"] = "missing"
				tab2[which(tab2[,"admin1"]=="Gough_Island"),c("latitude","longitude")] = c(NA, NA)	
				tab2[which(tab2[,"admin1"]=="Hope_Bay"),"precision"] = "precise (Google_Earth)"
				tab2[which(tab2[,"admin1"]=="Hope_Bay"),c("latitude","longitude")] = c(-63.3833, -56.9833)				
				tab2[which(tab2[,"admin1"]=="King_George_Island"),"precision"] = "precise (Google_Earth)"
				tab2[which(tab2[,"admin1"]=="King_George_Island"),c("latitude")] = c(-61.9882)
				tab2[which(tab2[,"admin1"]=="King_George_Island"),c("longitude")] = c(-58.0196)			
				tab2[which(tab2[,"admin1"]=="Livingston_Island"),"precision"] = "precise (Google_Earth)"
				tab2[which(tab2[,"admin1"]=="Livingston_Island"),c("latitude","longitude")] = c(-62.6306,-60.2045)				
				tab2[which(tab2[,"admin1"]=="Robert_Island"),"precision"] = "precise (Google_Earth)"
				tab2[which(tab2[,"admin1"]=="Robert_Island"),c("latitude","longitude")] = c(-62.3998,-59.5074)				
				tab2[which(tab2[,"admin1"]=="South_Shetlands_Islands"),"precision"] = "precise (Google_Earth)"
				tab2[which(tab2[,"admin1"]=="South_Shetlands_Islands"),c("latitude","longitude")] = c(-62.5902,-59.9111)				
				tab2[which(tab2[,"admin1"]=="Torgersen_Island"),"precision"] = "precise (Google_Earth)"
				tab2[which(tab2[,"admin1"]=="Torgersen_Island"),c("latitude","longitude")] = c(-64.7667,-64.0833)
				for (i in 1:dim(tab2)[1])
					{
						# tab2[i,]
						tab2[i,"admin1"] = gsub(" ","_",tab2[i,"admin1"])
						tab2[i,"admin2"] = gsub(" ","_",tab2[i,"admin2"])
						if ((!is.na(tab2[i,"precision"]))&(tab2[i,"precision"] != ""))
							{
								tab2[i,"latitude"] = tab2[i,"latitude"]; tab2[i,"longitude"] = tab2[i,"longitude"]
							}	else	{
								tab2[i,c("latitude","longitude")] = NA
								if ((tab2[i,"admin1"] == "")|(tab2[i,"admin1"] == "Unknown")|(tab2[i,"admin1"] == "unknown"))
									{
										index0 = which(countries_admins0==gsub("_"," ",tab2[i,"country"])); maxArea = 0; index1 = 0; index2 = 0
										for (j in 1:length(gadms0[[index0]]@polygons))
											{
												for (k in 1:length(gadms0[[index0]]@polygons[[j]]@Polygons))
													{
														if (maxArea < gadms0[[index0]]@polygons[[j]]@Polygons[[k]]@area)
															{
																maxArea = gadms0[[index0]]@polygons[[j]]@Polygons[[k]]@area; index1 = j; index2 = k
															}
													}
											}
										pol = gadms0[[index0]]@polygons[[index1]]@Polygons[[index2]]
										if (!file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons/",gsub(" ","_",tab2[i,"country"]),".kml")))
											{
												sink(file=paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons/",gsub(" ","_",tab2[i,"country"]),".kml"))
												cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"); cat("\n")
												cat("<kml xmlns=\"http://earth.google.com/kml/2.2\">"); cat("\n")
												cat(paste("\t<polygon id=\"",gsub(" ","_",tab2[i,"country"]),"\" samplingProbability=\"",1,"\">",sep="")); cat("\n")
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
										pol = sps; proj4string(pol) = gadms2[[index1]]@proj4string
										random_point = spsample(pol, 1, type="random")@coords
										tab2[i,"latitude"] = random_point[1,2]; tab2[i,"longitude"] = random_point[1,1]
									}
								if ((tab2[i,"admin1"] != "")&(tab2[i,"admin1"] != "Unknown")&(tab2[i,"admin1"] != "unknown")&(tab2[i,"admin2"] != "")&(tab2[i,"admin2"] != "Unknown")&(tab2[i,"admin2"] != "unknown"))
									{
										index1 = which(countries_admins2==gsub("_"," ",tab2[i,"country"])); maxArea = 0; index3 = 0
										if (tab2[i,"country"] == "USA")
											{
												index2 = which((grepl(paste0("US.",tab2[i,"admin1"]),gadms2[[index1]]@data[,"HASC_2"]))&(gadms2[[index1]]@data[,"NAME_2"]==gsub("_"," ",tab2[i,"admin2"])))
											}	else	{
												index2 = which((gadms2[[index1]]@data[,"NAME_1"]==gsub("_"," ",tab2[i,"admin1"]))&(gadms2[[index1]]@data[,"NAME_2"]==gsub("_"," ",tab2[i,"admin2"])))
											}
										for (j in 1:length(gadms2[[index1]]@polygons[[index2]]@Polygons))
											{
												if (maxArea < gadms2[[index1]]@polygons[[index2]]@Polygons[[j]]@area)
													{
														maxArea = gadms2[[index1]]@polygons[[index2]]@Polygons[[j]]@area; index3 = j
													}
											}
										pol = gadms2[[index1]]@polygons[[index2]]@Polygons[[index3]]
										if (!file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons/",paste0(tab2[i,"country"],"_",tab2[i,"admin1"],"_",tab2[i,"admin2"]),".kml")))
											{
												sink(file=paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons/",paste0(tab2[i,"country"],"_",tab2[i,"admin1"],"_",tab2[i,"admin2"]),".kml"))
												cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"); cat("\n")
												cat("<kml xmlns=\"http://earth.google.com/kml/2.2\">"); cat("\n")
												cat(paste("\t<polygon id=\"",paste0(tab2[i,"country"],"_",tab2[i,"admin1"],"_",tab2[i,"admin2"]),"\" samplingProbability=\"",1,"\">",sep="")); cat("\n")
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
										pol = sps; proj4string(pol) = gadms2[[index1]]@proj4string
										random_point = spsample(pol, 1, type="random")@coords
										tab2[i,"latitude"] = random_point[1,2]; tab2[i,"longitude"] = random_point[1,1]
									}
								if ((tab2[i,"admin1"] != "")&(tab2[i,"admin1"] != "Unknown")&(tab2[i,"admin1"] != "unknown")&((tab2[i,"admin2"] == "")|(tab2[i,"admin2"] == "Unknown")|(tab2[i,"admin2"] == "unknown")))
									{
										index1 = which(countries_admins1==gsub("_"," ",tab2[i,"country"])); maxArea = 0; index3 = 0
										index2 = which(grepl(paste0("\\.",tab2[i,"admin1"]),gadms1[[index1]]@data[,"HASC_1"]))
										if ((length(index2) == 0)&(nchar(tab2[i,"admin1"]) > 2))
											{
												index2 = which(gadms1[[index1]]@data[,"NAME_1"] == gsub("_"," ",tab2[i,"admin1"]))
											}
										if ((length(index2) == 0)&(tab2[i,"country"] == "Canada")&(tab2[i,"admin1"] == "NL"))
											{
												index2 = which(grepl(paste0("\\.","NF"),gadms1[[index1]]@data[,"HASC_1"]))
											}
										for (j in 1:length(gadms1[[index1]]@polygons[[index2]]@Polygons))
											{
												if (maxArea < gadms1[[index1]]@polygons[[index2]]@Polygons[[j]]@area)
													{
														maxArea = gadms1[[index1]]@polygons[[index2]]@Polygons[[j]]@area; index3 = j
													}
											}
										pol = gadms1[[index1]]@polygons[[index2]]@Polygons[[index3]]
										if (!file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons/",paste0(tab2[i,"country"],"_",tab2[i,"admin1"]),".kml")))
											{
												sink(file=paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons/",paste0(tab2[i,"country"],"_",tab2[i,"admin1"]),".kml"))
												cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"); cat("\n")
												cat("<kml xmlns=\"http://earth.google.com/kml/2.2\">"); cat("\n")
												cat(paste("\t<polygon id=\"",paste0(tab2[i,"country"],"_",tab2[i,"admin1"]),"\" samplingProbability=\"",1,"\">",sep="")); cat("\n")
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
					}
				if (!file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.txt")))
					{
						write.table(tab2, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.txt"), row.names=F, quote=F, sep="\t")
					}
			}
		if ((genotypes[h] == "B36")|(genotypes[h] == "C21")|(genotypes[h] == "D11"))
			{
				tab1 = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.csv"), head=T)
				tab2 = tab1[,c("trait","country","admin1","admin2","collection_date")]; temp = matrix(nrow=dim(tab2)[1], ncol=2)
				colnames(temp) = c("latitude","longitude"); tab2 = cbind(tab2,temp)
				if (genotypes[h] == "B36")
					{
						tab2[which(tab2[,"admin2"]=="BigStone"),"admin2"] = "Big Stone"
						tab2[which(tab2[,"admin2"]=="DeSoto"),"admin2"] = "Desoto"
						tab2[which(tab2[,"admin2"]=="Golden-Valley"),"admin2"] = "Golden Valley"
						tab2[which(tab2[,"admin2"]=="LaMoure"),"admin2"] = "Lamoure"
						tab2[which(tab2[,"admin2"]=="Lewis"),"admin2"] = "Lewis and Clark"
						tab2[which(tab2[,"admin2"]=="Livingston Parish"),"admin2"] = "Livingston"
						tab2[which(tab2[,"admin2"]=="Noble"),"admin2"] = "Nobles"
						tab2[which(tab2[,"admin2"]=="Norfolk City"),"admin2"] = "Norfolk"
						tab2[which(tab2[,"admin2"]=="Otter"),"admin2"] = "Otter Tail"
						tab2[which(tab2[,"admin2"]=="Quantico"),"admin2"] = "Prince William"
						tab2[which(tab2[,"admin2"]=="St. Louis"),"admin2"] = "Saint Louis"
						tab2[,"admin2"] = gsub(" ","_",tab2[,"admin2"])
					}
				if (genotypes[h] == "C21")
					{					
						tab2[which(tab2[,"admin2"]=="Plum_Island"),"admin2"] = "Suffolk"
						tab2[which(tab2[,"admin2"]=="Stafford"),"admin2"] = "Strafford"
					}
				if (genotypes[h] == "D11")
					{
						tab2[which(grepl("St_",tab2[,"admin2"])),"admin2"] = gsub("St_","Saint_",tab2[which(grepl("St_",tab2[,"admin2"])),"admin2"])
						tab2[which((tab2[,"admin1"]=="WI")&(tab2[,"admin2"]=="Burnet")),"admin2"] = "Burnett"
						tab2[which(tab2[,"admin2"]=="DuPage"),"admin2"] = "Dupage"
						tab2[which(tab2[,"admin2"]=="LaMoure"),"admin2"] = "Lamoure"
						tab2[which(tab2[,"admin2"]=="Redacted"),"admin2"] = "Unknown"
						tab2[which(tab2[,"admin2"]=="Yukon"),"admin2"] = "Yukon-Koyukuk"
						tab2[which(tab2[,"admin2"]=="Dorchester"),"admin1"] = "MD"
						tab2[which(tab2[,"admin2"]=="Harford"),"admin1"] = "MD"
					}
				for (i in 1:dim(tab2)[1])
					{
						# tab2[i,]
						if ((tab2[i,"admin2"] != "")&(tab2[i,"admin2"] != "Unknown")&(tab2[i,"admin2"] != "unknown"))
							{
								maxArea = 0; index3 = 0
								index2 = which((grepl(paste0("US.",tab2[i,"admin1"]),gadms2[[index1]]@data[,"HASC_2"]))&(gadms2[[index1]]@data[,"NAME_2"]==gsub("_"," ",tab2[i,"admin2"])))
								for (j in 1:length(gadms2[[index1]]@polygons[[index2]]@Polygons))
									{
										if (maxArea < gadms2[[index1]]@polygons[[index2]]@Polygons[[j]]@area)
											{
												maxArea = gadms2[[index1]]@polygons[[index2]]@Polygons[[j]]@area; index3 = j
											}
									}
								pol = gadms2[[index1]]@polygons[[index2]]@Polygons[[index3]]
								if (!file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons/",paste0(tab2[i,"country"],"_",tab2[i,"admin1"],"_",tab2[i,"admin2"]),".kml")))
									{
										sink(file=paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons/",paste0(tab2[i,"country"],"_",tab2[i,"admin1"],"_",tab2[i,"admin2"]),".kml"))
										cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"); cat("\n")
										cat("<kml xmlns=\"http://earth.google.com/kml/2.2\">"); cat("\n")
										cat(paste("\t<polygon id=\"",paste0(tab2[i,"country"],"_",tab2[i,"admin1"],"_",tab2[i,"admin2"]),"\" samplingProbability=\"",1,"\">",sep="")); cat("\n")
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
								pol = sps; proj4string(pol) = gadms2[[index1]]@proj4string
								random_point = spsample(pol, 1, type="random")@coords
								tab2[i,"latitude"] = random_point[1,2]; tab2[i,"longitude"] = random_point[1,1]
							}
						if ((tab2[i,"admin2"] == "")|(tab2[i,"admin2"] == "Unknown")|(tab2[i,"admin2"] == "unknown"))
							{
								index1 = which(countries_admins1==gsub("_"," ",tab2[i,"country"])); maxArea = 0; index3 = 0
								index2 = which(grepl(paste0("\\.",tab2[i,"admin1"]),gadms1[[index1]]@data[,"HASC_1"]))
								if ((length(index2) == 0)&(nchar(tab2[i,"admin1"]) > 2))
									{
										index2 = which(grepl(gsub("_"," ",tab2[i,"admin1"]),gadms1[[index1]]@data[,"NAME_1"]))
									}
								if ((length(index2) == 0)&(tab2[i,"country"] == "Canada")&(tab2[i,"admin1"] == "NL"))
									{
										index2 = which(grepl(paste0("\\.","NF"),gadms1[[index1]]@data[,"HASC_1"]))
									}
								for (j in 1:length(gadms1[[index1]]@polygons[[index2]]@Polygons))
									{
										if (maxArea < gadms1[[index1]]@polygons[[index2]]@Polygons[[j]]@area)
											{
												maxArea = gadms1[[index1]]@polygons[[index2]]@Polygons[[j]]@area; index3 = j
											}
									}
								pol = gadms1[[index1]]@polygons[[index2]]@Polygons[[index3]]
								if (!file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons/",paste0(tab2[i,"country"],"_",tab2[i,"admin1"]),".kml")))
									{
										sink(file=paste0("Genotype_",genotypes[h],"/",genotypes[h],"_all_polygons/",paste0(tab2[i,"country"],"_",tab2[i,"admin1"]),".kml"))
										cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"); cat("\n")
										cat("<kml xmlns=\"http://earth.google.com/kml/2.2\">"); cat("\n")
										cat(paste("\t<polygon id=\"",paste0(tab2[i,"country"],"_",tab2[i,"admin1"]),"\" samplingProbability=\"",1,"\">",sep="")); cat("\n")
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
				if (!file.exists(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.txt")))
					{
						write.table(tab2, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.txt"), row.names=F, quote=F, sep="\t")
					}
			}
		if ((genotypes[h] == "A1")|(genotypes[h] == "A3"))
			{
				tab1 = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_metadata.csv"), head=T)
				tab2 = matrix(nrow=dim(tab1), ncol=5); colnames(tab2) = c("trait","country","admin1","latitude","longitude")
				tab2[,"trait"] = tab1[,"trait"]; tab2[,"country"] = tab1[,"country"]; tab2[,"admin1"] = tab1[,"admin1"]
				for (i in 1:dim(tab2)[1])
					{
						index1 = which(countries_admins1==tab2[i,"country"]); maxArea = 0; index3 = 0
						index2 = which(grepl(paste0(".",tab2[i,"admin1"]),gadms1[[index1]]@data[,"HASC_1"]))
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

	# 2.4. Editing the XML file with the polygons

for (h in 1:length(genotypes))
	{
		if (!genotypes[h]%in%c("AC","MM"))
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
												if ((!is.na(tab2[j,"precision"]))&(tab2[j,"precision"] != ""))
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
												if ((!is.na(tab2[j,"precision"]))&(tab2[j,"precision"] != ""))
													{
														usePolygon = FALSE
													}
											}
										if (usePolygon)
											{
												fileName = "MISSING_FILE_NAME"
												if ((is.na(tab[j,"admin1"]))|(tab[j,"admin1"] == "Unknown"))
													{
														fileName = tab[j,"country"]
													}
												if ((!is.na(tab[j,"admin1"]))&(tab[j,"admin1"] != "Unknown")&((is.na(tab[j,"admin2"]))|(tab[j,"admin2"] == "Unknown")))
													{
														fileName = paste0(tab[j,"country"],"_",gsub(" ","",tab[j,"admin1"]))
													}
												if ((!is.na(tab[j,"admin2"]))&(tab[j,"admin2"] != "Unknown"))
													{
														fileName = paste0(tab[j,"country"],"_",gsub(" ","",tab[j,"admin1"]),"_",gsub(" ","",tab[j,"admin2"]))
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
												if ((!is.na(tab2[j,"precision"]))&(tab2[j,"precision"] != ""))
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
												if ((!is.na(tab2[j,"precision"]))&(tab2[j,"precision"] != ""))
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
	}

# 3. Extracting the spatio-temporal information embedded in posterior trees

for (h in 1:length(genotypes))
	{
		if (!genotypes[h]%in%c("AC","MM"))
			{
				trees = scan(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_pols.trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
			}	else	{
				trees = scan(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_run.trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
			}
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
		if (!genotypes[h]%in%c("AC","MM"))
			{
				for (i in 1:dim(tab)[1]) collection_dates[i] = unlist(strsplit(tab[i,"trait"],"\\|"))[length(unlist(strsplit(tab[i,"trait"],"\\|")))-2]
			}	else	{
				collection_dates = tab[,"collection_date"]
			}
		mostRecentSamplingDatum = max(decimal_date(ymd(collection_dates))); print(date_decimal(mostRecentSamplingDatum))
		mcc_tre = readAnnotatedNexus(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000.tree"))
		mcc = mccTreeExtractions(mcc_tre, mostRecentSamplingDatum)
		write.csv(mcc, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000.csv"), row.names=F, quote=F)
		dir.create(file.path(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000_ext")), showWarnings=F)
		for (i in 1:length(trees))
			{
				tab = postTreeExtractions(trees[[i]], mostRecentSamplingDatum)
				write.csv(tab, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000_ext/TreeExtractions_",i,".csv"), row.names=F, quote=F)
			}
		if (genotypes[h] == "B32") # to discard a sample associated NA latitude and longitude values
			{
				mcc = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000.csv"), head=T)
				index = which(grepl("Gough_Island",mcc[,"tipLabel"]))
				if (length(index) > 0)
					{
						write.csv(mcc[-index,], paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000.csv"), row.names=F, quote=F)
					}
				for (i in 1:nberOfExtractionFiles)
					{
						tab = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000_ext/TreeExtractions_",i,".csv"), head=T)
						index = which(grepl("Gough_Island",tab[,"tipLabel"]))
						if (length(index) > 0)
							{
								write.csv(tab[-index,], paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000_ext/TreeExtractions_",i,".csv"), row.names=F, quote=F)
							}
					}
			}
	}
dir.create(file.path(paste0("Genotype_C21/C21_without_LDD")), showWarnings=F)
if (!file.exists(paste0("Genotype_C21/C21_without_LDD/TreeExtractions_1.csv")))
	{
		for (i in 1:nberOfExtractionFiles)
			{
				tab1 = read.csv(paste0("Genotype_C21/C21_RRW_1000_ext/TreeExtractions_",i,".csv"), head=T, sep=",")
				tab2 = tab1[which(tab1[,"endLon"]>-105),] # to discard one eastern LLD event
				write.csv(tab2, paste0("Genotype_C21/C21_without_LDD/TreeExtractions_",i,".csv"), row.names=F, quote=F)
			}
	}

# 4. Estimating lineage dispersal statistics from the continuous phylogeographic analyses

timeSlices = 100; onlyTipBranches = F; showingPlots = F; nberOfCores = 10; slidingWindow = 1/24
for (h in 1:length(genotypes))
	{
		localTreesDirectory = paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000_ext")
		outputName = paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_statistics/",genotypes[h],"_all") # "all" = all phylogenetic branches
		dir.create(file.path(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_statistics/")), showWarnings=F)
		if (genotypes[h] != "MM")
			{
				spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)
			}	else	{
				source("statisticsMM.r")
				statisticsMM(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)
			}
		if (genotypes[h] == "C21")
			{
				localTreesDirectory = paste0("Genotype_",genotypes[h],"/",genotypes[h],"_without_LDD")
				outputName = paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_statistics/",genotypes[h],"_woL") # "woL" = without LDD event
				spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)
			}
	}
t_start_wave_list = list(); t_max_dist_090_list = list()
if (!genotypes[1]%in%c("AC","MM"))
	{
		pdf(paste0("Wavefront_distances_NEW.pdf"), width=12, height=12)
		par(mfrow=c(4,1), oma=c(0,0,0,0), mar=c(2.5,4,0.5,0.5), mgp=c(0,0.1,0), lwd=0.3, bty="o", col="gray30")
		yMax = 16000; ats = seq(0,16000,5000)
	}	else	{
		pdf(paste0("Wavefront_distances_NEW.pdf"), width=12, height=3)
		par(mfrow=c(1,2), oma=c(0,0,0,0), mar=c(2.5,4,0.5,0.5), mgp=c(0,0.1,0), lwd=0.3, bty="o", col="gray30")
		yMax = 8000; ats = seq(0,10000,2000)
	}
colour_scale_1 = met.brewer(name="Hiroshige", n=111, type="continuous")[1:101]; minYear = 9999; maxYear = -9999
for (h in 1:length(genotypes))
	{
		mcc = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000.csv"), head=T, sep=",")
		if (minYear > min(mcc[,"startYear"])) minYear = min(mcc[,"startYear"])
		if (maxYear < max(mcc[,"endYear"])) maxYear = max(mcc[,"endYear"])	
	}
for (h in 1:length(genotypes))
	{
		swf_median = read.table(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_statistics/",genotypes[h],"_all_median_spatial_wavefront_distance.txt"), header=T)
		swf_95pHPD = read.table(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_statistics/",genotypes[h],"_all_95%HPD_spatial_wavefront_distance.txt"), header=T)
		if (genotypes[h] == "C21")
			{
				swf_median = read.table(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_statistics/",genotypes[h],"_woL_median_spatial_wavefront_distance.txt"), header=T)
				swf_95pHPD = read.table(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_statistics/",genotypes[h],"_woL_95%HPD_spatial_wavefront_distance.txt"), header=T)
			}
		swf = cbind(swf_median, swf_95pHPD[,2:3]); colnames(swf) = c("time","median","95pHDP_lower","95pHDP_upper") # "swf" is for "spatial wavefront"
		swf = swf[which((swf[,"time"]>=minYear)&(swf[,"time"]<maxYear)),]; max_dist_090 = max(swf_median[,"distance"])*0.90		
		t_start_wave = swf_median[which(swf_median[,"distance"]>0)[1],"time"]; t_start_wave_list[[h]] = t_start_wave
		t_max_dist_090 = swf_median[which(swf_median[,"distance"]>=max_dist_090)[1],"time"]; t_max_dist_090_list[[h]] = t_max_dist_090
		xMin = min(swf[,"time"]); xMax = max(swf[,"time"]); timeSlice = swf[1,"time"]-swf[2,"time"]; xMin = minYear; xMax = maxYear
		yMin = min(swf[,"95pHDP_lower"], na.rm=T); yMin = 0; # yMax = max(swf[,"95pHDP_upper"], na.rm=T)
		colours = paste0(colour_scale_1[(((swf[,c("time")]-minYear)/(maxYear-minYear))*(length(colour_scale_1)-1))+1],60)
		plot(swf[,"time"], swf[,"median"], lwd=0.7, type="l", cex.axis=0.8, cex.lab=0.8, col="gray30", axes=F, xlab=NA, ylab=NA, xlim=c(xMin,xMax), ylim=c(yMin,yMax))
		xx_l = c(swf[,c("time")],rev(swf[,c("time")])); yy_l = c(swf[,"95pHDP_lower"],rev(swf[,"95pHDP_upper"]))
		getOption("scipen"); opt = options("scipen"=20); polygon(xx_l,yy_l,col=rgb(187/255,187/255,187/255,0.25),border=0)
		for (j in 1:length(swf[,"time"]))
			{
				x1 = swf[j,"time"]-(timeSlice/2); x2 = swf[j,"time"]+(timeSlice/2)
				y1 = swf[j,"95pHDP_lower"]-1000; y2 = swf[j,"95pHDP_upper"]+1000
				polygon(c(x1,x2,x2,x1), c(y1,y1,y2,y2), col="gray90", border=NA)
				polygon(c(x1,x2,x2,x1), c(y1,y1,y2,y2), col=colours[j], border=NA)
			}
		getOption("scipen"); opt = options("scipen"=20); polygon(xx_l,yy_l,col=NA,border="gray30")
		lines(swf[,"time"], swf[,"median"], lwd=1.0, type="l", cex.axis=0.8, cex.lab=0.8, col="gray30")
		for (j in c(2021:2024)) abline(v=t_max_dist_090, lty=2, col="gray30", lwd=0.3)
		axis(side=1, lwd.tick=0.3, cex.axis=1.2, lwd=0.3, tck=-0.025, col="gray30", col.axis="gray30", col.tick="gray30", mgp=c(0,0.85,0), at=seq(2020,2026))
		axis(side=2, lwd.tick=0.3, cex.axis=1.2, lwd=0.3, tck=-0.025, col="gray30", col.axis="gray30", col.tick="gray30", mgp=c(0,0.75,0), at=ats)
		if (h == 1) mtext("Distance from epidemic origin (km)", side=2, col="gray30", cex=0.9, line=2.5, las=3)
	}
dev.off()
for (h in 1:length(genotypes))
	{
		wavefront_velocities_invasion_phase = rep(NA, nberOfExtractionFiles)
		spatial_wavefront_distances = read.table(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_statistics/",genotypes[h],"_all_spatial_wavefront_distances.txt"), head=T)
		if (genotypes[h] == "C21")
			{
				spatial_wavefront_distances = read.table(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_statistics/",genotypes[h],"_woL_spatial_wavefront_distances.txt"), head=T)
			}
		spatial_wavefront_distances = spatial_wavefront_distances[which((spatial_wavefront_distances[,"time"]>=t_start_wave_list[[h]])&
																		(spatial_wavefront_distances[,"time"]<=t_max_dist_090_list[[h]])),]
		for (i in 1:nberOfExtractionFiles)
			{
				duration = spatial_wavefront_distances[dim(spatial_wavefront_distances)[1],"time"]-spatial_wavefront_distances[1,"time"]
				wavefront_velocities_invasion_phase[i] = spatial_wavefront_distances[dim(spatial_wavefront_distances)[1],1+i]/duration
			}
		median = round(median(wavefront_velocities_invasion_phase)); hds = round(HDInterval::hdi(wavefront_velocities_invasion_phase)[1:2])
		cat("\t",genotypes[h]," wavefront velocity (invasion phase): ",median," km/year [",hds[1],"-",hds[2],"]\n",sep="") # 13315 km/y [11202-16676] when also considering the LDD for C2.1
	}
for (h in 1:length(genotypes))
	{
		dir.create(file.path(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_invasion_phase")), showWarnings=F)
		dir.create(file.path(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_after_invasion")), showWarnings=F)
		for (i in 1:nberOfExtractionFiles)
			{
				tab1 = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_",nberOfExtractionFiles,"_ext/TreeExtractions_",i,".csv"), head=T)
				tab2 = tab1[which(tab1[,"endYear"]<=t_max_dist_090_list[[h]]),]; tab3 = tab1[which(tab1[,"startYear"]>t_max_dist_090_list[[h]]),]
				write.csv(tab2, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_invasion_phase/TreeExtractions_",i,".csv"), row.names=F, quote=F)
				write.csv(tab3, paste0("Genotype_",genotypes[h],"/",genotypes[h],"_after_invasion/TreeExtractions_",i,".csv"), row.names=F, quote=F)
			}
	}
if ((showingPlots)&(!grepl("A",genotypes[1])))
	{
		dev.new(width=9, height=2.2); par(mfrow=c(1,4), mgp=c(0,0,0), oma=c(0.5,0.5,0,0), mar=c(3.0,3.0,0.5,0.5), lwd=0.3, col="gray30")
		for (h in 1:length(genotypes))
			{
				tab = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000_ext/TreeExtractions_",nberOfExtractionFiles,".csv"), head=T)
				d = rdist.earth(tab[,c("startLon","startLat")], tab[,c("endLon","endLat")], miles=F); d = diag(d); d2 = d^2
				t = tab[,"length"]; tx4 = 4*t; plot(t, log(d2/tx4), col="#217B9950", pch=16, cex=0.8, axes=F, ann=F, frame=T, xlim=c(0,0.35), ylim=c(0,log(230000000)))
				axis(side=1, lwd.tick=0.3, cex.axis=0.8, mgp=c(0,0.25,0), lwd=0.0, tck=-0.026, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(0,0.3,0.1))
				axis(side=2, lwd.tick=0.3, cex.axis=0.8, mgp=c(0,0.40,0), lwd=0.0, tck=-0.026, col.tick="gray30", col.axis="gray30", col="gray30")
				title(ylab=expression("Branch DC (km"^2*"/year, log)"), cex.lab=1, mgp=c(1.6,0,0), col.lab="gray30")
				title(xlab=expression("Time (years)"), cex.lab=1, mgp=c(1.5,0,0), col.lab="gray30")
			}
	}
timeSlices = 100; onlyTipBranches = F; showingPlots = F; nberOfCores = 10; slidingWindow = 1/24
for (h in 1:length(genotypes))
	{
		localTreesDirectory = paste0("Genotype_",genotypes[h],"/",genotypes[h],"_invasion_phase")
		outputName = paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_statistics/",genotypes[h],"_dIP") # "dIP" = during invasion phase
		spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)
		localTreesDirectory = paste0("Genotype_",genotypes[h],"/",genotypes[h],"_after_invasion")
		outputName = paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_statistics/",genotypes[h],"_aIP") # "aIP" = after invasion phase
		spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)
	 }
tab = matrix(nrow=length(genotypes), ncol=5); tab[,1] = paste0("H5N1 - ",genotypes)
colnames(tab) = c("Genotype","WDC (km2/day)","WDC during invasion phase","WDC after invasion phase","IBD signal")
for (h in 1:length(genotypes))
	{
		estimates = read.table(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_statistics/",genotypes[h],"_all_estimated_dispersal_statistics.txt"), head=T, sep="\t")
		vS = estimates[,"weighted_diffusion_coefficient"]; median = round(median(vS)/365.25); quantiles = round(HDInterval::hdi(vS)[1:2]/365.25)
		tab[h,"WDC (km2/day)"] = paste0(median," [",quantiles[1],"-",quantiles[2],"]")
		estimates = read.table(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_statistics/",genotypes[h],"_dIP_estimated_dispersal_statistics.txt"), head=T, sep="\t")
		vS = estimates[,"weighted_diffusion_coefficient"]; median = round(median(vS)/365.25); quantiles = round(HDInterval::hdi(vS)[1:2]/365.25)
		tab[h,"WDC during invasion phase"] = paste0(median," [",quantiles[1],"-",quantiles[2],"]")
		estimates = read.table(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_statistics/",genotypes[h],"_aIP_estimated_dispersal_statistics.txt"), head=T, sep="\t")
		vS = estimates[,"weighted_diffusion_coefficient"]; median = round(median(vS)/365.25); quantiles = round(HDInterval::hdi(vS)[1:2]/365.25)
		tab[h,"WDC after invasion phase"] = paste0(median," [",quantiles[1],"-",quantiles[2],"]")
		estimates = read.table(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_statistics/",genotypes[h],"_all_estimated_dispersal_statistics.txt"), head=T, sep="\t")
		vS = estimates[,"isolation_by_distance_signal_rP2"]; median = round(median(vS),3); quantiles = round(HDInterval::hdi(vS)[1:2],3)
		tab[h,"IBD signal"] = paste0(median," [",quantiles[1],"-",quantiles[2],"]")
	}
write.csv(tab, "Statistics_NEW.csv", row.names=F, quote=F)

# 5. Visualising the continuous phylogeographic reconstructions based on the polygons

e1 = extent(-175, -25, -77, 80); e2 = extent(-168, -52, 10, 74)
countries_1 = crop(shapefile("All_shapefiles/World_borders/World_countries_shapefile.shp"), e1)
borders_1 = crop(shapefile("All_shapefiles/Natural_Earth/International_borders/Only_international_borders.shp"), e1)
coasts_1 = crop(shapefile("All_shapefiles/Natural_Earth/Coast_lines_borders/Only_coast_lines_borders.shp"), e1)
countries_2 = crop(shapefile("All_shapefiles/World_borders/World_countries_shapefile.shp"), e2)
borders_2 = crop(shapefile("All_shapefiles/Natural_Earth/International_borders/Only_international_borders.shp"), e2)
coasts_2 = crop(shapefile("All_shapefiles/Natural_Earth/Coast_lines_borders/Only_coast_lines_borders.shp"), e2)
colour_scale_1 = met.brewer(name="Hiroshige", n=111, type="continuous")[1:101]
colour_scale_2 = met.brewer(name="Hokusai1", n=121, type="continuous")[11:111]

	# 5.1. Individual figures maps per each genotype

for (h in 1:length(genotypes))
	{
		mcc = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000.csv"), head=T, sep=",")
		if (genotypes[h] == "MM")
			{
				ghosts = matrix(0, nrow=dim(mcc)[1], ncol=2); colnames(ghosts) = c("startGhost","endGhost")
				ghosts[which(grepl("ghost",mcc[,"tipLabel"])),"endGhost"] = 1
				ghosts[which(grepl("ghost",mcc[,"tipLabel"])),"startGhost"] = 1
				ghostClades = TRUE
				while (ghostClades)
					{
						buffer = FALSE
						for (i in 1:dim(mcc)[1])
							{
								if ((ghosts[i,"endGhost"] == 0)&(mcc[i,"node2"]%in%mcc[,"node1"]))
									{
										indices = which(mcc[,"node1"]==mcc[i,"node2"])
										if (sum(ghosts[indices,"startGhost"]) == 2)
											{
												buffer = TRUE; ghosts[i,"endGhost"] = 1; ghosts[i,"startGhost"] = 1
											}
										if (sum(ghosts[indices,"startGhost"]) == 1)
											{
												buffer = TRUE; ghosts[i,"endGhost"] = 1
											}
									}
							}
						if (buffer == FALSE) ghostClades = FALSE
					}
			}
		minYear = min(mcc[,"startYear"]); maxYear = max(mcc[,"endYear"])
		endYears_indices = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
		endYears_colours = colour_scale_1[endYears_indices]
		if (!genotypes[h]%in%c("AC","MM"))
			{
				localTreesDirectory = paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000_ext")
				prob = 0.80; precision = 1/12; startDatum = minYear
				polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))				
				if (genotypes[h] == "C21") # to discard empty HPD region #7
					{
						polygons[[7]] = polygons[[8]]; polygons[[8]] = polygons[[9]]; polygons = polygons[1:8]
					}
			}	else	{
				mcc_tre = readAnnotatedNexus(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000.tree"))
				mcc_tab = mcc; polygons = suppressWarnings(spreadGraphic1(mcc_tre, mcc_tab))				
			}
		polygons_colours = rep(NA, length(polygons))
		for (i in 1:length(polygons))
			{
				date = as.numeric(names(polygons[[i]]))
				polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
				if (!genotypes[h]%in%c("AC","MM"))
					{
						polygons_colours[i] = paste0(colour_scale_1[polygon_index],"20")
					}	else	{
						polygons_colours[i] = paste0(colour_scale_1[polygon_index],"59")
					}
			}
		cutOffs = c(maxYear); g = 1; croppingPolygons = FALSE; plottingAllNodes = FALSE
		if ((genotypes[h] == "B32")|(genotypes[h] == "AC")|(genotypes[h] == "MM"))
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
						plot(rast, legend.only=T, add=T, col=colour_scale_1, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.15,0.50,0.280,0.287),
							 legend.args=list(text="", cex=0.7, col="gray30"), horizontal=T,
							 axis.args=list(cex.axis=0.5, lwd=0, lwd.tick=0.3, col.tick="gray30", tck=-1, col="gray30", col.axis="gray30", line=0, mgp=c(0,-0.15,0)))
		 			}
				for (i in 1:length(polygons))
					{
						if (as.numeric(names(polygons[[i]])) <= cutOffs[g]) plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
					}
				for (i in 1:dim(mcc)[1])
					{
						LTY = 1
						if (mcc[i,"endYear"] <= cutOffs[g])
							{
								if ((genotypes[h] == "MM")&&(ghosts[i,"startGhost"] == 1))
									{
										curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
						  		  		    		arr.width=0, lwd=0.2, lty=2, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
						  		  	}	else	{
										curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
						  		  		    		arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)						  		  		
						  		  	}
						  	}
					}
				for (i in dim(mcc)[1]:1)
					{
						if ((mcc[i,"startYear"] <= cutOffs[g])&(!mcc[i,"node1"]%in%mcc[,"node2"]))
							{
								startYears_index = (((mcc[i,"startYear"]-minYear)/(maxYear-minYear))*100)+1
								points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=colour_scale_1[startYears_index], cex=0.4)
								points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", lwd=0.2, cex=0.4)
							}
						if (mcc[i,"endYear"] <= cutOffs[g])
							{
								if ((genotypes[h] == "MM")&&(ghosts[i,"endGhost"] == 1))
									{
										points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col="white", cex=0.4)
										points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", lwd=0.2, cex=0.4)
									}	else	{
										points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.4)
										points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", lwd=0.2, cex=0.4)
									}
							}
					}
				rect(e1@xmin, e1@ymin, e1@xmax, e1@ymax, lwd=0.3, border="gray30")
				dev.off()
			}
		if ((genotypes[h] == "B36")|(genotypes[h] == "C21")|(genotypes[h] == "D11")|(genotypes[h] == "A1")|(genotypes[h] == "A3"))
			{	
				pdf(paste0(genotypes[h],"_RRW.pdf"), width=5.3, height=4.2)
				par(oma=c(0,0,0,0), mar=c(0.0,0.4,0.2,0.4), mgp=c(0,0.1,0), lwd=0.3, bty="o")
				plot(countries_2, col="gray90", border=NA, ann=F, axes=F)
				plot(borders_2, col="white", lwd=0.40, add=T)
				plot(coasts_2, col="gray70", lwd=0.40, add=T)
				rast = raster(matrix(nrow=1, ncol=2))
				rast[1] = minYear; rast[2] = maxYear
				if (g == length(cutOffs))
					{
						plot(rast, legend.only=T, add=T, col=colour_scale_1, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.08,0.40,0.130,0.139),
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
								points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=colour_scale_1[startYears_index], cex=0.4)
								points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", lwd=0.2, cex=0.4)
							}
						if (mcc[i,"endYear"] <= cutOffs[g])
							{
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.4)
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", lwd=0.2, cex=0.4)
							}
					}
				rect(e2@xmin, e2@ymin, e2@xmax, e2@ymax, lwd=0.3, border="gray30")
				dev.off()
			}
	}

	# 5.2. Composite figure with one map/genotype

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
		endYears_colours_list[[h]] = colour_scale_1[endYears_indices]
		localTreesDirectory = paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000_ext")
		prob = 0.80; precision = 1/12; startDatum = minYear
		polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))
		if (genotypes[h] == "C21") # to discard empty HPD region #9
			{
				polygons[[9]] = polygons[[10]]; polygons[[10]] = polygons[[11]]; polygons = polygons[1:10]
			}
		polygons_colours = rep(NA, length(polygons))
		for (i in 1:length(polygons))
			{
				date = as.numeric(names(polygons[[i]]))
				polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
				polygons_colours[i] = paste0(colour_scale_1[polygon_index],"20")
			}
		polygons_list[[h]] = polygons; polygons_colours_list[[h]] = polygons_colours
	}
rast = raster(matrix(nrow=1, ncol=2)); rast[1] = minYear; rast[2] = maxYear
for (h in 1:length(genotypes))
	{
		mcc = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000.csv"), head=T, sep=",")
		pdf(paste0("All_BCD_",genotypes[h],".pdf"), width=5.5, height=5.7)
		par(oma=c(0,0,0,0), mar=c(1,0.4,1,0.4), mgp=c(0,0.1,0), lwd=0.3, bty="o")
		plot(countries_1, col="gray90", border=NA, ann=F, axes=F)
		plot(borders_1, col="white", lwd=0.3, add=T)
		plot(coasts_1, col="gray70", lwd=0.5, add=T)
		plot(rast, legend.only=T, add=T, col=colour_scale_1, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.15,0.50,0.280,0.287),
			 legend.args=list(text="", cex=0.7, col="gray30"), horizontal=T,
			 axis.args=list(cex.axis=0.5, lwd=0, lwd.tick=0.3, col.tick="gray30", tck=-1, col="gray30", col.axis="gray30", line=0, mgp=c(0,-0.15,0)))
		for (i in 1:length(polygons_list[[h]]))
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
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=colour_scale_1[startYears_index], cex=0.4)
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", lwd=0.2, cex=0.4)
					}
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours_list[[h]][i], cex=0.4)
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", lwd=0.2, cex=0.4)
			}
		rect(e1@xmin, e1@ymin, e1@xmax, e1@ymax, lwd=0.3, border="gray30")
		dev.off()
	}

polygon_list = list(); coords_list = list(); prob = 0.80
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
		names(spdf) = genotypes[h]; polygon_list[[h]] = spdf; coords_list[[h]] = coords
	}
cols1 = list(); cols2 = list()
cols1[[1]] = rgb(222,67,39,220,maxColorValue=255); cols2[[1]] = rgb(222,67,39,100,maxColorValue=255) # red
cols1[[2]] = rgb(250,165,33,220,maxColorValue=255); cols2[[2]] = rgb(250,165,33,100,maxColorValue=255) # orange
cols1[[3]] = rgb(106,61,154,255,maxColorValue=255); cols2[[3]] = rgb(106,61,154,100,maxColorValue=255) # purple
cols1[[4]] = rgb(70,118,187,220,maxColorValue=255); cols2[[4]] = rgb(70,118,187,100,maxColorValue=255) # blue
pdf(paste0("MRCA_HPDs_NEW.pdf"), width=5.5, height=5.7)
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
for (h in 1:length(genotypes))
	{
		pdf(paste0(genotypes[h],"_ancestor.pdf"), width=5.5, height=5.7)
		par(oma=c(0,0,0,0), mar=c(1,0.4,1,0.4), mgp=c(0,0.1,0), lwd=0.3, bty="o")
		plot(countries_1, col="gray90", border=NA, ann=F, axes=F)
		plot(borders_1, col="white", lwd=0.3, add=T)
		plot(coasts_1, col="gray70", lwd=0.5, add=T)
		points(coords_list[[h]], axes=F, col=cols2[[h]], border=NA, add=T, pch=16, cex=0.3)
		rect(e1@xmin, e1@ymin, e1@xmax, e1@ymax, lwd=0.3, border="gray30")
		dev.off()
	}

	# 5.3. Successive snapshots for D1.1

h = 4; croppingPolygons = FALSE; plottingAllNodes = FALSE; mccHPDpolygons = TRUE
mcc = read.csv(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000.csv"), head=T, sep=",")
minYear = min(mcc[,"startYear"]); maxYear = max(mcc[,"endYear"])
endYears_indices = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
endYears_colours = colour_scale_2[endYears_indices]
if (!genotypes[h]%in%c("AC","MM"))
	{
		if (mccHPDpolygons == TRUE)
			{
				mcc_tre = readAnnotatedNexus(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000.tree"))
				mcc_tab = mcc; polygons = suppressWarnings(spreadGraphic1(mcc_tre, mcc_tab))				
			}
		if (mccHPDpolygons == FALSE)
			{
				localTreesDirectory = paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000_ext")
				prob = 0.80; precision = 1/12; startDatum = minYear
				polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))				
			}
		if (genotypes[h] == "C21") # to discard empty HPD region #7
			{
				polygons[[7]] = polygons[[8]]; polygons[[8]] = polygons[[9]]; polygons = polygons[1:8]
			}
	}	else	{
		mcc_tre = readAnnotatedNexus(paste0("Genotype_",genotypes[h],"/",genotypes[h],"_RRW_1000.tree"))
		mcc_tab = mcc; polygons = suppressWarnings(spreadGraphic1(mcc_tre, mcc_tab))				
	}
polygons_colours = rep(NA, length(polygons))
for (i in 1:length(polygons))
	{
		date = as.numeric(names(polygons[[i]]))
		polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
		if (!genotypes[h]%in%c("AC","MM"))
			{
				polygons_colours[i] = paste0(colour_scale_2[polygon_index],"20")
			}	else	{
				polygons_colours[i] = paste0(colour_scale_2[polygon_index],"59")
			}
		if (mccHPDpolygons == TRUE)
			{
				polygons_colours[i] = paste0(colour_scale_2[polygon_index],"0D")
			}
	}
dates1a = c("2024-07-23","2024-08-25","2024-10-19","2025-04-13"); dates1b = decimal_date(ymd(dates1a))
dates2a = c("2024-08-22","2024-09-25","2024-11-18","2025-05-13"); dates2b = decimal_date(ymd(dates2a))

pdf(paste0(genotypes[h],"_snapshots.pdf"), width=5.3*2, height=4.2*2)
par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(0.0,0.0,0.0,0.0), mgp=c(0,0.1,0), lwd=0.3, bty="o")
rast = raster(matrix(nrow=1, ncol=2)); rast[1] = minYear; rast[2] = maxYear
labels = c("2024-10-01","2025-01-01","2025-04-01","2025-07-01","2025-10-01")
ats = decimal_date(ymd(labels))
for (g in 1:length(dates1a))
	{	
		plot(countries_2, col="gray90", border=NA, ann=F, axes=F)
		plot(borders_2, col="white", lwd=0.40, add=T)
		plot(coasts_2, col="gray70", lwd=0.40, add=T)
		if (g == length(dates1a))
			{
				plot(rast, legend.only=T, add=T, col=colour_scale_2, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.08,0.40,0.130,0.139),
					 legend.args=list(text="", cex=0.7, col="gray30"), horizontal=T,
					 axis.args=list(cex.axis=0.45, lwd=0, lwd.tick=0.3, col.tick="gray30", tck=-1.0, col="gray30", col.axis="gray30", line=0,
					 mgp=c(0,-0.20,0), at=ats, label=labels))
 			}
		for (i in 1:length(polygons))
			{
				if ((as.numeric(names(polygons[[i]])) >= dates1b[g])&(as.numeric(names(polygons[[i]])) <= dates2b[g]))
					{
						plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
					}
			}
		for (i in 1:dim(mcc)[1])
			{
				if ((mcc[i,"endYear"] >= dates1b[g])&(mcc[i,"endYear"] <= dates2b[g]))
					{
						curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
				  		  		    arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
				  	}
			}
		for (i in 1:dim(mcc)[1])
			{
				if ((mcc[i,"endYear"] >= dates1b[g])&(mcc[i,"endYear"] <= dates2b[g])&(!mcc[i,"node1"]%in%mcc[,"node2"]))
					{
						startYears_index = (((mcc[i,"startYear"]-minYear)/(maxYear-minYear))*100)+1
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=colour_scale_2[startYears_index], cex=0.5)
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", lwd=0.2, cex=0.5)
					}
				if ((mcc[i,"endYear"] >= dates1b[g])&(mcc[i,"endYear"] <= dates2b[g]))
					{
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.5)
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", lwd=0.2, cex=0.5)
					}
			}
		rect(e2@xmin, e2@ymin, e2@xmax, e2@ymax, lwd=0.3, border="gray30")
	}
dev.off()

# 6. Analysing the impact of waterfowl flyways on the dispersal of D1.1 lineages in North America

flyways = shapefile("USA_flyways/Flyways_wo_Alaska.shp")

	# 6.1. Preparation of a discrete phylogeographic analysis (BSSVS analysis) based on the waterfowl flyways

flyway_assignation = read.csv("USA_flyways/Flyways_states_list.csv", head=T)
metadata = read.csv("Genotype_D11/D11_metadata.csv", head=T)
metadata$flyway = rep(NA, dim(metadata)[1])
for (i in 1:dim(metadata)[1])
	{
		index = which(flyway_assignation==paste0(metadata[i,"country"],"-",metadata[i,"admin1"]))
		if (paste0(metadata[i,"country"],"-",metadata[i,"admin1"]) == "Cayman_Islands-West_Bay")
			{
				index = which(flyway_assignation=="CaymanIslands")
			}
		metadata[i,"flyway"] = flyway_assignation[index,"flyway"]
	}
write.table(metadata[,c("trait","flyway")], "Genotype_D11/D11_flyways.txt", row.names=F, quote=F, sep="\t")

	# 6.2. Investigating the impact of waterfowl flyways on the dispersal frequency of viral lineages

for (i in 1:nberOfExtractionFiles)
	{
		tab = read.csv(paste0("Genotype_D11/D11_RRW_1000_ext/TreeExtractions_",i,".csv"), head=T)
		tab$startFlyway = rep(NA, dim(tab)[1]); tab$endFlyway = rep(NA, dim(tab)[1])
		for (j in 1:length(flyways))
			{
				for (k in 1:length(flyways@polygons[[j]]@Polygons))
					{
						pol_coords = flyways@polygons[[j]]@Polygons[[k]]@coords
						indices = which(point.in.polygon(tab[,"startLon"],tab[,"startLat"],pol_coords[,"x"],pol_coords[,"y"]) == 1)
						if (length(indices) > 0) tab[indices,"startFlyway"] = gsub(" Flyway","",flyways@data[j,"NAME"])
						indices = which(point.in.polygon(tab[,"endLon"],tab[,"endLat"],pol_coords[,"x"],pol_coords[,"y"]) == 1)
						if (length(indices) > 0) tab[indices,"endFlyway"] = gsub(" Flyway","",flyways@data[j,"NAME"])
					}
			}
		tab = tab[which((!is.na(tab$startFlyway))&(!is.na(tab$endFlyway))),]
		write.csv(tab, paste0("Genotype_D11/D11_RRW_USA_ext/TreeExtractions_",i,".csv"), row.names=F, quote=F)
	}

localTreesDirectory = "Genotype_D11/D11_RRW_USA_ext/"; randomisationDirectory = localTreesDirectory
envVariable = raster("Null_raster.tif"); randomProcedure = 3; repulsion = NULL; overwrite = FALSE; showingPlots = FALSE
treesRandomisation(localTreesDirectory, randomisationDirectory, nberOfExtractionFiles, envVariable, randomProcedure, 
				   repulsion, overwrite, showingPlots)

proportionOfFlywayChanges = matrix(nrow=nberOfExtractionFiles, ncol=2); colnames(proportionOfFlywayChanges) = c("obs","ran")
for (i in 1:nberOfExtractionFiles)
	{
		obs = read.csv(paste0("Genotype_D11/D11_RRW_USA_ext/TreeExtractions_",i,".csv"), head=T)
		ran = read.csv(paste0("Genotype_D11/D11_RRW_USA_ext/TreeRandomisation_",i,".csv"), head=T)
		proportionOfFlywayChanges[i,"obs"] = sum(obs[,"startFlyway"]!=obs[,"endFlyway"])/dim(obs)[1]
		proportionOfFlywayChanges[i,"ran"] = sum(ran[,"startFlyway"]!=ran[,"endFlyway"])/dim(ran)[1]
	}
p = sum(proportionOfFlywayChanges[,"obs"]<proportionOfFlywayChanges[,"ran"])/dim(proportionOfFlywayChanges)[1]
BF = (p/(1-p)) # 4.15

	# 6.3. Estimating the posterior probabilities for the flyway of origin of D1.1

flyways_of_origin = rep(NA, nberOfExtractionFiles)
for (i in 1:nberOfExtractionFiles)
	{
		tab = read.csv(paste0("Genotype_D11/D11_RRW_1000_ext/TreeExtractions_",i,".csv"), head=T)
		root_index = which(!tab[,"node1"]%in%tab[,"node2"])[1]
		for (j in 1:length(flyways))
			{
				for (k in 1:length(flyways@polygons[[j]]@Polygons))
					{
						pol_coords = flyways@polygons[[j]]@Polygons[[k]]@coords
						if (point.in.polygon(tab[root_index,"startLon"],tab[root_index,"startLat"],pol_coords[,"x"],pol_coords[,"y"]) == 1)
							{
								flyways_of_origin[i] = flyways@data[j,"NAME"]
							}
					}
			}
	}
different_flyways = unique(flyways_of_origin); different_flyways = different_flyways[order(different_flyways)]
flyways_of_origin = flyways_of_origin[which(!is.na(flyways_of_origin))]; different_flyways = different_flyways[which(!is.na(different_flyways))]
for (i in 1:length(different_flyways))
	{
		cat(different_flyways[i],", posterior probability = ",round(sum(flyways_of_origin==different_flyways[i])/nberOfExtractionFiles,3),"\n",sep="")
	}

