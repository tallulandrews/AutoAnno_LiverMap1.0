metadata <- read.table("Cell_clusterID_cycle.txt",  header=T)
mat <- read.table("FilteredCounts.csv", sep=",", header=T)
rownames(mat) <- mat[,1]
mat <- mat[,-1]
require("Matrix")
mat <- as(mat, "Matrix")

mat <- mat[!grepl("^MT-", rownames(mat)),] #remove mt genes

sf <- Matrix::colSums(mat);
norm <- t(t(mat)/sf*median(sf));

binary <- mat
binary[binary > 1] <- 1;

require(CellTypeProfiles)
cluster_means <- CellTypeProfiles::my_row_mean_aggregate(mat, metadata$Cluster)
cluster_norm <- CellTypeProfiles::my_row_mean_aggregate(norm, metadata$Cluster)
cluster_detection <- CellTypeProfiles::my_row_mean_aggregate(binary, metadata$Cluster)

table(metadata$Cluster, metadata$Sample)

max_diff <- function(mat) {
	on_off <- matrix(0, ncol=ncol(mat), nrow=nrow(mat));
	my_split_max_gap <- function(x) {
		x <- sort(x)
		jumps <- diff(x);
		br_pt <- which(jumps == max(jumps))
		return(c(x[br_pt], max(jumps)));
	}
	thresh <- apply(mat, 1, my_split_max_gap);
	on_off <- t(sapply(1:ncol(thresh), function(i) {mat[i,] > thresh[1,i]}))
	return(list(score=thresh[2,], on_off=on_off));
}

# New genes
mean_diff <- max_diff(cluster_norm)
detect_diff <- max_diff(cluster_detection)

good <- mean_diff$score > 0.3 | detect_diff$score > 0.1
unique <- rowSums(mean_diff$on_off) == 1 & rowSums(detect_diff$on_off) == 1
agree <- apply((mean_diff$on_off+detect_diff$on_off),1,function(x){sum(x==1)==0})

gene_lab <- rep("None", length(good));
T_mark <- (rowSums(mean_diff$on_off[,c(2,9,18)]) > 0 & rowSums(mean_diff$on_off[,-c(2,9,18)])==0) | (rowSums(detect_diff$on_off[,c(2,9,18)]) > 0 & rowSums(detect_diff$on_off[,-c(2,9,18)])==0)
B_mark <- (rowSums(mean_diff$on_off[,c(7,16)]) > 0 & rowSums(mean_diff$on_off[,-c(7,16)])==0) | (rowSums(detect_diff$on_off[,c(7,16)]) > 0 & rowSums(detect_diff$on_off[,-c(7,16)])==0)
mac_mark <- (rowSums(mean_diff$on_off[,c(4,10)]) > 0 & rowSums(mean_diff$on_off[,-c(4,10)])==0) | (rowSums(detect_diff$on_off[,c(4,10)]) > 0 & rowSums(detect_diff$on_off[,-c(4,10)])==0)
lsecs_mark <- (rowSums(mean_diff$on_off[,c(11,12)]) > 0 & rowSums(mean_diff$on_off[,-c(11,12)])==0) | (rowSums(detect_diff$on_off[,c(11,12)]) > 0 & rowSums(detect_diff$on_off[,-c(11,12)])==0)
centhep_mark <- (rowSums(mean_diff$on_off[,c(1,3)]) > 0 & rowSums(mean_diff$on_off[,-c(1,3)])==0) | (rowSums(detect_diff$on_off[,c(1,3)]) > 0 & rowSums(detect_diff$on_off[,-c(1,3)])==0)
porthep_mark <- (rowSums(mean_diff$on_off[,c(5,14)]) > 0 & rowSums(mean_diff$on_off[,-c(5,14)])==0) | (rowSums(detect_diff$on_off[,c(5,14)]) > 0 & rowSums(detect_diff$on_off[,-c(5,14)])==0)
hep_mark <- (rowSums(mean_diff$on_off[,c(1,3,5,6,14,15)]) > 0 & rowSums(mean_diff$on_off[,-c(1,3,5,6,14,15)])==0) | (rowSums(detect_diff$on_off[,c(1,3,5,6,14,15)]) > 0 & rowSums(detect_diff$on_off[,-c(1,3,5,6,14,15)])==0)


gene_lab[hep_mark & good] <- "Hep"
gene_lab[centhep_mark & good] <- "CentralHep"
gene_lab[porthep_mark & good] <- "PortHep"
gene_lab[lsecs_mark & good] <- "LSECs"
gene_lab[mac_mark & good] <- "Mac"
gene_lab[T_mark & good] <- "any_T_cell"
gene_lab[B_mark & good] <- "any_B_cell"

# cluster specific
clust_labs <- read.delim("Cluster_names.txt", sep="\t", header=F)

gene_lab[good & unique & agree] <- as.character(clust_labs[apply(mean_diff$on_off[good & unique & agree,],1,function(x){colnames(mean_diff$on_off)[x]}),2])

out <- data.frame(gene=rownames(cluster_means), label=gene_lab);
write.table(out, file="my_marker_genes.txt", row.names=F, col.names=T)

write.table(cluster_means, file="raw_cluster_means.txt", row.names=T, col.names=T)
write.table(cluster_detection, file="cluster_detect.txt", row.names=T, col.names=T)
write.table(cluster_norm, file="norm_cluster_means.txt", row.names=T, col.names=T)


# old genes:

mean_score <- apply(cluster_means, 1, max)- apply(cluster_means, 1, function(x){x[x==max(x)]<- -100; max(x)})
norm_score <- apply(cluster_norm, 1, max)- apply(cluster_norm, 1, function(x){x[x==max(x)]<- -100; max(x)})
detect_score <- apply(cluster_detection, 1, max)- apply(cluster_detection, 1, function(x){x[x==max(x)]<- -100; max(x)})

mean_id <- apply(cluster_means, 1, function(x){colnames(cluster_means)[which(x==max(x))[1]]})
norm_id <- apply(cluster_norm, 1, function(x){colnames(cluster_norm)[which(x==max(x))[1]]})
detect_id <- apply(cluster_detect, 1, function(x){colnames(cluster_detect)[which(x==max(x))[1]]})


summary(score)
score[score > 1]
table(cluster_id[score > 1])
table(cluster_id[score > 0.5])
savehistory("Get_ref_markers_original_map.Rhistory")
table(cluster_id[score > 0.25])
out <- data.frame(gene=names(cluster_id)[score > 0.25], cluster=cluster_id[score > 0.25], score=score[score > 0.25])
write.table(out, "ref_markers_original_map.out", row.names=F, col.names=F)
savehistory("Get_ref_markers_original_map.Rhistory")
