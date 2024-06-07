suppressPackageStartupMessages(
c(
    library(Seurat),
    library(SingleCellExperiment),
    library(aricode),
    library(stringr),
    library(optparse),
    library(data.table)
 )
)
set.seed(102)

getLouvainClusters <- function(scData, numClusters = 6, resolution_start=0.1){
    currNumClusters <- 0
    resolution <- resolution_start
    
    scData <- Seurat::FindNeighbors(scData, reduction = "pca")
    
    while (currNumClusters < numClusters) {
        clusters <- Seurat::FindClusters(scData, resolution = resolution, algorithm=1)$seurat_clusters
        resolution <- resolution + 0.02
        currNumClusters <- length(unique(clusters))
    }
    print(paste0("Found: ", currNumClusters, " at resolution: ", resolution))
    return(clusters)
}


parser <- OptionParser(usage = "%prog [options] sample RDS",
                       option_list=list())
opt <- parse_args(parser, positional_arguments=4)
args <- opt$args

sample <- args[[1]]
model <- args[[2]]
cross_validation_combination <- args[[3]]
out_folder <- args[[4]]

rds_path <- paste0(out_folder, "/data/rds/", sample, ".rds")
cluster_out_path = paste0(out_folder, "/", cross_validation_combination, "/", model, "_evaluate/clusters_default/model-", sample, "-best_param.csv")


scData <- readRDS(rds_path)

clusters <- getLouvainClusters(scData, numClusters = length(unique(scData$ground_truth)))

dt <- data.table(
    spot_barcode=colnames(scData),
    louvain_cluster=clusters
)

fwrite(dt, cluster_out_path)