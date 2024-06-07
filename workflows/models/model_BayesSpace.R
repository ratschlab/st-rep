suppressPackageStartupMessages(
c(
    library(Seurat),
    library(SingleCellExperiment),
    library(aricode),
    library(stringr),
    library(optparse),
    library(data.table),
    library(mclust),
    library(BayesSpace),
    library(yaml)
 )
)


parser <- OptionParser(usage = "%prog [options] sample RDS",
                       option_list=list())
opt <- parse_args(parser, positional_arguments=7)
args <- opt$args

sample <- args[[1]]
param_name <- args[[2]]
img_format <- args[[3]]
cross_validation_combination <- args[[4]]
model <- args[[5]]
mode <- args[[6]]
out_folder <- args[[7]]


rds_path <- paste0(out_folder, "/data/rds/", sample, ".rds")
parameter_path <- paste0(out_folder, "/", cross_validation_combination, "/", model, "_fine_tune/parameters/", param_name, ".yaml")
cluster_out_path <- paste0(out_folder, "/", cross_validation_combination, "/", model, "_", mode, "/clusters/model-", sample, "-", param_name, ".csv")


parameters <- yaml.load_file(parameter_path)

scData <- readRDS(rds_path)
scData <- as.SingleCellExperiment(scData)


scData@colData$row = scData@colData$x_array
scData@colData$col = scData@colData$y_array

scData <- spatialPreprocess(scData, platform="Visium")

parameters$sce <- scData
parameters$q <- length(unique(scData$ground_truth))


scData <- do.call(spatialCluster, parameters)


dt <- data.table(
    spot_barcode=colnames(scData),
    bayesspace_cluster=scData$spatial.cluster
)

fwrite(dt, cluster_out_path)