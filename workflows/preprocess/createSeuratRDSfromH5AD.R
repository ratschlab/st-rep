Sys.setenv(RETICULATE_PYTHON = "/cluster/customapps/biomed/grlab/users/knonchev/miniconda3/envs/nonchev/bin/python")
suppressPackageStartupMessages(
c(
    library(Seurat),
    library(SingleCellExperiment),
    library(SeuratData),
    #library(SeuratDisk),
    #library(zellkonverter),
    library(reticulate),
    library(BiocParallel),
    library(scuttle),
    library(scran),
    library(scater),
    library(aricode),
    library(stringr),
    library(optparse)
    
 )
)

set.seed(2023)



parser <- OptionParser(usage = "%prog [options] sample RDS",
                       option_list=list())
opt <- parse_args(parser, positional_arguments=2)
args <- opt$args

SAMPLE <- args[[1]]
out_folder <- args[[2]]


h5ad_preprocess_in <- paste0(out_folder, "/data/h5ad/", SAMPLE, ".h5ad")
rds_preprocess_out <- paste0(out_folder, "/data/rds/", SAMPLE, ".rds")

use_condaenv("nonchev")
py_discover_config()


createSeuratObjectFromAndata <- function(andata_path, rds_out){
    sc <- import("scanpy")
    aestetik <- import("aestetik")
    adata <- sc$read_h5ad(andata_path)
    
    adata <- aestetik$utils$utils_transcriptomics$preprocess_adata(adata)
    adata_raw = adata["raw"]
    
    colData <- as.data.frame(adata["obs"])
    rowData <- as.data.frame(adata_raw["var"])

    print("Extract matrix...")
    counts <- t(adata_raw["X"])
    rownames(counts) <- rowData$gene_ids
    colnames(counts) <- rownames(colData)
    print("Create sce...")
    scData <- CreateSeuratObject(counts=counts)
    scData <- AddMetaData(
          object = scData,
          metadata = colData
        )

    emb <- adata["obsm"]["X_pca"]
    rownames(emb) <- colnames(scData)
    scData[['pca']] <- CreateDimReducObject(embeddings = emb, key = 'pca_', assay = 'RNA')



    print("Save scData...")
    saveRDS(scData, rds_out)
}

createSeuratObjectFromAndata(h5ad_preprocess_in, rds_preprocess_out)