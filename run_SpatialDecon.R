library(SpatialDecon)
library(SeuratObject)
library(Matrix)
library(optparse)
library(yaml)
source("data_loader.R")

preprocess=function(data){
    st_counts_norm = sweep(data$st_counts, 2, colSums(data$st_counts), "/") * mean(colSums(data$st_counts))
    st_object=CreateSeuratObject(counts=st_counts_norm,assay="Spatial")
    stopifnot(setequal(colnames(st_object),rownames(data$st_location)))
    st_object=AddMetaData(st_object,data$st_locations[colnames(st_object),1],col.name="x")
    st_object=AddMetaData(st_object,data$st_locations[colnames(st_object),2],col.name="y")
    
    stopifnot(all(colnames(data$sc_counts)==names(data$sc_labels)))
    
    sc_counts_matrix=as.matrix(data$sc_counts)
    sc_counts_matrix=Matrix::Matrix((sc_counts_matrix),sparse=TRUE)
    sc_labels_df=data.frame(cell_barcodes=names(data$sc_labels),sc_labels=data$sc_labels)
    sc_matrix <- create_profile_matrix(
        mtx = sc_counts_matrix,            # cell x gene count matrix
        cellAnnots = sc_labels_df,  # cell annotations with cell type and cell name as columns 
        cellTypeCol = "sc_labels",  # column containing cell type
        cellNameCol = "cell_barcodes",           # column containing cell ID/name
        matrixName = "custom_cell_type_matrix", # name of final profile matrix
        outDir = NULL,                    # path to desired output directory, set to NULL if matrix should not be written
        normalize = TRUE,                # Should data be normalized? 
        minCellNum = 1,
        minGenes = 1
    ) 
    
    return(
        list(
            st_object=st_object,
            sc_matrix=sc_matrix
        )
    )
}

option_list = list(
    make_option(c("--dataset"), type = "character", default=NULL),
    make_option("--dataset_params", type = "integer", default=-1)
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
print(opt)
stopifnot(opt$dataset %in% c("MERFISH","seqFISH","SlideseqV2","stereo_seq","10Xvisium","ST_PDAC"))
result_list=opt
data_loader_name=sprintf("load_%s",opt$dataset)
result_list$data_loader_name=data_loader_name
data_loader_fn=(environment()[[data_loader_name]])

out_dir=file.path("./results_SpatialDecon",opt$dataset)
dir.create(out_dir,showWarnings=FALSE,recursive=TRUE)

if(opt$dataset %in% c("MERFISH","seqFISH","SlideseqV2")){
    stopifnot(opt$dataset_params!=-1)
    data=data_loader_fn(opt$dataset_params)
    out_matrix_unnorm_fp=file.path(out_dir,sprintf("%s_%d.unnorm.csv",opt$dataset,opt$dataset_params))
    out_matrix_norm_fp=file.path(out_dir,sprintf("%s_%d.norm.csv",opt$dataset,opt$dataset_params))
    out_yaml_fp=file.path(out_dir,sprintf("%s_%d.yaml",opt$dataset,opt$dataset_params))
} else {
    data=data_loader_fn()
    out_matrix_unnorm_fp=file.path(out_dir,sprintf("%s.unnorm.csv",opt$dataset))
    out_matrix_norm_fp=file.path(out_dir,sprintf("%s.norm.csv",opt$dataset))
    out_yaml_fp=file.path(out_dir,sprintf("%s.yaml",opt$dataset))
}
processed_data=preprocess(data)
start_time <- Sys.time()
res = runspatialdecon(object = processed_data$st_object,
                      bg = 0.01,
                      X = processed_data$sc_matrix,
                      align_genes = TRUE)
end_time <- Sys.time()

weights=t(res$beta)
norm_weights=sweep(weights, 1, rowSums(weights), "/")
result_list$running_wtime=as.numeric(end_time - start_time, units="secs")
write_yaml(result_list,out_yaml_fp)
write.csv(as.matrix(weights),out_matrix_unnorm_fp)
write.csv(as.matrix(norm_weights),out_matrix_norm_fp)