library(spacexr)
library(Matrix)
library(optparse)
library(yaml)
source("data_loader.R")

option_list = list(
    make_option(c("--max-cores","-t"), type = "integer", default=NULL, help="# threads"),
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
out_dir=file.path("./results_RCTD",opt$dataset)
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

sc_reference=Reference(
    counts=data$sc_counts,
    cell_types=data$sc_labels
)

st_data=SpatialRNA(
    counts=data$st_counts,
    coords=data$st_location,
    require_int=FALSE
)

start_time <- Sys.time()

myRCTD <- create.RCTD(st_data, sc_reference, max_cores = opt$`max-cores`, CELL_MIN_INSTANCE = 1)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')

end_time <- Sys.time()

weights=myRCTD@results$weights
norm_weights=normalize_weights(weights)
result_list$running_wtime=as.numeric(end_time - start_time, units="secs")

write_yaml(result_list,out_yaml_fp)
write.csv(as.matrix(weights),out_matrix_unnorm_fp)
write.csv(as.matrix(norm_weights),out_matrix_norm_fp)
