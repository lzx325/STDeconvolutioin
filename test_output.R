source("data_loader.R")
dataset_names=c("MERFISH","seqFISH","ST_PDAC","stereo_seq","10Xvisium","SlideseqV2")
out_dir="results_SpatialDecon"

for(dataset_name in dataset_names){
    data_loader_fn=environment()[[sprintf("load_%s",dataset_name)]]
    method_dir=file.path(out_dir,dataset_name)
    print(dataset_name)
    if(dataset_name == "MERFISH"){
        for(param in c(20,50,100)){
            data=data_loader_fn(param)
            out_matrix_norm_fp=file.path(method_dir,sprintf("%s_%d.norm.csv",dataset_name,param))
            out_matrix_norm=read.csv(out_matrix_norm_fp,row.names=1)
            print(dim(out_matrix_norm))
            stopifnot(all(rownames(out_matrix_norm) %in% rownames(data$st_locations)))

        }
        
    } else if(dataset_name == "seqFISH"){
        for(param in c(3000,6000,10000)){
            data=data_loader_fn(param)
            out_matrix_norm_fp=file.path(method_dir,sprintf("%s_%d.norm.csv",dataset_name,param))
            out_matrix_norm=read.csv(out_matrix_norm_fp,row.names=1)
            print(dim(out_matrix_norm))

            stopifnot(all(rownames(out_matrix_norm) %in% rownames(data$st_locations)))
        }
    } else if(dataset_name == "SlideseqV2"){
        for(param in c(11,17)){
            data=data_loader_fn(param)
            out_matrix_norm_fp=file.path(method_dir,sprintf("%s_%d.norm.csv",dataset_name,param))
            out_matrix_norm=read.csv(out_matrix_norm_fp,row.names=1)
            print(dim(out_matrix_norm))

            stopifnot(all(rownames(out_matrix_norm) %in% rownames(data$st_locations)))
        }
    } else {
        data=data_loader_fn()
        out_matrix_norm_fp=file.path(method_dir,sprintf("%s.norm.csv",dataset_name))
        out_matrix_norm=read.csv(out_matrix_norm_fp,row.names=1)
        print(dim(out_matrix_norm))
        stopifnot(all(rownames(out_matrix_norm) %in% rownames(data$st_locations)))
    }
}