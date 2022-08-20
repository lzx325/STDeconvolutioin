source("data_loader.R")
dataset_names=c("MERFISH","seqFISH","ST_PDAC","stereo_seq","10Xvisium","SlideseqV2")
out_dir="results_RCTD"

pad_result_df=function(data,result_df){
    colnames_dropped=setdiff(levels(data$sc_labels),colnames(result_df))
    if(length(colnames_dropped)>0){
        cat("colnames not found: ",colnames_dropped,"\n",length(colnames_dropped),"in total","\n")
        result_df[,colnames_dropped]=0
    }
    rownames_dropped=setdiff(rownames(data$st_locations),rownames(result_df))
    if(length(rownames_dropped)>0){
        cat("padding rownames: ",length(rownames_dropped),"in total","\n")
        result_df[rownames_dropped,]=0
    }

    result_df=result_df[rownames(data$st_locations),levels(data$sc_labels)]
    return(result_df)
}
for(dataset_name in dataset_names){
    data_loader_fn=environment()[[sprintf("load_%s",dataset_name)]]
    method_dir=file.path(out_dir,dataset_name)
    print(dataset_name)
    if(dataset_name == "MERFISH"){
        for(param in c(20,50,100)){
            data=data_loader_fn(param)
            out_matrix_norm_fp=file.path(method_dir,sprintf("%s_%d.norm.csv",dataset_name,param))
            out_matrix_norm=read.csv(out_matrix_norm_fp,row.names=1,check.names=FALSE)
            print(dim(out_matrix_norm))
            stopifnot(all(rownames(out_matrix_norm) %in% rownames(data$st_locations)))
            print(colnames(out_matrix_norm))
            print(levels(data$sc_labels))
            out_matrix_pad_fp=file.path(method_dir,sprintf("%s_%d.pad.csv",dataset_name,param))
            out_matrix_pad=pad_result_df(data,out_matrix_norm)
            write.csv(out_matrix_pad,out_matrix_pad_fp)

        }
        
    } else if(dataset_name == "seqFISH"){
        for(param in c(3000,6000,10000)){
            data=data_loader_fn(param)
            out_matrix_norm_fp=file.path(method_dir,sprintf("%s_%d.norm.csv",dataset_name,param))
            out_matrix_norm=read.csv(out_matrix_norm_fp,row.names=1,check.names=FALSE)
            print(dim(out_matrix_norm))

            stopifnot(all(rownames(out_matrix_norm) %in% rownames(data$st_locations)))

            print(colnames(out_matrix_norm))
            print(levels(data$sc_labels))

            out_matrix_pad_fp=file.path(method_dir,sprintf("%s_%d.pad.csv",dataset_name,param))
            out_matrix_pad=pad_result_df(data,out_matrix_norm)
            write.csv(out_matrix_pad,out_matrix_pad_fp)

        }
    } else if(dataset_name == "SlideseqV2"){
        for(param in c(11,17)){
            data=data_loader_fn(param)
            out_matrix_norm_fp=file.path(method_dir,sprintf("%s_%d.norm.csv",dataset_name,param))
            out_matrix_norm=read.csv(out_matrix_norm_fp,row.names=1,check.names=FALSE)
            print(dim(out_matrix_norm))

            stopifnot(all(rownames(out_matrix_norm) %in% rownames(data$st_locations)))
            print(colnames(out_matrix_norm))
            print(levels(data$sc_labels))
            out_matrix_pad_fp=file.path(method_dir,sprintf("%s_%d.pad.csv",dataset_name,param))
            out_matrix_pad=pad_result_df(data,out_matrix_norm)
            write.csv(out_matrix_pad,out_matrix_pad_fp)

        }
    } else {
        data=data_loader_fn()
        out_matrix_norm_fp=file.path(method_dir,sprintf("%s.norm.csv",dataset_name))
        out_matrix_norm=read.csv(out_matrix_norm_fp,row.names=1,check.names=FALSE)
        print(dim(out_matrix_norm))
        stopifnot(all(rownames(out_matrix_norm) %in% rownames(data$st_locations)))
        print(colnames(out_matrix_norm))
        print(levels(data$sc_labels))
        out_matrix_pad_fp=file.path(method_dir,sprintf("%s.pad.csv",dataset_name))
        out_matrix_pad=pad_result_df(data,out_matrix_norm)
        write.csv(out_matrix_pad,out_matrix_pad_fp)
    }
    


}