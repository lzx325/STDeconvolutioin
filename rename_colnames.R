ST_PDAC_namemap=c(
"Ductal - APOL1 high/hypoxic",
"Ductal - CRISP3 high/centroacinar like"
)



names(ST_PDAC_namemap)=c(
"Ductal - APOL1 high or hypoxic",
"Ductal - CRISP3 high or centroacinar like"
)

stereo_seq_namemap=c(
    "n15-M/TC-1",
    "n16-M/TC-2",
    "n17-M/TC-3"
)

names(stereo_seq_namemap)=c(
    "n15-M or TC-1",
    "n16-M or TC-2",
    "n17-M or TC-3"
)

out_dir="results_RCTD"
stereo_seq_fp_list=c(
    "stereo_seq/stereo_seq.norm.csv.old",
    "stereo_seq/stereo_seq.unnorm.csv.old"
)

ST_PDAC_fp_list=c(
    "ST_PDAC/ST_PDAC.norm.csv.old",
    "ST_PDAC/ST_PDAC.unnorm.csv.old"
)

for(fp in stereo_seq_fp_list){
    fp=file.path(out_dir,fp)
    table=read.csv(fp,row.names=1,check.names=FALSE)
    for(i in 1:ncol(table)){
        if(colnames(table)[i] %in% names(stereo_seq_namemap)){
            colnames(table)[i]=stereo_seq_namemap[colnames(table)[i]]
        }
    }

    newfp=gsub("\\.old$","",fp)
    write.csv(table,newfp)
}

for(fp in ST_PDAC_fp_list){
    fp=file.path(out_dir,fp)
    table=read.csv(fp,row.names=1,check.names=FALSE)
    for(i in 1:ncol(table)){
        if(colnames(table)[i] %in% names(ST_PDAC_namemap)){
            colnames(table)[i]=ST_PDAC_namemap[colnames(table)[i]]
        }
    }
    newfp=gsub("\\.old$","",fp)
    write.csv(table,newfp)
}