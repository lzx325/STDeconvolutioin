library(Matrix)
print_shape=function(data){
    for(nm in names(data)){
        v=data[[nm]]
        if(is.factor(v)||is.vector(v)){
            cat(nm,length(v),"n_levels:",length(unique(v)),"\n")
        } else if(is.matrix(v)||is.data.frame(v)||class(v)=="dgCMatrix"){
            cat(nm,dim(data[[nm]]),"\n")
        }
        
    }
    stopifnot(setequal(colnames(data$st_counts),rownames(data$st_locations)))
    stopifnot(setequal(colnames(data$sc_counts),names(data$sc_labels)))
}

load_MERFISH=function(res){
    st_obj_fp=sprintf("MERFISH/simMERFISH_%d.RDS",res)
    sc_counts_fp="MERFISH/raw_somatosensory_sc_exp.txt"
    sc_labels_fp="MERFISH/somatosensory_sc_labels.txt"
    expr=readRDS(st_obj_fp)
    st_counts=sparseMatrix(
        i=expr$sim$i,
        j=expr$sim$j,
        x=expr$sim$v,
        dims=c(expr$sim$nrow,expr$sim$ncol),
        dimnames=expr$sim$dimnames
    )
    st_counts=t(st_counts)
    st_locations=expr$st_location
    colnames(st_locations)=c("x","y")
    
    sc_counts=read.csv(sc_counts_fp,sep="\t",row.names=1)
    sc_labels=read.csv(sc_labels_fp,header=FALSE)$V1
    
    names(sc_labels)=colnames(sc_counts)
    
    ret_list=list(
        st_counts=st_counts,
        st_locations=st_locations,
        sc_counts=sc_counts,
        sc_labels=sc_labels
    )
    print_shape(ret_list)
    return(ret_list)
}

load_seqFISH=function(n_genes){
    st_counts_fp=sprintf("seqFISH+/Out_gene_expressions_%dgenes.csv",n_genes)
    st_locations_fp="seqFISH+/Out_rect_locations.csv"
    sc_counts_fp="seqFISH+/raw_somatosensory_sc_exp.txt"
    sc_labels_fp="seqFISH+/somatosensory_sc_labels.txt"
    
    st_counts=read.csv(st_counts_fp,sep=",",row.names=1)
    st_counts=t(st_counts)
    st_locations=read.csv(st_locations_fp,sep=",",row.name=1)
    st_locations=st_locations[,c("X","Y")]
    
    sc_counts=read.csv(sc_counts_fp,sep="\t",row.names=1)
    sc_labels=read.csv(sc_labels_fp,header=FALSE)$V1
    names(sc_labels)=colnames(sc_counts)
    
    ret_list=list(
            st_counts=st_counts,
            st_locations=st_locations,
            sc_counts=sc_counts,
            sc_labels=sc_labels
    )
    print_shape(ret_list)
    return(ret_list)
}

load_SlideseqV2=function(n_clusters){
    st_counts_fp="SlideseqV2/Hippocampus_MappedDGEForR.csv"
    st_locations_fp="SlideseqV2/Hippocampus_BeadLocationsForR.csv"
    sc_counts_fp="SlideseqV2/sc_count.csv"
    sc_labels_fp=sprintf("SlideseqV2/sc_%d_celltype.csv",n_clusters)
    
    st_counts=read.csv(st_counts_fp,sep=",",row.names=1)
    st_locations=read.csv(st_locations_fp,sep=",",row.names=1)
    
    sc_counts=read.csv(sc_counts_fp,sep=",",row.name=1)
    sc_labels_df=read.csv(sc_labels_fp,row.names=1)
    sc_labels=sc_labels_df[,1]
    names(sc_labels)=rownames(sc_labels_df)
    
    ret_list=list(
            st_counts=st_counts,
            st_locations=st_locations,
            sc_counts=sc_counts,
            sc_labels=sc_labels
    )
    print_shape(ret_list)
    return(ret_list)
}

load_stereo_seq=function(){
    st_obj_fp="stereo_seq/stereoseq_ob.RDS"
    sc_obj_fp="stereo_seq/scRNA_seq_ob.RDS"
    st_obj=readRDS(st_obj_fp)
    st_counts=t(st_obj$st_count)
    st_locations=st_obj$st_location
    
    sc_obj=readRDS(sc_obj_fp)
    sc_counts=sc_obj$sc_count

    sc_labels=as.factor(sc_obj$cell_type_anno$cell_type_anno)
    names(sc_labels)=rownames(sc_obj$cell_type_anno)
    
    ret_list=list(
            st_counts=st_counts,
            st_locations=st_locations,
            sc_counts=sc_counts,
            sc_labels=sc_labels
    )
    print_shape(ret_list)
    return(ret_list)
    
}

load_10Xvisium=function(){
    st_counts_fp="10Xvisium/st_count.csv"
    st_locations_fp="10Xvisium/st_location.csv"
    sc_counts_fp="10Xvisium/sc_mousebrain.csv"
    sc_labels_fp="10Xvisium/sc_celltype.csv"
    
    st_counts=read.csv(st_counts_fp,sep=",",row.names=1,check.names = FALSE)
    st_locations=read.csv(st_locations_fp,sep=",",row.names=1)
    
    sc_counts=read.csv(sc_counts_fp,sep=",",row.names=1)
    sc_counts=t(sc_counts)
    sc_labels_df=read.csv(sc_labels_fp,sep=",",row.names=1)
    sc_labels=sc_labels_df$annotation_1
    names(sc_labels)=rownames(sc_labels_df)
    
    ret_list=list(
            st_counts=st_counts,
            st_locations=st_locations,
            sc_counts=sc_counts,
            sc_labels=sc_labels
    )
    print_shape(ret_list)
    return(ret_list)
}

load_ST_PDAC=function(){
    data=readRDS("ST_PDAC/PDAC_GSM4100721.rds")
    sc_counts=data$sc_count

    sc_labels=as.factor(data$cell_type)
    names(sc_labels)=colnames(sc_counts)

    st_counts=data$st_count
    st_locations=data$st_location
    ret_list=list(
            st_counts=st_counts,
            st_locations=st_locations,
            sc_counts=sc_counts,
            sc_labels=sc_labels
    )
    print_shape(ret_list)
    return(ret_list)
    
}

