#-------------for GBM landscape analysis of hypoxia related genes:
library(survival)
library(KMsurv)
library(RColorBrewer)
library(glmnet)
library(survivalROC)
library(survminer)
library(pheatmap)
library(survcomp)
library(parallel)
library(clustree)# for hclust to plot tree
#------
library(limma)
library(data.table)
library(hgu133plus2.db)
library(hgu133a.db);
library(hgug4112a.db);#Agilent-014850 Whole Human Genome Microarray 4x44K G4112F 
library(annotate)
library(affy);
library("illuminaHumanv3.db")#for illuminaHumanWG-v3,Illumina HumanHT12v3 
library("illuminaHumanv2.db")#for Illumina human-6 v2.0
library(illuminaHumanv4.db)#for Illumina HumanHT-12 V4.0 expression beadchip
process_cel_data<-function(celpath,method){
	ReadAffy(celfile.path=celpath)->celdat;
	if(method=="mas5"){
		mas5(celdat)->celdat_normalized;
	}else if(method=="rma"){
		rma(celdat)->celdat_normalized;
	}
	exprs(celdat_normalized)->celdat_normalized_exp
	colnames(celdat_normalized_exp)->expd_names;
	unlist(lapply(expd_names,function(x){unlist(strsplit(x,split="\\."))->res;res[1]}))->expd_names;
	expd_names->colnames(celdat_normalized_exp);
	data.frame("DATA"=rownames(celdat_normalized_exp),celdat_normalized_exp)->celdat_normalized_exp;
	#scale(celdat_normalized_exp,center=rep(1000,ncol(celdat_normalized_exp)))->celdat_mas5_exp_scaled;
	return(celdat_normalized_exp);
}
map_illuminaHuman2SYMBOL<-function(expd,db){
	as.character(expd[,1])->expd_probes;
	if(db=="illuminaHumanv3SYMBOL"){
		unlist(as.list(illuminaHumanv3SYMBOL))->illuminaHumanv4SYMBOL_list;
	}else if(db=="illuminaHumanv4SYMBOL"){
		unlist(as.list(illuminaHumanv4SYMBOL))->illuminaHumanv4SYMBOL_list;
	}else{
		unlist(as.list(illuminaHumanv2SYMBOL))->illuminaHumanv4SYMBOL_list;
	}
	illuminaHumanv4SYMBOL_list[as.character(expd_probes)]->expd_probes_symbol;
	#select(illuminaHumanv4.db,keys = expd_probes,columns=c("SYMBOL"),keytype="PROBEID")->expd_probes_symbol;
	which(is.na(expd_probes_symbol))->na_index;
	expd_probes[na_index]->expd_probes_symbol[na_index]
	data.frame("Probe"=expd_probes,"gName"=expd_probes_symbol,expd[,-1])->expd_t;
	return(expd_t);
}
map_agilent_probe2SYMBOL<-function(expd,db){
	as.character(expd[,1])->expd_probes;
	if(db=="hgug4112aSYMBOL"){
		unlist(as.list(hgug4112aSYMBOL))->hgug4112aSYMBOL_list;
	}else{
		unlist(as.list(hgug4112aSYMBOL))->hgug4112aSYMBOL_list;
	}
	hgug4112aSYMBOL_list[as.character(expd_probes)]->expd_probes_symbol;
	#select(illuminaHumanv4.db,keys = expd_probes,columns=c("SYMBOL"),keytype="PROBEID")->expd_probes_symbol;
	which(is.na(expd_probes_symbol))->na_index;
	expd_probes[na_index]->expd_probes_symbol[na_index]
	data.frame("Probe"=expd_probes,"gName"=expd_probes_symbol,expd[,-1])->expd_t;
	return(expd_t);
}
map_affy_probe2SYMBOL<-function(expd,db){
	as.character(expd$DATA)->expd_probes;	
	getSYMBOL(expd_probes,db)->expd_probes_substr_SYMBOL;
	which(is.na(expd_probes_substr_SYMBOL))->na_index;
	names(expd_probes_substr_SYMBOL)[na_index]->expd_probes_substr_SYMBOL[na_index]
	data.frame("Probe"=expd[,1],"gName"=expd_probes_substr_SYMBOL,expd[,-1])->expd;
	return(expd)
}
read_geo_data<-function(geof,pattern){
	fread(geof,header=T,sep="\t",stringsAsFactors=F,fill=T)->geof_factors;
	as.data.frame(geof_factors)->geof_factors;#log2 transformed
	grep(pattern,colnames(geof_factors))->probe_index;
	if(length(probe_index)>0){
		geof_factors[,-probe_index]->geof_factors;
	}
	return(geof_factors);
}
preprocess_expd_v2<-function(expd,gStartColumn,aliveEvent,log_trans){
	#-----------remove OS ==NA
	expd[!is.na(expd$A1_OS),]->expd;
	expd[!is.na(expd$A2_Event),]->expd;
	expd[expd$A1_OS!=0,]->expd;
	expd[expd$A1_OS>=30,]->expd;
	#----------remove gene value==1 or 0 in all samples-------
	c()->filter_colums;
	for(i in gStartColumn:ncol(expd)){
		#length(expd[expd[,i]<1,i])->failed_samples;#for TCGA
		length(expd[expd[,i]==0,i])->failed_samples;#for GEO
		if(failed_samples/nrow(expd)<0.5){
			c(filter_colums,i)->filter_colums
		}
	}	
	#----change colnames: '-' to '.'
	gsub("-",".",colnames(expd)[filter_colums])->names(expd)[filter_colums];
	#--------------------------------
	if(log_trans=="yes"){
		data.frame(expd[,1:(gStartColumn-1)],log2(expd[,filter_colums]+1))->expd.filter;
	}else{
		data.frame(expd[,1:(gStartColumn-1)],expd[,filter_colums])->expd.filter;
	}
	print(length(filter_colums));
	flush.console();
	#---------status: 0->alive,1->death---------
	c()->status;
	for(i in 1:nrow(expd.filter)){
		if(expd.filter$A2_Event[i]==aliveEvent){
			c(status,0)->status;
		}else{
			c(status,1)->status; 
		}
	}
	status->expd.filter$Status
	return(expd.filter);
}
#--------------------------------------------------------------------
setwd("D:/TCGA-GBM/hypoxia");
c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(8,"Dark2"),brewer.pal(12,"Set3"),brewer.pal(8,"Accent"),brewer.pal(11,"Spectral")[c(1,4,10,11)])->myd_colors;
#######################################################################
#--------------TCGA GBM affy expression data:
preprocess_expd<-function(expd,gStartColumn,aliveEvent){
	#-----------remove OS ==NA
	expd[!is.na(expd$A1_OS),]->expd;
	expd[!is.na(expd$A2_Event),]->expd;
	expd[expd$A1_OS!=0,]->expd;
	expd[expd$A1_OS>=30,]->expd;
	#----------remove IRGP value==1 or 0 in all samples-------
	c()->filter_colums;
	for(i in gStartColumn:ncol(expd)){
		length(expd[expd[,i]<1,i])->failed_samples;
		if(failed_samples/nrow(expd)<0.5){
			c(filter_colums,i)->filter_colums
		}
	}
	expd[,c(1:(gStartColumn-1),filter_colums)]->expd.filter;
	print(length(filter_colums));
	flush.console();
	#---------status: 0->alive,1->death---------
	c()->status;
	for(i in 1:nrow(expd.filter)){
		if(expd.filter$A2_Event[i]==aliveEvent){
			c(status,0)->status;
		}else{
			c(status,1)->status; 
		}
	}
	status->expd.filter$Status
	return(expd.filter);
}
10->g_start_column
read.table("D:/TCGA-GBM/TME/HT_HG-U133A.merged.factors.txt",header=T,sep="\t",stringsAsFactors=F)->myd;
preprocess_expd(myd,g_start_column,"Alive")->myd_processed;
#--------------------
#-----for GSE13041: Affymetrix Human Genome U133A Array+Affymetrix Human Genome U95 Version 2 Array,#can't use this function;
process_cel_data("D:/TCGA-GBM/GSE13041_family/GSE13041_RAW","rma")->GSE13041_expd;#Affymetrix Human Genome U133A Array
map_affy_probe2SYMBOL(GSE13041_expd,"hgu133a.db")->GSE13041_expd_symbol
write.table(GSE13041_expd_symbol,"GSE13041_expd_symbol.txt",quote=F,sep="\t",row.names=F)
#
read_geo_data("GSE13041_expd_symbol.merged.factors.txt","_at")->GSE13041_expd_factors;
preprocess_expd(GSE13041_expd_factors,9,"ALIVE")->GSE13041_expd_processed;


#---for GSE45301: Illumina HumanHT-12 V3.0 expression beadchip
#read.ilmn("data-set/GSE45301_U87MG_cells_non-normalized.txt",probeid = "ID_REF",expr="SAMPLE",other.columns=c("Detection Pval"),sep="\t")->GSE45301_non_expd;
read.table("data-set/GSE45301_family/MergeExpro_contrib1-GPL6947.txt",header=T,sep="\t",stringsAsFactors=F)->GSE45301_expd;
map_illuminaHuman2SYMBOL(GSE45301_expd,"illuminaHumanv3.db")->GSE45301_expd_symbol;
write.table(GSE45301_expd_symbol,"data-set/GSE45301_expd_symbol.txt",quote=F,sep="\t",row.names=F)
#---for GSE118683: Illumina HumanHT-12 V4.0 expression beadchip
read.table("data-set/GSE118683_family/MergeExpro_contrib1-GPL10558.txt",header=T,sep="\t",stringsAsFactors=F)->GSE118683_expd;
map_illuminaHuman2SYMBOL(GSE118683_expd,"illuminaHumanv4.db")->GSE118683_expd_symbol;
write.table(GSE118683_expd_symbol,"data-set/GSE118683_expd_symbol.txt",quote=F,sep="\t",row.names=F)
################################################################################
#---------------------------------for hypoxia DEGs analysis: use limma packages;
read_geo_data("GSE45301_expd_symbol.merged.factors.txt","ILMN_")->GSE45301_exp_factors;#not log10 transformed, q-quantile normalized
read_geo_data("GSE101799_S16281_normcounts.factors.txt","ILMN_")->GSE101799_exp_factors;#RNA-seq: readcount
read_geo_data("GSE118683_expd_symbol.merged.factors.txt","ILMN_")->GSE118683_exp_factors;#log10 transformed, VSN normalized

################################------limma for DEGs using RMA affy data: GSE45301
GSE45301_exp_factors[grep("U87MG_",GSE45301_exp_factors$Title),]->GSE45301_exp_processed;
as.matrix(GSE45301_exp_processed[,7:ncol(GSE45301_exp_processed)])->geo_raw_matrix;
GSE45301_exp_processed$A0_Samples->rownames(geo_raw_matrix)
t(geo_raw_matrix)->geo_raw_matrix;
#--------------------------------------
GSE45301_exp_processed[,c("A0_Samples","Title")]->GSE45301_condition_group
gsub("U87MG_","",GSE45301_condition_group$Title)->GSE45301_condition_group$Title;
c("SampleID","Condition")->colnames(GSE45301_condition_group);
GSE45301_condition_group->iC_filter_group;
geo_raw_matrix[,as.character(iC_filter_group$SampleID)]->geo_expd_matrix;##
library(parallel)
library(limma)
detectCores()->no_cores;
makeCluster(no_cores-2)->c1;	
clusterExport(c1,c("geo_expd_matrix"));#
parSapply(c1,1:nrow(geo_expd_matrix),function(i){length(geo_expd_matrix[i,geo_expd_matrix[i,]==0])->failed_samples;if(failed_samples/ncol(geo_expd_matrix)<0.5){i}})->filter_colums;
stopCluster(c1);
unlist(filter_colums)->filter_colums;
if(length(filter_colums)>0){
	geo_expd_matrix[filter_colums,]->geo_expd_matrix;
}
#---------------------------------do rank test:
do_ranktest_diff<-function(expd,groups){
	#method: wilcoxon or kruskal wallis
	if(class(expd)!="matrix"){
		print("input data is not matrix!");
		return(NULL);
	}
	names(table(groups$Condition))->condition_table;
	detectCores()->no_cores;
	makeCluster(no_cores-2)->c1;
	c()->test_results;
	if(length(condition_table)>2){
		####
		clusterExport(c1,c("expd","groups","condition_table"),envir=environment());#
		parSapply(c1,rownames(expd),function(g){
			expd[g,as.character(groups$SampleID)]->g_values;
			data.frame("F"=groups$Condition,"V"=g_values)->f_level_df;
			factor(f_level_df[,1],levels=condition_table)->f_level_df[,1];
			paste(c("V","F"),collapse="~")->tmp_formula;
			kruskal.test(as.formula(tmp_formula),data=f_level_df)->g_ks_test;
			c(mean(g_values,na.rm=T),median(g_values,na.rm=T),sd(g_values,na.rm=T),g_ks_test$statistic,g_ks_test$p.value)->tmp_res;
			tmp_res;
		})->test_results;
		unlist(test_results)->test_results;
		matrix(test_results,ncol=5,byrow=T)->test_results_matrix;
		c("MeanExpr","MedianExpr","VarExpr","W","P.value")->colnames(test_results_matrix);
	}else if(length(condition_table)==2){
		as.character(groups[groups$Condition==condition_table[1],1])->group1_samples;
		as.character(groups[groups$Condition==condition_table[2],1])->group2_samples;
		clusterExport(c1,c("expd","group1_samples","group2_samples"),envir=environment());#
		parSapply(c1,rownames(expd),function(g){
			expd[g,group1_samples]->group1_values;
			expd[g,group2_samples]->group2_values;
			mean(group1_values,na.rm=T)->group1_mean;
			mean(group2_values,na.rm=T)->group2_mean;
			wilcox.test(group1_values,group2_values)->group1_group2_test;
			c(group2_mean/group1_mean,group2_mean,group1_group2_test$statistic,group1_group2_test$p.value,group1_mean)->tmp_res;
			tmp_res;
		})->test_results;
		unlist(test_results)->test_results;
		matrix(test_results,ncol=5,byrow=T)->test_results_matrix;
		c("logFC","AveExpr","W","P.value","B")->colnames(test_results_matrix);
	}
	stopCluster(c1);
	data.frame("gName"=rownames(expd),test_results_matrix)->res;
	res[!is.na(res$P.value),]->res;
	res[order(res$P.value),]->res;
	p.adjust(res$P.value)->res$FDR;
	return(res);
}
do_ttest_diff<-function(expd,groups){
	#method: t test for two groups;
	if(class(expd)!="matrix"){
		print("input data is not matrix!");
		return(NULL);
	}
	names(table(groups$Condition))->condition_table;
	detectCores()->no_cores;
	makeCluster(no_cores-2)->c1;
	c()->test_results;
	#-----
	as.character(groups[groups$Condition==condition_table[1],1])->group1_samples;
	as.character(groups[groups$Condition==condition_table[2],1])->group2_samples;
	clusterExport(c1,c("expd","group1_samples","group2_samples"),envir=environment());#
	parSapply(c1,rownames(expd),function(g){
		expd[g,group1_samples]->group1_values;
		expd[g,group2_samples]->group2_values;
		mean(group1_values,na.rm=T)->group1_mean;
		mean(group2_values,na.rm=T)->group2_mean;
		t.test(group1_values,group2_values)->group1_group2_test;
		c(group2_mean/group1_mean,group2_mean,group1_group2_test$statistic,group1_group2_test$p.value,group1_mean)->tmp_res;
		tmp_res;
	})->test_results;
	unlist(test_results)->test_results;
	matrix(test_results,ncol=5,byrow=T)->test_results_matrix;
	c("logFC","AveExpr","t","P.value","B")->colnames(test_results_matrix);
	#------------
	stopCluster(c1);
	data.frame("gName"=rownames(expd),test_results_matrix)->res;
	res[!is.na(res$P.value),]->res;
	res[order(res$P.value),]->res;
	p.adjust(res$P.value)->res$FDR;
	return(res);

}
#do_ranktest_diff(geo_expd_matrix,iC_filter_group)->GSE45301_diff_ranktest;
do_ttest_diff(geo_expd_matrix,iC_filter_group)->GSE45301_diff_ttest;#normoxia/hypoxia
#---------------------------------do limma:
do_limma<-function(expd,colData){
	as.matrix(expd)->expd;
	data.frame("Condition"=colData$Condition)->phenoData;
	colData[,1]->rownames(phenoData);
	if(all(colnames(expd)==rownames(phenoData))){
		new("AnnotatedDataFrame",data = phenoData)->phenoData;
		ExpressionSet(assayData=expd,phenoData = phenoData)->experimentData;
		pData(experimentData)$Condition->condition;
		model.matrix(~condition)->design;
		c("Intercept","Condition_Compare")->colnames(design)
		lmFit(experimentData,design = design)->fit;
		eBayes(fit)->fit_eBayes;
		return(fit_eBayes);
	}else{
		intersect(colnames(expd),rownames(phenoData))->shared_samples;
		phenoData[shared_samples,]->phenoData;
		expd[,shared_samples]->expd;
		print(shared_columns);
		flush.console();
		new("AnnotatedDataFrame",data = phenoData)->phenoData;
		ExpressionSet(assayData=expd,phenoData = phenoData)->experimentData;
		pData(experimentData)$Condition->condition;
		model.matrix(~condition)->design;
		c("Intercept","Condition_Compare")->colnames(design)
		lmFit(experimentData,design = design)->fit;
		eBayes(fit)->fit_eBayes;
		return(fit_eBayes);
	}
}
do_limma(geo_expd_matrix,iC_filter_group)->fit_eBayes;
topTable(fit_eBayes,coef=2,p.value=1,number=nrow(geo_expd_matrix))->GSE45301_ebayes.res;#topTable or topTreat both are ok,
data.frame("gName"=rownames(GSE45301_ebayes.res),GSE45301_ebayes.res)->GSE45301_ebayes.res;#t test does work in this, FDR<0.1,
write.table(GSE45301_ebayes.res,"GSE45301_ebayes.txt",quote=F,sep="\t",row.names=F)

###################################################################################------DESeq2 for DEGs  GSE101799: 
library(DESeq2)
#-----------------------------------------for GSE45301 DEGs: some problems, so this is deprecated
GSE101799_exp_factors[grep("neurosphere",GSE101799_exp_factors$CellType),]->GSE101799_exp_processed;
as.matrix(GSE101799_exp_processed[,6:ncol(GSE101799_exp_processed)])->geo_raw_matrix;
GSE101799_exp_processed$A0_Samples->rownames(geo_raw_matrix)
t(geo_raw_matrix)->geo_raw_matrix;
apply(geo_raw_matrix,2,function(x){round(x,0)})->geo_raw_matrix;
#--------------------------------------
GSE101799_exp_processed[,c("A0_Samples","Treatment")]->GSE101799_condition_group
c("SampleID","Condition")->colnames(GSE101799_condition_group);
GSE101799_condition_group->iC_filter_group;
geo_raw_matrix[,as.character(iC_filter_group$SampleID)]->geo_expd_matrix;##
#--
detectCores()->no_cores;
makeCluster(no_cores-2)->c1;	
clusterExport(c1,c("geo_expd_matrix"));#
parSapply(c1,1:nrow(geo_expd_matrix),function(i){length(geo_expd_matrix[i,geo_expd_matrix[i,]==0])->failed_samples;if(failed_samples/ncol(geo_expd_matrix)<0.5){i}})->filter_colums;
stopCluster(c1);
unlist(filter_colums)->filter_colums;
if(length(filter_colums)>0){
	geo_expd_matrix[filter_colums,]->geo_expd_matrix;
}
#---------------------------------do rank test:
do_ttest_diff(geo_expd_matrix,iC_filter_group)->GSE101799_diff_ttest;#normoxia/hypoxia, not work for readcount;
#---------------------------------do DESeq
do_DESeq<-function(expd,colData){
	DESeqDataSetFromMatrix(countData=expd,colData=colData,design=~Condition)->dds
	keep<-(rowMeans(counts(dds))>=ncol(colData));
	dds[keep,]->dds.keep
	DESeq(dds.keep)->dds.keep.deseq
	results(dds.keep.deseq)->dds_res
	dds_res[order(dds_res$pvalue),]->dds_res_order;
	as.data.frame(dds_res_order)->dds_res_order_df;
	data.frame("TransID"=rownames(dds_res_order_df),dds_res_order_df)->dds_res_order_df;##log2((iC2 exp)/(iC1 exp))
	return(dds_res_order_df);
}
sample_distance_estimate<-function(expd,colData){
	DESeqDataSetFromMatrix(countData=expd,colData=colData,design=~Condition)->dds
	keep<-(rowMeans(counts(dds))>=ncol(colData));
	dds[keep,]->dds.keep
	assay(vst(dds.keep),blind = FALSE)->dds_keep_vst;
	dist(t(dds_keep_vst))->sampleDists
	as.matrix(sampleDists)->sampleDistMatrix
	colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)->colors
	data.frame("Condition"=colData$Condition)->annot_col;
	colData$SampleID->rownames(annot_col)
	brewer.pal(9,"Set1")[c(1,2)]->annot_colors;
	names(table(annot_col))->names(annot_colors)
	list("Condition"=annot_colors)->annot_colors;
	pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors,show_colnames=F,show_rownames=F,annotation_col=annot_col,annotation_row=annot_col,annotation_names_row=F,annotation_names_col=F,annotation_colors=annot_colors)->p;
	return(p);
}
do_DESeq(geo_expd_matrix,iC_filter_group)->GSE101799_deseq.res;
sample_distance_estimate(geo_expd_matrix,iC_filter_group)->GSE101799_expd_matrix_distance;
#print(myd_readcount_matrix_distance)
write.table(GSE101799_deseq.res,"GSE101799_deseq.txt",quote=F,sep="\t",row.names=F)
#-----------------------------------

################################------limma for DEGs using RMA affy data: GSE118683:
GSE118683_exp_factors[grep("IDH1wt",GSE118683_exp_factors$IDH1_Status),]->GSE118683_exp_processed;
as.matrix(GSE118683_exp_processed[,9:ncol(GSE118683_exp_processed)])->geo_raw_matrix;
GSE118683_exp_processed$A0_Samples->rownames(geo_raw_matrix)
t(geo_raw_matrix)->geo_raw_matrix;
#--------------------------------------
GSE118683_exp_processed[,c("A0_Samples","Treatment")]->GSE118683_condition_group
c("SampleID","Condition")->colnames(GSE118683_condition_group);
GSE118683_condition_group->iC_filter_group;
geo_raw_matrix[,as.character(iC_filter_group$SampleID)]->geo_expd_matrix;##
#------
detectCores()->no_cores;
makeCluster(no_cores-2)->c1;	
clusterExport(c1,c("geo_expd_matrix"));#
parSapply(c1,1:nrow(geo_expd_matrix),function(i){length(geo_expd_matrix[i,geo_expd_matrix[i,]==0])->failed_samples;if(failed_samples/ncol(geo_expd_matrix)<0.5){i}})->filter_colums;
stopCluster(c1);
unlist(filter_colums)->filter_colums;
if(length(filter_colums)>0){
	geo_expd_matrix[filter_colums,]->geo_expd_matrix;
}
#---------------------------------do rank test:
do_ttest_diff(geo_expd_matrix,iC_filter_group)->GSE118683_diff_ttest;#normoxia/hypoxia
#---------------------------------do limma:
do_limma(geo_expd_matrix,iC_filter_group)->fit_eBayes;
topTable(fit_eBayes,coef=2,p.value=1,number=nrow(geo_expd_matrix))->GSE118683_ebayes.res;#topTable or topTreat both are ok,
data.frame("gName"=rownames(GSE118683_ebayes.res),GSE118683_ebayes.res)->GSE118683_ebayes.res;#t test does work in this, FDR<0.1,
write.table(GSE118683_ebayes.res,"GSE118683_ebayes.txt",quote=F,sep="\t",row.names=F)

########################################################################
#------------------------compare: 
#--------------------draw volcano figures for DL and LMS, UPS and MFS:
draw_volcano_deseq_figure<-function(deseq,cut_p,cut_lfc,p_title){
	brewer.pal(9,"Set1")->myd_colors;
	1->deseq[is.na(deseq$padj),"padj"];
	c()->filter_index;
	c()->deseq_color;
	for(i in 1:nrow(deseq)){
		if(deseq$log2FoldChange[i]<(-cut_lfc) && deseq$padj[i]<cut_p){
			c(filter_index,i)->filter_index;
			c(deseq_color,myd_colors[2])->deseq_color;
		}else if(deseq$log2FoldChange[i]>cut_lfc && deseq$padj[i]<cut_p){
			c(filter_index,i)->filter_index;
			c(deseq_color,myd_colors[1])->deseq_color;
		}else{
			c(deseq_color,myd_colors[9])->deseq_color
		}
	}
	#print(table(deseq_color));flush.console();
	plot(x=deseq$log2FoldChange,y=-log(deseq$padj),pch=20,cex=0.5,col=deseq_color,xlim=c(-5,5),xlab="Log2 fold change",ylab="-Log10(padj)")
	legend("topright",legend=c("Down","Up"),pch=20,cex=1,col=myd_colors[c(2,1)],title=p_title,border=NA)
	return(deseq[filter_index,]);
}
draw_volcano_limma_figure<-function(resd,cut_p,cut_lfc,p_title){
	brewer.pal(9,"Set1")->myd_colors;
	c()->filter_index;
	c()->deseq_color;
	for(i in 1:nrow(resd)){
		if(resd$logFC[i]<(-cut_lfc) && resd$adj.P.Val[i]<cut_p){
			c(filter_index,i)->filter_index;
			c(deseq_color,myd_colors[2])->deseq_color;
		}else if(resd$logFC[i]>cut_lfc && resd$adj.P.Val[i]<cut_p){
			c(filter_index,i)->filter_index;
			c(deseq_color,myd_colors[1])->deseq_color;
		}else{
			c(deseq_color,myd_colors[9])->deseq_color
		}
	}
	#print(table(deseq_color));flush.console();
	plot(x=resd$logFC,y=-log(resd$adj.P.Val),pch=20,cex=0.5,col=deseq_color,xlab="Log2 fold change",ylab="-Log10(padj)")
	legend("topright",legend=c("Down","Up"),pch=20,cex=1,col=myd_colors[c(2,1)],title=p_title,border=NA)
	return(resd[filter_index,]);
}
par(mfrow=c(1,2))
draw_volcano_limma_figure(GSE45301_ebayes.res,0.2,0,"normoxia/hypoxia")->GSE45301_ebayes.res_filter;
mtext("GSE45301")
draw_volcano_limma_figure(GSE118683_ebayes.res,0.2,0,"normoxia/hypoxia")->GSE118683_ebayes.res_filter
mtext("GSE118683")

######################################################################################
draw_survial_curve_custom<-function(myd,column,bk,myd_colors){
	myd[myd[,column]!="",]->myd.rm;
	myd.rm[!is.na(myd.rm[,column]),]->myd.rm;
	if(length(myd.rm$A1_OS)>1){
		Surv(as.numeric(myd.rm$A1_OS),as.numeric(myd.rm$Status))->myd.surv;
	}else if(length(myd.rm$OS)>1){
		Surv(as.numeric(myd.rm$OS),as.numeric(myd.rm$Status))->myd.surv;
	}
	survfit(formula=myd.surv~myd.rm[,column])->myd.fit;
	survdiff(formula=myd.surv~myd.rm[,column],rho=0)->myd.diff;
	table(myd.rm[,column])->myd.table;
	max(myd.rm$A1_OS)+100->max_xlim;
	plot(myd.fit,col=myd_colors,xlab="Time(days)",ylab="Overall Survival(%)",lwd=2,axes=F,main=paste("KM(",colnames(myd.rm)[column],")",sep=""),xlim=c(-100,max_xlim*1.1));
	axis(side=1,at=seq(0,max_xlim,bk),labels=seq(0,max_xlim,bk),pos=0);
	rug(x=seq(0,max_xlim,bk)+bk/2,ticksize=-0.01,side=1,pos=0);
	axis(side=2,at=seq(0,1,0.2),labels=seq(0,100,20),pos=0);
	rug(x=seq(0,0.9,0.2)+0.1,ticksize=-0.01,side=2,pos=0);
	#abline(h=seq(0.2,1,0.2),col=brewer.pal(9,"Greys")[3],lty=3)
	1-pchisq(myd.diff$chisq,df=length(myd.diff$n)-1)->pvalue;
	legend("topright",legend=paste(names(myd.table),paste("(N=",myd.table,")",sep="")),fill=myd_colors,bty="n",cex=1.2);
	if(pvalue<1e-5){
		legend(x=bk/2,y=0.15,legend="p < 1e-5",bty="n",cex=1.2)
	}else{
		legend(x=bk/2,y=0.15,legend=paste("p=",round(pvalue,5),sep=""),bty="n",cex=1.2)
	}
	return(c(pvalue,myd.table));
}
cox_univariant_gene_regr<-function(myd,colums){
	Surv(as.numeric(myd$A1_OS),as.numeric(myd$Status))->myd.surv
	c()->univar_anova_p;
	c()->univar_coxph_HR;
	c()->univar_coxph_low95;
	c()->univar_coxph_high95;
	c()->univar_coxph_logtest;
	c()->univar_coxph_sctest;
	c()->univar_coxph_waldtest;
	c()->fpkm.mean;
	c()->fpkm.median;
	colnames(myd)[colums]->myd.names;
	for(i in myd.names){
		mean(myd[,i],na.rm=T)->tmp.mean;
		median(myd[,i],na.rm=T)->tmp.median;
		c(fpkm.mean,tmp.mean)->fpkm.mean;
		c(fpkm.median,tmp.median)->fpkm.median;
		as.formula(paste("myd.surv~",i))->tmp.formula;
		coxph(formula=tmp.formula,data=myd)->tmp.coxph;
		summary(tmp.coxph)->tmp.coxph.summary;
		c(univar_anova_p,tmp.coxph.summary$coefficients[,5])->univar_anova_p;
		c(univar_coxph_HR,tmp.coxph.summary$coefficients[,2])->univar_coxph_HR;
		c(univar_coxph_low95,tmp.coxph.summary$conf.int[,3])->univar_coxph_low95;
		c(univar_coxph_high95,tmp.coxph.summary$conf.int[,4])->univar_coxph_high95;
		c(univar_coxph_logtest,tmp.coxph.summary$logtest[3])->univar_coxph_logtest;
		c(univar_coxph_sctest,tmp.coxph.summary$sctest[3])->univar_coxph_sctest;
		c(univar_coxph_waldtest,tmp.coxph.summary$waldtest[3])->univar_coxph_waldtest;	
	}
	data.frame("gName"=myd.names,"Pvalue"=univar_anova_p,"HR"=univar_coxph_HR,"Low(95%CI)"=univar_coxph_low95,"High(95%CI)"=univar_coxph_high95,"Logrank"=univar_coxph_logtest,"Sctest"=univar_coxph_sctest,"Waldtest"=univar_coxph_waldtest,"fpkm_median"=fpkm.median,"fpkm_mean"=fpkm.mean)->myd.coxph.df;
	#myd.coxph.df[myd.coxph.df$fpkm_median>=1,]->myd.coxph.df;
	myd.coxph.df[!is.na(myd.coxph.df$Logrank),]->myd.coxph.df;
	myd.coxph.df[order(myd.coxph.df$Logrank),]->myd.coxph.df;
	return(myd.coxph.df);
}
#---------------------------------for hypoxia related genes from reactome+kegg+go:
readLines("merged.hypoxia.genes")->merged_hypoxia.genes;
#-------------------------compare hypoxia genes: 
data.frame("Source"=rep(c("General","GSE45301","GSE118683"),c(length(merged_hypoxia.genes),nrow(GSE45301_ebayes.res_filter),nrow(GSE118683_ebayes.res_filter))),"gName"=c(merged_hypoxia.genes,as.character(GSE45301_ebayes.res_filter$gName),as.character(GSE118683_ebayes.res_filter$gName)))->diff_genes.df
table(diff_genes.df$gName,diff_genes.df$Source)->diff_genes.df_table;
vennCounts(diff_genes.df_table)->deseq_list.vennCount;
vennDiagram(deseq_list.vennCount,names=c("General","GSE118683","GSE45301"),cex=1.5,lwd=1.5,circle.col=brewer.pal(9,"Set1"));
######################################################################################
#---for GEO and KEGG enrichment:
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(gridExtra);
library(enrichplot);#upsetplot,enrichMap... #http://www.webgestalt.org
draw_GO_KEGG_plot<-function(myd,showCategory){
	if(length(myd)>0){
		c(brewer.pal(9,"Set1"),brewer.pal(11,"Set3"),brewer.pal(11,"Spectral"))->myd.colors;
		ifelse(nrow(myd)<showCategory,nrow(myd),showCategory)->showCategory;
		myd[1:showCategory,]->myd;
		myd$Count->myd.count;
		myd.count->myd_count_color;
		if(max(myd.count)>30){
			round(myd.count/2)->myd_count_color;
		}
		-log10(myd$pvalue)->myd.pvalue;
		plot(x=myd.count,y=myd.pvalue,pch=19,cex=-log10(myd$pvalue),col=myd.colors[myd_count_color],xlab="Count",ylab="-log10(P value)",xlim=c(min(myd.count)-1,max(myd.count)+1),ylim=c(min(myd.pvalue)-0.2,max(myd.pvalue)+0.2));
		text(x=myd.count,y=myd.pvalue,labels=myd$Description,offset=0.5,pos=3,cex=0.8,adj=c(0,1))
	}else{
		plot.new();
	}
}
do_GO_KEGG<-function(genes,gtitle){
	as.character(genes)->GO_KEGG.selected_genes;
	groupGO(gene=GO_KEGG.selected_genes,OrgDb=org.Hs.eg.db,ont="BP",keyType='SYMBOL')->GO_KEGG.selected_genes.group;
	as.data.frame(GO_KEGG.selected_genes.group)->GO_KEGG.selected_genes.group;
	write.table(GO_KEGG.selected_genes.group,paste(gtitle,".groupGO.txt",sep=""),quote=F,row.names=F,sep="\t");
	enrichGO(gene=GO_KEGG.selected_genes,OrgDb=org.Hs.eg.db,keyType='SYMBOL',ont="BP",pAdjustMethod="BH",pvalueCutoff=0.1)->GO_KEGG.selected_genes.GO;
	GO_KEGG.selected_genes.GO->GO.enrichment;
	#simplify(GO_KEGG.selected_genes.GO)->GO_KEGG.selected_genes.GO;
	as.data.frame(GO_KEGG.selected_genes.GO)->GO_KEGG.selected_genes.GO
	write.table(GO_KEGG.selected_genes.GO,paste(gtitle,".GOEnrichment.txt",sep=""),quote=F,row.names=F,sep="\t")
	bitr(GO_KEGG.selected_genes,fromType="SYMBOL",toType=c("ENSEMBL","ENTREZID"),OrgDb=org.Hs.eg.db)->myd.gName.df;
	enrichKEGG(gene=myd.gName.df[,3],organism="hsa",keyType="ncbi-geneid",pvalueCutoff=0.5)->myd.gName.df.KEGG
	write.table(as.data.frame(myd.gName.df.KEGG),paste(gtitle,".KEGGEnrichment.txt",sep=""),quote=F,row.names=F,sep="\t")
	#if(nrow(GO_KEGG.selected_genes.GO)>0 && nrow(myd.gName.df.KEGG)>0){
	#	par(mfrow=c(1,2))
	#	draw_GO_KEGG_plot(GO_KEGG.selected_genes.GO,10)
	#	draw_GO_KEGG_plot(myd.gName.df.KEGG,10)
	#}else if(nrow(GO_KEGG.selected_genes.GO)>0 && nrow(myd.gName.df.KEGG)==0){
	#	par(mfrow=c(1,1))
	#	draw_GO_KEGG_plot(GO_KEGG.selected_genes.GO,10)
	#}else if(nrow(GO_KEGG.selected_genes.GO)==0 && nrow(myd.gName.df.KEGG)>0){
	#	par(mfrow=c(1,1))
	#	draw_GO_KEGG_plot(myd.gName.df.KEGG,10)
	#}
	return(list("enrichGO"=GO.enrichment,"enrichKEGG"=myd.gName.df.KEGG));
}
unique(GSE45301_ebayes.res_filter$gName,GSE118683_ebayes.res_filter$gName)->geo_ebayes.genes;
intersect(merged_hypoxia.genes,geo_ebayes.genes)->geo_ebayes.shared_genes;
do_GO_KEGG(geo_ebayes.shared_genes,"geo_ebayes.shared_genes")->geo_ebayes.shared_genes.enrichment;
dotplot(geo_ebayes.shared_genes.enrichment$enrichGO)->GO.p1
dotplot(geo_ebayes.shared_genes.enrichment$enrichKEGG)->KEGG.p2
grid.arrange(p1,p2,ncol=2)

#--------------coxph for TCGA:
which(colnames(myd_processed)%in%geo_ebayes.shared_genes)->geo_ebayes.shared_genes_index;
cox_univariant_gene_regr(myd_processed,geo_ebayes.shared_genes_index)->myd_processed.coxph;
myd_processed.coxph[!is.na(myd_processed.coxph$Pvalue),]->myd_processed.coxph
#--------------coxph for GSE13041:
which(colnames(GSE13041_expd_processed)%in%geo_ebayes.shared_genes)->geo_ebayes.shared_genes_index;
cox_univariant_gene_regr(GSE13041_expd_processed,geo_ebayes.shared_genes_index)->GSE13041_expd_processed.coxph;
GSE13041_expd_processed.coxph[!is.na(GSE13041_expd_processed.coxph$Pvalue),]->GSE13041_expd_processed.coxph
#---
prepare_superheatmap<-function(coxphd,cutp){
	c()->datasets;
	list()->genes;
	list()->pvalues;
	c()->x_length;
	for(x in coxphd){
		x[x$Pvalue<cutp,"gName"]->x_genes;
		x[x$Pvalue<cutp,"HR"]->x_pvalues;
		c(x_length,length(x_genes))->x_length;
		c(genes,list(x_genes))->genes;
		c(pvalues,list(x_pvalues))->pvalues;
	}
	as.character(unique(unlist(genes)))->uniq_genes;
	names(coxphd)->coxphd_names;
	c()->uniq_pvalues;
	coxphd_names->names(genes);coxphd_names->names(pvalues);
	for(i in coxphd_names){
		genes[[i]]->i_genes;
		pvalues[[i]]->i_pvalues;
		setdiff(uniq_genes,i_genes)->i_diff_genes;
		c(i_pvalues,rep(0,length(i_diff_genes)))->i_pvalues;
		c(as.character(i_genes),i_diff_genes)->names(i_pvalues);
		i_pvalues[uniq_genes]->i_pvalues;
		c(uniq_pvalues,i_pvalues)->uniq_pvalues;
	}
	matrix(uniq_pvalues,ncol=length(uniq_genes),byrow=T)->res;
	uniq_genes->colnames(res);
	coxphd_names->rownames(res);
	return(res);
}
#------for TCGA+GSE13041 shared DEGs and coxph filter:
list("TCGA"=myd_processed.coxph,"GSE13041"=GSE13041_expd_processed.coxph)->coxph_list;
prepare_superheatmap(coxph_list,0.05)->coxph_list.merged_df;
unlist(lapply(1:ncol(coxph_list.merged_df),function(x){
	which(coxph_list.merged_df[,x]!=0)->x0;
	which(coxph_list.merged_df[,x]>1)->x1;
	which(coxph_list.merged_df[,x]<1)->x2;
	if(length(x0)>1){
		if(length(x1)>1 || length(x2)>1){
			x;
		}
	}
}))->coxph_list.index;
colnames(coxph_list.merged_df)[coxph_list.index]->TCGA_GSE13041_shared_coxph_genes;

#----------------------use GBM specific hypoxia related genes to cluster: 
#-------------------
prepare_pheatmap_cluster<-function(pheatmapd,cuts,margin){
	NULL->cuts_d;
	if(margin==1){
		#pheatmapd$tree_row$labels->row_names;
		cutree(pheatmapd$tree_row,cuts)->cuts_d;
	}else if(margin==2){
		cutree(pheatmapd$tree_col,cuts)->cuts_d;
	}
	if(!is.null(cuts_d)){
		data.frame("Sample"=names(cuts_d),"KCluster"=paste("KC",cuts_d,sep=""))->res;
		return(res);
	}else{
		return(NULL);
	}
}
intersect(colnames(tmp_exp),geo_ebayes.shared_genes)->selected_geneSet;#
c(colorRampPalette(rev(brewer.pal(11,"RdYlBu")[7:11]))(200),colorRampPalette(rev(brewer.pal(11,"RdYlBu")[1:6]))(250))->fill_colors;
tmp_exp[,as.character(selected_geneSet)]->tmp_exp.matrix;
tmp_exp$A0_Samples->rownames(tmp_exp.matrix);
pheatmap(t(tmp_exp.matrix),show_colnames=F,show_rownames=F,color=fill_colors)->pheatmap.p1
prepare_pheatmap_cluster(pheatmap.p1,3,2)->tmp_exp.kcluster;
"KC2"->tmp_exp.kcluster[tmp_exp.kcluster$KCluster=="KC3","KCluster"]
merge(tmp_exp,tmp_exp.kcluster,by.x="A0_Samples",by.y="Sample")->tmp_expd_cluster;
par(mfrow=c(1,2));
draw_survial_curve_custom(tmp_expd_cluster,ncol(tmp_expd_cluster),500,myd_colors)#p=0.0134
#---
GSE13041_expd_processed[,geo_ebayes.shared_genes_index]->GSE13041_expd_processed.matrix;
GSE13041_expd_processed$A0_Samples->rownames(GSE13041_expd_processed.matrix)
pheatmap(t(GSE13041_expd_processed.matrix),show_colnames=F,show_rownames=F,color=fill_colors)->pheatmap.p2
prepare_pheatmap_cluster(pheatmap.p2,2,2)->GSE13041_expd_processed.kcluster;
merge(GSE13041_expd_processed,GSE13041_expd_processed.kcluster,by.x="A0_Samples",by.y="Sample")->GSE13041_expd_processed.kcluster;
draw_survial_curve_custom(GSE13041_expd_processed.kcluster,ncol(GSE13041_expd_processed.kcluster),500,myd_colors)#p=0.13
#--
grid.arrange(pheatmap.p1$gtable,pheatmap.p2$gtable,nrow=2)
#-----use mclust to select the the best clusters: 
library(mclust)
Mclust(tmp_exp.matrix,G=2:50)->tmp_exp.mclust_d;#for Mclust: row correspond to observations(n), and columns correspond to variables(m)
plot(tmp_exp.mclust_d$BIC)#the best gene clusters is 2: ,kmeans_k=2
#
Mclust(GSE13041_expd_processed.matrix,G=2:50)->GSE13041_expd_processed.mclust_d;#for Mclust: row correspond to observations(n), and columns correspond to variables(m)
plot(GSE13041_expd_processed.mclust_d$BIC)#the best gene clusters is 2: ,kmeans_k=2
#--for CGGA_301array_expd:
which(colnames(CGGA_301array_expd_processed)%in%geo_ebayes.shared_genes)->CGGA_301array.shared_genes_index;
CGGA_301array_expd_processed[,CGGA_301array.shared_genes_index]->CGGA_301array_expd_processed.matrix;
CGGA_301array_expd_processed$A0_Samples->rownames(CGGA_301array_expd_processed.matrix)
pheatmap(t(CGGA_301array_expd_processed.matrix),show_colnames=F,show_rownames=F,color=fill_colors)->pheatmap.p3
prepare_pheatmap_cluster(pheatmap.p3,2,2)->CGGA_301array_expd_processed.kcluster;
merge(CGGA_301array_expd_processed,CGGA_301array_expd_processed.kcluster,by.x="A0_Samples",by.y="Sample")->CGGA_301array_expd_processed.kcluster;
draw_survial_curve_custom(CGGA_301array_expd_processed.kcluster,ncol(CGGA_301array_expd_processed.kcluster),500,myd_colors)#p=0.13


#--------------------------------------------------------------------------use IRG method: Result1
calculate_cluster_centroids<-function(expd,expd_cluster,colCentroid){
	unique(expd_cluster$KCluster)->clusters;
	if(colCentroid){
		t(expd)->expd;
	}
	c()->centroid_values;
	for(cl in clusters){
		which(expd_cluster$KCluster==cl)->cl_index;
		as.character(expd_cluster$Sample[cl_index])->cl_samples;
		expd[cl_samples,]->expd_cl;
		#print(cl_samples);flush.console();
		apply(expd_cl,2,function(x){mean(x)})->means;
		c(centroid_values,means)->centroid_values;
	}
	matrix(centroid_values,nrow=length(clusters),byrow=T)->centroid_values.matrix;
	clusters->rownames(centroid_values.matrix);
	colnames(expd)->colnames(centroid_values.matrix);
	return(centroid_values.matrix)
}
#--------------calculate GSE13041 IGP values;
calculate_cluster_centroids(tmp_exp.matrix,tmp_exp.kcluster,FALSE)->tmp_expd_cluster.centroids;
#------------
library(clusterRepro)
t(GSE13041_expd_processed.matrix)->GSE13041_expd_processed.matrix;
t(tmp_expd_cluster.centroids)->tmp_expd_cluster.centroids;
clusterRepro(tmp_expd_cluster.centroids, GSE13041_expd_processed.matrix, Number.of.permutations = 1000)->Result1;
#--------------
IGP.clusterRepro(Data=GSE13041_expd_processed.matrix,tmp_expd_cluster.centroids)->GSE13041_expd_processed.result;
data.frame("Sample"=names(GSE13041_expd_processed.result$Class),"IRGCluster"=paste("C",GSE13041_expd_processed.result$Class,sep=""))->GSE13041_expd_processed.clusterRepro;
merge(GSE13041_expd_processed,GSE13041_expd_processed.clusterRepro,by.x="A0_Samples",by.y="Sample")->GSE13041_expd_processed.IGPcluster;
#draw_survial_curve_custom(GSE13041_expd_processed.IGPcluster,ncol(GSE13041_expd_processed.IGPcluster),500,myd_colors)#p<0.374
#----------draw TCGA_GSE13041_shared_coxph_genes boxplot in TCGA and GSE13041:
#--------------------------plot gene expression in boxplot
draw_boxplot_genes_by_factors<-function(expd,genes,f,myd_colors){
	#-----remove NA:
	expd[which(expd[,f]!="Un"),]->expd;
	#-
	expd[,f]->f_values;
	as.character(genes)->genes;
	intersect(genes,colnames(expd))->genes;
	table(expd[,f])->f_table;
	names(f_table)->f_names;
	paste(f_names,f_table,sep="(")->f_names_legend;
	paste(f_names_legend,")",sep="")->f_names_legend;
	lapply(f_names,function(x){which(f_values==x)->res;res})->f_names_index;
	lapply(genes,function(x){expd[,x]->x_values;lapply(f_names_index,function(y){x_values[y]->res;res})})->genes_f_values;
	list()->genes_f_values_list;
	c()->rank_p_values;
	for(i in 1:length(genes)){
		c()->f_levels;
		c()->f_level_values;
		for(j in 1:length(f_names)){
			c(f_levels,rep(f_names[j],length(genes_f_values[[i]][[j]])))->f_levels;
			c(f_level_values,genes_f_values[[i]][[j]])->f_level_values;
			c(genes_f_values_list,list(genes_f_values[[i]][[j]]))->genes_f_values_list
		}
		data.frame("F"=f_levels,"V"=f_level_values)->f_level_df;
		factor(f_level_df[,1],levels=f_names)->f_level_df[,1];
		paste(c("V","F"),collapse="~")->tmp_formula;
		c(rank_p_values,kruskal.test(as.formula(tmp_formula),data=f_level_df)$p.value)->rank_p_values;
	}
	c()->sig_symbol;
	for(x in rank_p_values){
		if(x <= 1e-5){
			c(sig_symbol,"***")->sig_symbol
		}else if(x > 1e-5 && x <= 0.01){
			c(sig_symbol,"**")->sig_symbol
		}else if(x > 0.01 && x <= 0.05){
			c(sig_symbol,"*")->sig_symbol
		}else{
			c(sig_symbol,"")->sig_symbol
		}
	}
	data.frame("Factors"=genes,"KruskalP"=rank_p_values)->kruskal_test_res;
	seq(length(f_names)+1,(length(genes)+1)*(length(f_names)+1)-1,length(f_names)+1)->seg_index;
	seg_index-1->seg_end_index;
	seg_index-length(f_names)->seg_start_index;
	c()->x_at;
	c()->x_label_at;
	for(i in 1:length(seg_start_index)){
		seq(seg_start_index[i],seg_end_index[i])->res;
		c(x_label_at,mean(res))->x_label_at;
		c(x_at,res)->x_at;
	}
	boxplot(genes_f_values_list,plot=F)->genes_f_values.boxplot;
	min(genes_f_values.boxplot$stats)->boxplot_min;
	max(genes_f_values.boxplot$stats)->boxplot_max;
	boxplot(genes_f_values_list,at=x_at,axes=F,ylab="Values",boxwex=0.5,pch=20,cex=0.5,col=myd_colors[1:length(f_names)],ylim=c(round(boxplot_min),boxplot_max)*1.5)#->genes_f_values.boxplot;
	(boxplot_max-round(boxplot_min))/5->bk;
	axis(side=1,at=x_label_at,labels=genes,tick=T,las=2)
	if(is.wholenumber(boxplot_max) && abs(boxplot_max) < 1){
		axis(side=2,at=seq(round(boxplot_min),boxplot_max*1.1,bk),labels=seq(round(boxplot_min),boxplot_max*1.1,bk))
	}else{
		axis(side=2,at=seq(round(boxplot_min),boxplot_max*1.1,bk),labels=round(seq(round(boxplot_min),boxplot_max*1.1,bk),3))
	}
	text(x=x_label_at,y=boxplot_max*1.3,labels=sig_symbol,cex=1.3)
	print(bk);
	flush.console();
	#legend(x=mean(seg_index),y=round(boxplot_max)*1.4,horiz=T,fill=myd_colors[1:length(f_names)],legend=f_names,title=f,inset=0,xjust=0.5,yjust=0.5)
	legend("top",horiz=T,fill=myd_colors[1:length(f_names)],legend=f_names_legend,title=f,inset=0,xjust=0.5,yjust=0.5)
	legend("topright",horiz=T,legend=c("*** p<1e-5","** p<0.01","* p<0.05"),title="Kruskal-Wallis Test")
	box();
	return(kruskal_test_res);
}
draw_boxplot_genes_by_factors_v2<-function(expd,genes,f,myd_colors,log_trans){
	#-----remove NA:
	expd[which(expd[,f]!="Un"),]->expd;
	#-
	expd[,f]->f_values;
	as.character(genes)->genes;
	intersect(genes,colnames(expd))->genes;
	table(expd[,f])->f_table;
	names(f_table)->f_names;
	paste(f_names,f_table,sep="(")->f_names_legend;
	paste(f_names_legend,")",sep="")->f_names_legend;
	lapply(f_names,function(x){which(f_values==x)->res;res})->f_names_index;
	if(log_trans!="YES"){
		lapply(genes,function(x){expd[,x]->x_values;lapply(f_names_index,function(y){x_values[y]->res;res})})->genes_f_values;
	}else{
		lapply(genes,function(x){log2(expd[,x]+1)->x_values;lapply(f_names_index,function(y){x_values[y]->res;res})})->genes_f_values;
	}
	list()->genes_f_values_list;
	c()->rank_p_values;
	for(i in 1:length(genes)){
		c()->f_levels;
		c()->f_level_values;
		for(j in 1:length(f_names)){
			c(f_levels,rep(f_names[j],length(genes_f_values[[i]][[j]])))->f_levels;
			c(f_level_values,genes_f_values[[i]][[j]])->f_level_values;
			c(genes_f_values_list,list(genes_f_values[[i]][[j]]))->genes_f_values_list
		}
		data.frame("F"=f_levels,"V"=f_level_values)->f_level_df;
		factor(f_level_df[,1],levels=f_names)->f_level_df[,1];
		paste(c("V","F"),collapse="~")->tmp_formula;
		c(rank_p_values,kruskal.test(as.formula(tmp_formula),data=f_level_df)$p.value)->rank_p_values;
	}
	c()->sig_symbol;
	for(x in rank_p_values){
		if(x <= 1e-5){
			c(sig_symbol,"***")->sig_symbol
		}else if(x > 1e-5 && x <= 0.01){
			c(sig_symbol,"**")->sig_symbol
		}else if(x > 0.01 && x <= 0.05){
			c(sig_symbol,"*")->sig_symbol
		}else{
			c(sig_symbol,"")->sig_symbol
		}
	}
	data.frame("Factors"=genes,"KruskalP"=rank_p_values)->kruskal_test_res;
	seq(length(f_names)+1,(length(genes)+1)*(length(f_names)+1)-1,length(f_names)+1)->seg_index;
	seg_index-1->seg_end_index;
	seg_index-length(f_names)->seg_start_index;
	c()->x_at;
	c()->x_label_at;
	for(i in 1:length(seg_start_index)){
		seq(seg_start_index[i],seg_end_index[i])->res;
		c(x_label_at,mean(res))->x_label_at;
		c(x_at,res)->x_at;
	}
	boxplot(genes_f_values_list,plot=F)->genes_f_values.boxplot;
	min(genes_f_values.boxplot$stats)->boxplot_min;
	max(genes_f_values.boxplot$stats)->boxplot_max;
	boxplot(genes_f_values_list,at=x_at,axes=F,ylab="Values",lty=1,boxcol=myd_colors[1:length(f_names)],boxwex=0.5,pch=20,cex=0.5,col=myd_colors[1:length(f_names)],ylim=c(round(boxplot_min),round(boxplot_max)*1.5))#->genes_f_values.boxplot;
	(round(boxplot_max)-round(boxplot_min))/5->bk;
	axis(side=1,at=x_label_at,labels=genes,tick=T,las=2)
	axis(side=2,at=seq(round(boxplot_min),round(boxplot_max)*1.1,bk),labels=seq(round(boxplot_min),round(boxplot_max)*1.1,bk))
	#legend(x=mean(seg_index),y=round(boxplot_max)*1.2,horiz=T,fill=myd_colors[1:length(f_names)],legend=f_names,title=f,inset=0,xjust=0.5,yjust=0.5)
	text(x=x_label_at,y=boxplot_max*1.3,labels=sig_symbol,cex=1.3)
	legend("top",horiz=T,fill=myd_colors[1:length(f_names)],legend=f_names_legend,title=f,inset=0,xjust=0.5,yjust=0.5)
	legend("topright",horiz=T,legend=c("*** p<1e-5","** p<0.01","* p<0.05"),title="Kruskal-Wallis Test")
	#abline(v=seg_index,lty=3)
	box();
	return(kruskal_test_res);
}
#------------
draw_boxplot_factors_by_gene_v1<-function(expd,clusters,f,f_column,myd_colors){
	which(expd[,f_column]==f)->f_index;
	expd[f_index,]->expd;
	list()->f_values;
	for(cl in clusters){
		c(f_values,list(log2(expd[,cl]+1)))->f_values;
	}
	clusters->names(f_values)
	data.frame("F"=rep(clusters,each=length(f_index)),"V"=unlist(f_values))->f_level_df;
	factor(f_level_df[,1],levels=clusters)->f_level_df[,1];
	paste(c("V","F"),collapse="~")->tmp_formula;
	kruskal.test(as.formula(tmp_formula),data=f_level_df)$p.value->rank_p_values;
	#------------------
	rank_p_values->x;
	if(x <= 1e-3){
		"p < 0.001"->sig_symbol
	}else{
		round(x,4)->sig_symbol
	}
	#-----------
	paste(f," Count: ",length(f_index),sep="")->f_count;
	paste("Rank Test ",sig_symbol,sep="")->f_pvalue;
	boxplot(f_values,col=myd_colors[1:length(clusters)],boxwex=0.4,ylab="log2(FPKM)",pch=20,cex=0.7,xaxt="n")
	#axis(1,at=seq(1,length(clusters),1),labels=clusters,tick=F,cex=0.3,las=2)
	legend("top",legend=clusters,fill=myd_colors[1:length(clusters)],bty='n',horiz=T)
	mtext(text=paste(f_count,f_pvalue,sep=","),font=3);
}
par(mar=c(5,4,4,4),mfrow=c(1,2))
draw_boxplot_genes_by_factors_v2(tmp_expd_cluster,TCGA_GSE13041_shared_coxph_genes,"KCluster",myd_colors,"NO")->x
draw_boxplot_genes_by_factors_v2(GSE13041_expd_processed.kcluster,TCGA_GSE13041_shared_coxph_genes,"KCluster",myd_colors,"NO")->x
draw_boxplot_genes_by_factors_v2(CGGA_301array_expd_processed.kcluster,TCGA_GSE13041_shared_coxph_genes,"KCluster",myd_colors,"NO")->x

##########################################################################################################################################
library(NMF)
library(iClusterPlus)
library(bnstruct)
##########################################################################################################################################
########read Methylation file: 27k and 450k 
#"GBM_450k_0.7cutoff.rmCrossMap.factors.res"->infile
read_methy_data<-function(infile,myd_samples){
	read.table(infile,header=T,sep="\t",stringsAsFactors=F)->myd.methy;
	#knn impute method to handle NA values
	t(myd.methy[,10:ncol(myd.methy)])->myd.methy_t;
	knn.impute(myd.methy_t,k=10,cat.var=1:ncol(myd.methy_t),to.impute=1:nrow(myd.methy_t),using=1:nrow(myd.methy_t))->myd.methy_t
	data.frame(myd.methy[,c(1:9)],t(myd.methy_t))->myd.methy;
	#-----------
	which(gsub("\\.","-",colnames(myd.methy))%in%myd_samples)->methy_index;
	myd.methy[,c(1:9,methy_index)]->myd.methy
	#write.table(myd.methy,"PAAD_450k_07.methy_filter.txt",quote=F,sep="\t",row.names=F)
	return(myd.methy)
}
paste(tmp_expd_cluster$A0_Samples,"-01",sep="")->myd_samples;
myd_samples->rownames(tmp_expd_cluster);
read_methy_data("GBM_450k_0.7cutoff.rmCrossMap.factors.res",myd_samples)->myd.methy_450k;
read_methy_data("GBM_27k_0.7cutoff.rmCrossMap.factors.res",myd_samples)->myd.methy_27k;
merge_two_methy_data<-function(methyd1,methyd2,genes){
	which(methyd1$Symbol%in%genes)->methyd1_index;
	which(methyd2$Symbol%in%genes)->methyd2_index;
	intersect(colnames(methyd1)[-(1:9)],colnames(methyd2)[-(1:9)])->shared_samples;
	setdiff(colnames(methyd1)[-(1:9)],colnames(methyd2)[-(1:9)])->methyd1_samples;
	setdiff(colnames(methyd2)[-(1:9)],colnames(methyd1)[-(1:9)])->methyd2_samples;
	if(length(methyd1_index)>0 && length(methyd2_index)>0){
		methyd1[methyd1_index,]->methyd1;
		methyd2[methyd2_index,]->methyd2;
		data.frame("Probe"=methyd2$Probe,methyd2[,10:ncol(methyd2)])->methyd2_matrix;
		merge(methyd1,methyd2_matrix,by.x="Probe",by.y="Probe")->methyd_merged;
		for(ss in shared_samples){
			paste(ss,"x",sep=".")->ss.x;
			paste(ss,"y",sep=".")->ss.y;
			unlist(lapply(1:nrow(methyd_merged),function(i){
				mean(methyd_merged[i,ss.x],methyd_merged[i,ss.y])->i_mean;
				i_mean;
			}))->ss_mean;
			data.frame(methyd_merged,"VX"=ss_mean)->methyd_merged;
			ss->names(methyd_merged)[ncol(methyd_merged)];
		}
		return(methyd_merged);
	}else{
		return(NULL);
	}
}
merge_two_methy_data(myd.methy_450k,myd.methy_27k,selected_geneSet)->myd.methy_merged;
cal_correlation_MET_EXP<-function(metd,expd,metd_start_column,expd_start_column){
	#--------map gene to probes
	c()->g_probes;
	names(table(metd$Symbol))->col1_table;
	for(s in col1_table){
		which(metd$Symbol%in%s)->tmp_probes;
		c(g_probes,paste(tmp_probes,collapse="_"))->g_probes;
	}
	col1_table->names(g_probes);
	#-------prepare matrix;
	t(metd[,metd_start_column:ncol(metd)])->metd_t;
	gsub("\\.","-",rownames(metd_t),perl=T)->rownames(metd_t);
	#-----------prepare same samples:
	intersect(expd$A0_Samples,rownames(metd_t))->shared_samples;
	if(length(shared_samples)==0){
		return(NULL);
	}
	expd[shared_samples,]->expd;
	metd_t[shared_samples,]->metd_t;
	metd$Probe->colnames(metd_t);
	#-------
	#-------calculate correlation
	colnames(expd)->expd_names;
	c()->metd_expd_cor;
	c()->metd_expd_cor_pvalue;
	c()->metd_probes;
	c()->metd_genes;
	c()->expd_genes;
	for(i in expd_start_column:ncol(expd)){
		expd_names[i]->g_name;
		g_probes[g_name]->g_name_probes;
		if(is.na(g_name_probes)){
			next;
		}
		as.numeric(unlist(strsplit(g_name_probes,"_")))->g_name_probes;
		for(p in g_name_probes){	
			#print(expd[,i]);
			#print(metd_t[,p]);flush.console();
			cor.test(expd[,i],metd_t[,p])->tmp_cor;
			c(metd_probes,as.character(colnames(metd_t)[p]))->metd_probes;
			c(metd_genes,g_name)->metd_genes;
			c(expd_genes,g_name)->expd_genes;
			c(metd_expd_cor,tmp_cor$estimate)->metd_expd_cor;
			c(metd_expd_cor_pvalue,tmp_cor$p.value)->metd_expd_cor_pvalue;
		}
	}
	data.frame("MetProbe"=metd_probes,"MetGene"=metd_genes,"ExpGene"=expd_genes,"Cor"=metd_expd_cor,"Cor.pvalue"=metd_expd_cor_pvalue)->res_df;
	return(res_df);
}
cal_correlation_MET_EXP(myd.methy_merged,tmp_expd_cluster,10,g_start_column)->MET_EXP.cor
prepare_corr_matrix<-function(cord,x,y,v,cutp){
	cord[,x]->x_names;
	cord[,y]->y_names;
	if(length(unique(x_names))!=nrow(cord)){
		paste(x_names,seq(1,length(x_names),1),sep=".")->x_names;
	}
	if(length(unique(y_names))!=nrow(cord)){
		paste(y_names,seq(1,length(y_names),1),sep=".")->y_names;
	}
	nrow(cord)->cord_nrow;
	for(i in 1:cord_nrow){
		
	}
}
draw_corr_pvalue_points<-function(cord,cutp,myd_colors,shows){
	cord[order(cord$Cor.pvalue),]->cord;
	unlist(lapply(cord$FDR,function(x){
		if(x<cutp){
			myd_colors[1];
		}else{
			myd_colors[2];
		}
	}))->col_colors;
	plot(x=cord$Cor,y=-log10(cord$FDR),col=col_colors,pch=20,xlab="Correlation",ylab="-log10(FDR)");
	abline(h=-log10(cutp),col=myd_colors[1],lty=3,lwd=0.8);
	text(x=cord$Cor[1:shows],y=-log10(cord$FDR[1:shows]),labels=cord[1:shows,2],cex=0.8)
	text(x=0,y=-log10(cutp),labels="FDR=0.05")
}
draw_corr_pvalue_points(MET_EXP.cor,0.05,c(myd_colors[1],"gray45"),22)
#for Methylation heatmap:
t(myd.methy_merged[,10:ncol(myd.methy_merged)])->myd.methy_merged_t;
myd.methy_merged$Probe->colnames(myd.methy_merged_t);
gsub("\\.","-",rownames(myd.methy_merged_t),perl=T)->rownames(myd.methy_merged_t);
data.frame("Sample"=rownames(myd.methy_merged_t),myd.methy_merged_t)->myd.methy_merged_t
merge(myd.methy_merged_t,tmp_exp.kcluster,by.x="Sample",by.y="Sample")->myd.methy_merged_t.kcluster;

#---sort by correlation p values:
MET_EXP.cor[order(MET_EXP.cor$Cor.pvalue),]->MET_EXP.cor;
p.adjust(MET_EXP.cor$Cor.pvalue)->MET_EXP.cor$FDR;
MET_EXP.cor[MET_EXP.cor$FDR<0.05,]->MET_EXP.cor_filter;
as.character(MET_EXP.cor_filter$MetProbe)->myd.methy_probes_sort;
generate_color<-function(x,myd_colors){
	names(table(x))->x_table;
	myd_colors[1:length(x_table)]->x_colors;
	x_table->names(x_colors);
	return(x_colors);
}
generate_gaps<-function(myd_clusters){
	table(myd_clusters)->cluster_count;
	c()->myd_gaps;
	0->tmp_gap;
	for(i in 1:(length(cluster_count)-1)){
		tmp_gap+cluster_count[i]->tmp_gap;
		c(myd_gaps,tmp_gap)->myd_gaps;
	}
	return(myd_gaps);
}
preorder_cnvcor_metcor_df_genes_heatmap<-function(myd,sample_clusters,clusters,geneOrder){
	myd[,as.character(geneOrder)]->myd_order;
	data.frame("Sample"=myd$Sample,myd_order)->myd_order;
	merge(myd_order,sample_clusters,by.x="Sample",by.y="Sample")->myd_merged;
	c()->order_index;
	for(i in clusters){
		which(myd_merged$KCluster==i)->tmp_index;
		hclust(dist(myd_merged[tmp_index,as.character(geneOrder)]))->tmp_index_hclust;
		tmp_index[tmp_index_hclust$order]->tmp_index;
		c(order_index,tmp_index)->order_index;
	}
	myd_merged[order_index,]->myd_merged;
	return(myd_merged)
}
#generate_gaps(tmp_exp.kcluster$KCluster)->gaps_col;
list("KCluster"=generate_color(tmp_exp.kcluster$KCluster,myd_colors))->col_colors;
c(colorRampPalette(rev(brewer.pal(11,"RdYlBu")[7:11]))(100),colorRampPalette(rev(brewer.pal(11,"RdYlBu")[1:6]))(250))->fill_colors;
preorder_cnvcor_metcor_df_genes_heatmap(myd.methy_merged_t,tmp_exp.kcluster,c("KC1","KC2"),myd.methy_probes_sort)->met_genes_changed_order;
pheatmap(t(met_genes_changed_order[,myd.methy_probes_sort]),color=fill_colors,cluster_rows=F,cluster_cols=F,gaps_col=249,show_colnames=F,labels_row=MET_EXP.cor_filter$MetGene,border_color=NA,annotation_col=met_genes_changed_order[,c("KCluster","KCluster")],annotation_colors=col_colors)->met_heatmap_gtable;

#----for m6A related genes:
library(corrplot)
read.table("m6a-genes.txt",header=T,sep="\t")->m6a_genes;
intersect(colnames(tmp_expd_cluster),m6a_genes$Symbols)->m6a_genes_shared;
tmp_expd_cluster[,c(as.character(unique(MET_EXP.cor_filter$MetGene)),m6a_genes_shared)]->m6a_genes_shared.matrix;
cor(m6a_genes_shared.matrix)->m6a_genes_shared.cor_matrix;
cor.mtest(m6a_genes_shared.matrix,conf.level=0.95)->corr_table.res1;
#----
change_cor_p2symbols<-function(cor_testd,cutps){
	cor_testd$p->cor_pd;
	c()->cor_p_symbols;
	for(i in 1:nrow(cor_pd)){
		for(j in 1:ncol(cor_pd)){
			if(cor_pd[i,j]<cutps[1]){
				c(cor_p_symbols,"***")->cor_p_symbols;
			}else if(cor_pd[i,j]<cutps[2] && cor_pd[i,j]>=cutps[1]){
				c(cor_p_symbols,"**")->cor_p_symbols;
			}else if(cor_pd[i,j]<cutps[3] && cor_pd[i,j]>=cutps[2]){
				c(cor_p_symbols,"*")->cor_p_symbols;
			}else{
				c(cor_p_symbols," ")->cor_p_symbols;
			}
		}
	}
	matrix(cor_p_symbols,ncol=ncol(cor_pd),byrow=T)->cor_p_symbols;
	#rownames(cor_pd)->rownames(cor_p_symbols);
	#colnames(cor_pd)->colnames(cor_p_symbols);
	return(cor_p_symbols);
}
change_cor_p2symbols(corr_table.res1,c(1e-5,1e-4,1e-3))->m6a_genes_shared.cor_p_symbols
colnames(m6a_genes_shared.matrix)->rownames(m6a_genes_shared.cor_p_symbols);
colnames(m6a_genes_shared.matrix)->colnames(m6a_genes_shared.cor_p_symbols);
#-----
prepare_text_colors<-function(catd,genes){
	unique(catd$m6a_process)->catd_types;
	myd_colors[1:length(catd_types)]->catd_types_col;
	catd_types->names(catd_types_col);
	#--
	as.character(catd$m6a_process)->m6a_process;
	catd$Symbol->names(m6a_process);
	m6a_process[genes]->m6a_process;
	catd_types_col[m6a_process]->catd_colors
	return(catd_colors)
}
prepare_text_colors(m6a_genes,m6a_genes_shared)->m6a_process_types_colors
#superheat(coxph_list.merged_df, bottom.label.text.angle = 90,bottom.label.text.size = 4,X.text=round(coxph_list.merged_df,2), heat.pal = c("white", "red"))
superheat(m6a_genes_shared.cor_matrix[as.character(unique(MET_EXP.cor_filter$MetGene)),m6a_genes_shared],X.text=m6a_genes_shared.cor_p_symbols[as.character(unique(MET_EXP.cor_filter$MetGene)),m6a_genes_shared],bottom.label.text.angle = 90,bottom.label.text.size = 4,heat.pal = c("blue","white", "red"),bottom.label.text.col=m6a_process_types_colors)
##########################################################################################################################################
#--------------------read CNV data
read.table("D:/TCGA-GBM/CNV/Copy Number Variation.nocnv.merge.txt",header=T,sep="\t",stringsAsFactors=F)->myd.cnv
myd.cnv[myd.cnv$Chromosome!="X",]->myd.cnv#***********
as.numeric(myd.cnv$Chromosome)->myd.cnv$Chromosome;#************
CNregions(seg=myd.cnv,epsilon=0,adaptive=FALSE,frac.overlap=0.5,rmSmallseg=TRUE,nProbes=5)->myd.cnv_merged
data.frame("Samples"=rownames(myd.cnv_merged),myd.cnv_merged)->myd.cnv_merged
myd.cnv_merged[grep("-10",myd.cnv_merged$Samples,invert=T),]->myd.cnv_merged_rm
fwrite(myd.cnv_merged_rm,"GBM.CNV_merged.txt",quote=F,sep="\t",row.names=F)
#----------
read.table("cnv_interval2genes.txt",header=T,sep="\t")->cnv_interval2genes
myd.cnv_merged[myd_samples,]->myd.cnv_merged;
#------------------------------------------------------------------for regions merged to one gene
fread("GBM.CNV_merged_regions2symbol.txt",header=T,sep="\t",stringsAsFactors=F)->myd.cnv_region2symbol;
as.data.frame(myd.cnv_region2symbol)->myd.cnv_region2symbol;
myd.cnv_region2symbol$Sample->rownames(myd.cnv_region2symbol);
myd.cnv_region2symbol[myd_samples,]->myd.cnv_region2symbol;

cal_correlation_CNV_EXP_region2symbol<-function(cnvd,expd,cnvd_start_column){
	#-------calculate correlation
	#cnvd$Sample->rownames(cnvd);
#	expd$A0_Samples->rownames(expd);
	cnvd[myd_samples,]->cnvd;
	expd[myd_samples,]->expd;
	colnames(cnvd)->cnvd_names;
	c()->cnvd_expd_cor;
	c()->cnvd_expd_cor_pvalue;
	c()->cnvd_genes;
	c()->cnvd_regions;
	c()->expd_genes;
	for(i in cnvd_start_column:ncol(cnvd)){
		cnvd_names[i]->g_name;
		as.character(g_name)->g_name;
		if(length(is.na(expd[1,g_name]))==0){
			next;
		}
		tryCatch({
			cor.test(cnvd[,i],expd[,g_name])->tmp_cor;
			c(cnvd_regions,cnvd_names[i])->cnvd_regions;
			c(cnvd_genes,as.character(g_name))->cnvd_genes;
			c(expd_genes,as.character(g_name))->expd_genes;
			c(cnvd_expd_cor,tmp_cor$estimate)->cnvd_expd_cor;
			c(cnvd_expd_cor_pvalue,tmp_cor$p.value)->cnvd_expd_cor_pvalue;
			},error=function(e){cat("Error:",g_name,"\n");}
		)
	}
	data.frame("CnvRegion"=cnvd_regions,"CnvGene"=cnvd_genes,"ExpGene"=expd_genes,"Cor"=cnvd_expd_cor,"Cor.pvalue"=cnvd_expd_cor_pvalue)->res_df;
	return(res_df);
}
t(myd.cnv_region2symbol[,-1])->myd.cnv_region2symbol_t;
apply(myd.cnv_region2symbol_t,1,function(x){2^x})->test;
data.frame("Samples"=myd.cnv_region2symbol$Sample,test)->myd.cnv_region2symbol_logTrans;
colnames(myd.cnv_region2symbol)->colnames(myd.cnv_region2symbol_logTrans);

cal_correlation_CNV_EXP_region2symbol(myd.cnv_region2symbol_logTrans,tmp_expd_cluster,2)->CNV_region2symbol_EXP.cor;
CNV_region2symbol_EXP.cor[!is.na(CNV_region2symbol_EXP.cor$Cor.pvalue),]->CNV_region2symbol_EXP.cor;
CNV_region2symbol_EXP.cor[order(CNV_region2symbol_EXP.cor$Cor.pvalue),]->CNV_region2symbol_EXP.cor;
p.adjust(CNV_region2symbol_EXP.cor$Cor.pvalue)->CNV_region2symbol_EXP.cor$FDR
#--------------draw  heatmap
CNV_region2symbol_EXP.cor[CNV_region2symbol_EXP.cor$FDR<0.05,]->CNV_region2symbol_EXP.cor_filter;
CNV_region2symbol_EXP.cor_filter$ExpGene[which(CNV_region2symbol_EXP.cor_filter$ExpGene%in%selected_geneSet)]->CNVCor_df_genes.shared_sort
as.character(CNVCor_df_genes.shared_sort)->CNVCor_df_genes.shared_sort;#35 genes

preorder_cnvcor_metcor_df_genes_heatmap(myd.cnv_region2symbol,tmp_exp.kcluster,c("KC1","KC2"),CNVCor_df_genes.shared_sort)->cnv_genes_changed_order;
generate_gaps(as.numeric(substr(cnv_genes_changed_order$KCluster,3,3)))->gaps_col;
list("KCluster"=generate_color(cnv_genes_changed_order$KCluster,myd_colors))->col_colors;
c(colorRampPalette(rev(brewer.pal(11,"RdYlBu")[7:11]))(270),colorRampPalette(rev(brewer.pal(11,"RdYlBu")[1:6]))(200))->fill_colors;
pheatmap(t(cnv_genes_changed_order[,CNVCor_df_genes.shared_sort]),col=fill_colors,cluster_rows=F,cluster_cols=F,gaps_col=gaps_col,show_colnames=F,show_rownames=T,border_color=NA,annotation_col=cnv_genes_changed_order[,c("KCluster","KCluster")],annotation_colors=col_colors)->cnv_heatmap_gtable;
#
draw_corr_pvalue_points(CNV_region2symbol_EXP.cor[which(CNV_region2symbol_EXP.cor$ExpGene%in%selected_geneSet),],0.05,c(myd_colors[1],"gray45"),35)
write.table(CNV_region2symbol_EXP.cor[which(CNV_region2symbol_EXP.cor$ExpGene%in%selected_geneSet),],"CNV_region2symbol_EXP.cor.txt",row.names=F,sep="\t",quote=F)
#--------------------------------------for SNV:
read.table("GBM.mutect2.gene_vcf2matrix.res",header=T,sep="\t",stringsAsFactors=F)->gene_snv_merged##
tmp_exp.kcluster->snv_nmf_3_cluster
gsub("-",".",snv_nmf_3_cluster$Sample)->snv_nmf_3_cluster$Sample;
snv_nmf_3_cluster$Sample->rownames(snv_nmf_3_cluster);
intersect(colnames(gene_snv_merged),snv_nmf_3_cluster$Sample)->exp_snv_shared_samples;
snv_nmf_3_cluster[exp_snv_shared_samples,]->snv_nmf_3_cluster;
cal_gene_snv_group_ftest<-function(snvd,groups){
	snvd->tmp_matrix;
	#--------------------------------
	c()->gName;
	c()->group1_Mut;
	c()->group1_Nor;
	c()->group2_Mut;
	c()->group2_Nor;
	c()->FisherP;
	#c()->padj;
	groups[groups$KCluster=="KC1",1]->group1;
	groups[groups$KCluster!="KC1",1]->group2;
	for(i in 1:nrow(tmp_matrix)){
		tmp_matrix[i,as.character(group1)]->group1_values;
		tmp_matrix[i,as.character(group2)]->group2_values;
		c(length(group1_values[group1_values!=0]),length(group1_values[group1_values==0]))->group1_count;
		c(length(group2_values[group2_values!=0]),length(group2_values[group2_values==0]))->group2_count;
		matrix(c(group1_count,group2_count),nrow=2,byrow=F)->group1_group2_table;
		#print(group1_group2_table);
		#flush.console();
		fisher.test(group1_group2_table)->group1_group2_table_test;
		#p.adjust(group1_group2_table_test$p.value,n=sum(group1_count,group2_count))->test_padj;
		c(FisherP,group1_group2_table_test$p.value)->FisherP;
		c(group1_Mut,group1_count[1])->group1_Mut;
		c(group1_Nor,group1_count[2])->group1_Nor;
		c(group2_Mut,group2_count[1])->group2_Mut;
		c(group2_Nor,group2_count[2])->group2_Nor;
		#c(padj,test_padj)->padj;
	}
	data.frame("gName"=tmp_matrix$gName,"FisherP"=FisherP,"Mut_C1"=group1_Mut,"Nor_C1"=group1_Nor,"Mut_C234"=group2_Mut,"Nor_C234"=group2_Nor)->res;
	p.adjust(res$FisherP,n=nrow(res))->res$padj;
	res[order(res$FisherP),]->res;
	return(res);
}
cal_gene_snv_group_ftest(gene_snv_merged,snv_nmf_3_cluster)->gene_snv_count;
write.table(gene_snv_count,"gene_snv_count.txt",quote=F,sep="\t",row.names=F)

which(gene_snv_merged$gName%in%selected_geneSet)->tmp_index;
gene_snv_merged[tmp_index,]->gene_snv_merged_filter;
#---------
#------for EGFR, TP53, KRAS genes mutation .
unlist(lapply(1:nrow(gene_snv_count),function(i){sum(gene_snv_count[i,c(3,5)])}))->gene_snv_count.sum;
gene_snv_count[order(gene_snv_count.sum,decreasing=T),]->gene_snv_count.by_sum;
c("EGFR","TP53","PTEN","NF1","MET","PIK3CA","RB1","IDH1")->target_genes;#https://www.stjude.org/disease/constitutional-mismatch-repair-deficiency.html
retrive_genes_mutation_samples<-function(snvd,genes){
	which(snvd$gName%in%genes)->genes_index;
	c()->mut_samples;
	for(i in genes_index){
		snvd[i,-1]->snv_values;
		which(snv_values!=0)->snv_index;
		colnames(snvd)[snv_index+1]->snv_index_samples;
		c(mut_samples,as.character(snv_index_samples))->mut_samples;
	}
	unique(mut_samples)->mut_samples;
	gsub("\\.","-",mut_samples)->mut_samples;
	data.frame("Sample"=colnames(snvd)[-1],"MMR_Status"=rep("WT",ncol(snvd)-1))->res;
	as.character(res$MMR_Status)->res$MMR_Status;
	gsub("\\.","-",res$Sample)->res$Sample;
	which(res$Sample%in%mut_samples)->samples_index;
	"Mut"->res[samples_index,"MMR_Status"];
	c("Sample",genes)->colnames(res);
	return(res);
}
prepare_genes_mutation_samples<-function(snvd,genes){
	for(i in 1:length(genes)){
		genes[i]->g;
		if(i==1){
			retrive_genes_mutation_samples(snvd,g)->g_samples;
		}else{
			retrive_genes_mutation_samples(snvd,g)->g2_samples;
			merge(g_samples,g2_samples,by.x="Sample",by.y="Sample")->g_samples
		}
	}
	paste(genes,"Mut",sep="_")->names(g_samples)[-1];
	return(g_samples);
}
prepare_genes_mutation_samples(gene_snv_merged,target_genes)->target_genes_samples;
#merge(tmp_expd_cluster,target_genes_samples,by.x="A0_Samples",by.y="Sample")->tmp_expd_cluster.MMR;

#----------
preoder_snv_gene_samples<-function(snvd,groups,clusters){
	data.frame("gName"=snvd$gName)->res;
	for(g in groups){
		unlist(lapply(g,function(x){which(clusters$KCluster==x)->i;i}))->tmp_index;
		gsub("-",".",as.character(clusters[tmp_index,1]))->g_samples
		snvd[,g_samples]->snvd_filter
		unlist(lapply(1:ncol(snvd_filter),function(x){which(snvd_filter[,x]!=0)->x_index;length(snvd_filter[x_index,x])}))->tmp_index_count;
		cbind(res,snvd_filter[,order(tmp_index_count)])->res;
	}
	res[,-1]->res;
	return(res);
}
preoder_snv_gene_samples(gene_snv_merged_filter,c("KC1","KC2"),snv_nmf_3_cluster)->gene_snv_merged_filter_sort;
unlist(lapply(c("KC1","KC2"),function(x){which(tmp_expd_cluster$KCluster==x)->i;i}))->tmp_index;
tmp_expd_cluster[,c(1:(g_start_column-1))]->tmp_annotation_col;
tmp_expd_cluster$KCluster->tmp_annotation_col$KCluster;
merge(tmp_annotation_col,target_genes_samples,by.x="A0_Samples",by.y="Sample")->tmp_annotation_col;#add mutation status
gsub("-","\\.",as.character(tmp_annotation_col$A0_Samples))->rownames(tmp_annotation_col)
#tmp_annotation_col[rownames(tmp_expd_cluster)[tmp_index],]->tmp_annotation_col
generate_gaps(as.numeric(substr(tmp_annotation_col$KCluster,3,3)))->gaps_col;

for(i in 1:nrow(gene_snv_merged_filter_sort)){
	which(gene_snv_merged_filter_sort[i,]!=0)->tmp_index;
	1->gene_snv_merged_filter_sort[i,tmp_index];
}    MET_Mut PIK3CA_Mut RB1_Mut IDH1_Mut
list("KCluster"=generate_color(cnv_genes_changed_order$KCluster,myd_colors),"EGFR_Mut"=generate_color(tmp_annotation_col$EGFR_Mut,brewer.pal(9,"Greys")[c(9,3)]),"TP53_Mut"=generate_color(tmp_annotation_col$TP53_Mut,brewer.pal(9,"Greys")[c(9,3)]),"PTEN_Mut"=generate_color(tmp_annotation_col$PTEN_Mut,brewer.pal(9,"Greys")[c(9,3)]),"NF1_Mut"=generate_color(tmp_annotation_col$NF1_Mut,brewer.pal(9,"Greys")[c(9,3)]),"MET_Mut"=generate_color(tmp_annotation_col$MET_Mut,brewer.pal(9,"Greys")[c(9,3)]),"PIK3CA_Mut"=generate_color(tmp_annotation_col$PIK3CA_Mut,brewer.pal(9,"Greys")[c(9,3)]),"RB1_Mut"=generate_color(tmp_annotation_col$RB1_Mut,brewer.pal(9,"Greys")[c(9,3)]),"IDH1_Mut"=generate_color(tmp_annotation_col$IDH1_Mut,brewer.pal(9,"Greys")[c(9,3)]),"A2_Event"=generate_color(tmp_annotation_col$A2_Event,brewer.pal(9,"Greys")[c(3,9)]),"A1_OS"=colorRampPalette(brewer.pal(9,"Purples"))(200))->col_colors;
pheatmap(gene_snv_merged_filter_sort,cluster_rows=F,cluster_cols=F,labels_row=gene_snv_merged$gName[which(gene_snv_merged$gName%in%selected_geneSet)],color=brewer.pal(9,"Greys")[c(2,6)],gaps_col=gaps_col,show_colnames=F,annotation_col=tmp_annotation_col[,c("KCluster","EGFR_Mut","TP53_Mut","PTEN_Mut","NF1_Mut","MET_Mut","PIK3CA_Mut","RB1_Mut","IDH1_Mut","A2_Event","A1_OS")],annotation_colors=col_colors,border_color="white")->gene_snv.p1
#as.grob(~barplot(gene_snv_count_filter$Mut_C3,horiz=T,))->gene_snv.p2;
#grid.arrange(gene_snv.p1$gtable,gene_snv.p2,ncol=2)
#data.frame("Symbol"=c(as.character(gene_snv_count_filter$gName),as.character(gene_snv_count_filter$gName)),"MutCount"=c(gene_snv_count_filter$Mut_C3,gene_snv_count_filter$Mut_C12),"KCluster"=rep(c("C3","C12"),each=nrow(gene_snv_count_filter)),"Index"=c(seq(nrow(gene_snv_count_filter),1,-1),seq(nrow(gene_snv_count_filter),1,-1)))->gene_snv_count.df;
#ggplot(gene_snv_count.df,aes(x=Index,y=MutCount))->p
#p+geom_bar(aes(fill=IRGCluster),stat="identity")+coord_flip()
which(gene_snv_count$gName%in%selected_geneSet)->gene_snv_count.selected_index;
as.matrix(gene_snv_count[,c("Mut_C1","Mut_C234")])->gene_snv_count.table
gene_snv_count$gName->rownames(gene_snv_count.table)
t(gene_snv_count.table)->gene_snv_count.table;
gene_snv_count.table[,gene_snv_count.selected_index]->gene_snv_count.table;
barplot(gene_snv_count.table[,seq(length(gene_snv_count.selected_index),1,-1)],col=c("red","blue"),beside=F,las=2,horiz=T,border=NA)
##########################################################################################################################################
##########################################################################################################################################
#------------HRGs risk model:TCGA_GSE13041_shared_coxph_genes
#-----------for TCGA set 
par(mfrow=c(1,2))
coxph_regression_prediction(tmp_expd_cluster,TCGA_GSE13041_shared_coxph_genes,c(1,3,5))->TCGA_shared_coxph_genes.coxph;
data.frame("gName"=TCGA_shared_coxph_genes.coxph$gName,"Coef"=TCGA_shared_coxph_genes.coxph$Coef)->TCGA_shared_coxph_genes.coxph_coef
calculate_risk_score(tmp_expd_cluster,TCGA_shared_coxph_genes.coxph_coef)->coxph.risk_score;
draw_survial_curve_custom(coxph.risk_score,5,500,brewer.pal(11,"Spectral")[c(4,10)])##p<0.001
#-----------for GSE13041 set 
coxph_regression_prediction(GSE13041_expd_processed,TCGA_GSE13041_shared_coxph_genes,c(1,3,5))->GSE13041_shared_coxph_genes.coxph;
calculate_risk_score(GSE13041_expd_processed,TCGA_GSE13041_shared_coxph_genes.coxph_coef)->GSE13041.risk_score;
draw_survial_curve_custom(GSE13041.risk_score,5,500,brewer.pal(11,"Spectral")[c(4,10)])##p=0.044
#----------for CGGA data:
preprocess_expd_CGGA<-function(expd,gStartColumn,aliveEvent){
	#-----------remove OS ==NA
	expd[!is.na(expd$A1_OS),]->expd;
	expd[!is.na(expd$A2_Event),]->expd;
	expd[expd$A1_OS!=0,]->expd;
	expd[expd$A1_OS>=30,]->expd;
	#----------remove IRGP value==1 or 0 in all samples-------
	c()->filter_colums;
	for(i in gStartColumn:ncol(expd)){
		length(expd[expd[,i]<(-10000),i])->failed_samples;
		if(failed_samples/nrow(expd)<0.5){
			c(filter_colums,i)->filter_colums
		}
	}
	expd[,c(1:(gStartColumn-1),filter_colums)]->expd.filter;
	print(length(filter_colums));
	flush.console();
	#---------status: 0->alive,1->death---------
	c()->status;
	for(i in 1:nrow(expd.filter)){
		if(expd.filter$A2_Event[i]==aliveEvent){
			c(status,0)->status;
		}else{
			c(status,1)->status; 
		}
	}
	status->expd.filter$Status
	return(expd.filter);
}
read.table("CGGA.mRNA_array_301.gene_level.20190701.factors.txt",header=T,sep="\t",stringsAsFactors=F)->CGGA_301array_expd;
preprocess_expd_CGGA(CGGA_301array_expd,15,0)->CGGA_301array_expd_processed;
coxph_regression_prediction(CGGA_301array_expd_processed,TCGA_GSE13041_shared_coxph_genes,c(1,3,5))->CGGA_301array_expd_processed.coxph;
calculate_risk_score(CGGA_301array_expd_processed,TCGA_GSE13041_shared_coxph_genes.coxph_coef)->CGGA_301array_expd_processed.risk_score;
draw_survial_curve_custom(CGGA_301array_expd_processed.risk_score,5,500,brewer.pal(11,"Spectral")[c(4,10)])##p<1e-5
#---------------draw RMS:
calc_rms<-function(myd_raw,f,tau,n.grid){
	myd_raw->myd;
	myd[myd[,f]!="",]->myd;
	myd[!is.na(myd[,f]),]->myd;
	#-------------
	myd$A1_OS/30->myd$A1_OS;
	#-------------
	myd[,f]->marker;
	seq(0,1,length=n.grid)->marker.pp;
	quantile(marker,marker.pp)->marker.qq;
	data.frame(marker=marker)->fitdat.df;
	data.frame(marker=marker.qq)->newdat.df;
	f->colnames(newdat.df)
	#-------------
	if(length(myd$A1_OS)>1){
		Surv(as.numeric(myd$A1_OS),as.numeric(myd$Status))->myd.surv;
	}else if(length(myd$OS)>1){
		Surv(as.numeric(myd$OS),as.numeric(myd$Status))->myd.surv;
	}
	as.formula(paste("myd.surv~",f))->tmp_formula;
	coxph(tmp_formula,data=myd)->cox.model;
	#---calculate log rank p AND HR:
	summary(cox.model)->tmp.coxph.summary;
	tmp.coxph.summary$logtest[3]->univar_coxph_logtest;
	tmp.coxph.summary$coefficients[,2]->univar_coxph_HR;
	tmp.coxph.summary$conf.int[,3]->univar_coxph_low95;
	tmp.coxph.summary$conf.int[,4]->univar_coxph_high95;
	paste(round(univar_coxph_low95,3),round(univar_coxph_high95,3),sep="-")->ci_95;
	paste(round(univar_coxph_HR,3),ci_95,sep="(")->univar_coxph_HR
	#---calculate rms: 
	summary(survfit(cox.model,newdata=newdat.df),rmean=tau)->rms.calc;
	rms.calc$table[,"*rmean"]->rms.mean;
	return(list("pp"=marker.pp,"rmean"=rms.mean,"p.value"=univar_coxph_logtest,"HR"=univar_coxph_HR));
}
draw_rms_curve_factors<-function(myd,factors,tau,n.grid,myd_colors){
	for(i in 1:length(factors)){	
		if(i==1){
			factors[i]->f;
			calc_rms(myd,f,tau,n.grid)->f_res;
			plot(f_res$pp,f_res$rmean,type='p',pch=20,cex=0.8,ylab="RMS,mo",xlab=paste("Percentile of Scores"),col=myd_colors[i]);
		}else{
			factors[i]->f;
			calc_rms(myd,f,tau,n.grid)->f_res;
			lines(f_res$pp,f_res$rmean,type='p',pch=20,cex=1,col=myd_colors[i])
		}
		paste("HR:",f_res$HR,")",sep="")->HR_labels;
		paste("p:",round(f_res$p.value,3))->p_labels;
		paste(HR_labels,p_labels,sep=",")->text_labels;
		text(x=0.3,y=15*i,labels=text_labels,font=3,col=myd_colors[i]);
	}
	abline(h=f_res$rmean[c(21,81)],col=myd_colors[i],lty=3,lwd=0.8);
	abline(v=c(0.2,0.8),col=myd_colors[i],lty=2,lwd=0.8)
	legend("topright",legend=factors,col=myd_colors[1:length(factors)],lty=1,cex=1.2);
}
merge(tmp_expd_cluster,coxph.risk_score[,c(1,4)],by.x="A0_Samples",by.y="A0_Samples")->tmp_expd_cluster.combined_riskscore;
#
merge(GSE13041_expd_processed,GSE13041.risk_score[,c(1,4)],by.x="A0_Samples",by.y="A0_Samples")->GSE13041.combined_riskscore;
#
merge(CGGA_301array_expd_processed,CGGA_301array_expd_processed.risk_score[,c(1,4)],by.x="A0_Samples",by.y="A0_Samples")->CGGA_301array_expd_processed.combined_riskscore;
par(mfrow=c(1,3))
draw_rms_curve_factors(tmp_expd_cluster.combined_riskscore,c("RiskScore","Age"),120,100,myd_colors)
draw_rms_curve_factors(GSE13041.combined_riskscore,c("RiskScore","AGE"),120,100,myd_colors)
draw_rms_curve_factors(CGGA_301array_expd_processed.combined_riskscore,c("RiskScore","Age"),120,100,myd_colors)


############################################################
#----------sample information summary：
Sample_info_summary<-function(myd,columns){
	list()->summary.res;
	for(i in columns){
		as.data.frame(table(myd[i]))->tmp.df;
		c(colnames(myd)[i],"Values")->colnames(tmp.df);
		print(tmp.df);
	}
}
Sample_info_summary(tmp_expd_cluster,c(3,8))
Sample_info_summary(GSE13041_expd_processed,c(3,7))
Sample_info_summary(CGGA_301array_expd_processed,c(3,8))

list(GSE13041_expd_processed,tmp_expd_cluster,CGGA_301array_expd_processed)->datasets;
lapply(datasets,function(x){
	summary(x[which(x$Status==0),"A1_OS"])->x1;
	summary(x[which(x$Status!=0),"A1_OS"])->x2;
	list(x1,x2);
});
#----------
draw_os_by_status<-function(myd){
	myd[myd$Status==0,]->myd.0;
	myd[myd$Status==1,]->myd.1;
	brewer.pal(9,"Greys")->myd_colors;
	seq(0,max(myd$A1_OS)*1.2,365)->myd.breaks
	max(table(cut(myd$A1_OS,breaks=myd.breaks)))*1.1->myd.max;
	seq(0,max(myd.0$A1_OS)*1.2,365)->myd.0.breaks;
	max(table(cut(myd.0$A1_OS,breaks=myd.0.breaks)))*1.1->myd.0.max;
	seq(0,max(myd.1$A1_OS)*1.2,365)->myd.1.breaks;
	max(table(cut(myd.1$A1_OS,breaks=myd.1.breaks)))*1.1->myd.1.max;
	par(mfrow=c(1,3),col="gray");
	hist(myd$A1_OS,breaks=myd.breaks,axes=F,xlab="",main="Overall",col=myd_colors[4])#,border=NA
	axis(side=1,at=myd.breaks[1:(length(myd.breaks)-1)]+365/2,labels=paste(myd.breaks[-1]/365," year",sep=""),las=2,pos=0)
	axis(side=2,at=seq(0,myd.max,20),labels=seq(0,myd.max,20),pos=0)
	hist(myd.0$A1_OS,breaks=myd.0.breaks,axes=F,xlab="",main="Alive Samples",col=myd_colors[4]);
	axis(side=1,at=myd.0.breaks[1:(length(myd.0.breaks)-1)]+365/2,labels=paste(myd.0.breaks[-1]/365," year",sep=""),las=2,pos=0)
	axis(side=2,at=seq(0,myd.0.max,10),labels=seq(0,myd.0.max,10),pos=0)
	hist(myd.1$A1_OS,breaks=myd.1.breaks,axes=F,xlab="",main="Dead Samples",col=myd_colors[4]);
	axis(side=1,at=myd.1.breaks[1:(length(myd.1.breaks)-1)]+365/2,labels=paste(myd.1.breaks[-1]/365," year",sep=""),las=2,pos=0)
	axis(side=2,at=seq(0,myd.1.max,10),labels=seq(0,myd.1.max,10),pos=0)
}
draw_os_by_status(tmp_expd_cluster)
draw_os_by_status(GSE13041_expd_processed)
draw_os_by_status(CGGA_301array_expd_processed)











