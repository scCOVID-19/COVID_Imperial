library(Seurat)
library(plyr)
library(dplyr)
library(patchwork)
library(harmony)
library(future)
library(ggplot2)
library(trqwe)
library(unix)
library(readr)
library(reshape2)
library(rmcorr)
library(pheatmap)
library(dplyr)
library(Seurat)
library(data.table)
library(numDeriv)
library(tidyr)
source("/rds/general/user/emacdona/projects/covid19-transcriptome/live/sc_rnaseq/scripts/functions/glmm_functions_4.R")

setwd("/rds/general/user/emacdona/projects/covid19-transcriptome/live/sc_rnaseq/DA/")
data<-mcreadRDS("../all_compartments_qc_new_annotation_291024_DA.rds",mc.core=4)

###################
#Run GLMM on cell type composition of all compartments combined over severity
#####################
meta<-data.frame(data@meta.data)
meta$cell_type<-meta$level_2_new
meta<-meta[meta$case_control=="POSITIVE",]
meta$WHO_temp_severity <- factor(meta$WHO_temp_severity, levels = c('mild', 'moderate', 'severe', 'critical'))
meta$WHO_temp_severity_group <- factor(meta$WHO_temp_severity, levels = c('mild', 'moderate', 'severe', 'critical'), labels = c('mild_moderate', 'mild_moderate', 'severe_critical', 'severe_critical'))


ranef_list_all_time <- list()
annot<-"cell_type"
meta$annotation <- meta[,annot]
Y = table(meta$sample_id,meta$annotation)

myMetadata <- meta[!duplicated(meta$sample_id),]
myMetadata$calc_age_log2 <- log2(myMetadata$calc_age)
myMetadata$sample_id_orig.ident <- paste(myMetadata$sample_id,myMetadata$orig.ident)

myMetadata$random_effect_term <- myMetadata$individual_id

metadata <- myMetadata[,c("sample_id","calc_age_log2","sex","corrected_ethnicity","centre","WHO_temp_severity_group","random_effect_term")] 
Y <- Y[rownames(Y)%in%metadata$sample_id,]
nsamples = nrow(Y)
ncells = ncol(Y)

metadataExp=cbind(metadata[rep(match(rownames(Y),as.character(metadata$sample_id)),ncells),],Celltype=rep(colnames(Y),rep(nsamples,ncells)))
res.prop=glmer(I(c(Y))~
                   (1|Celltype)
                 +calc_age_log2
                 +(1|sex)
                 +(1|corrected_ethnicity)
                 +(1|centre)
                 +(1|WHO_temp_severity_group)
                 +(1|random_effect_term)
                 
                 +(calc_age_log2-1|Celltype)
                 +(1|sex:Celltype)
                 +(1|corrected_ethnicity:Celltype)
                 +(1|centre:Celltype)
                 +(1|WHO_temp_severity_group:Celltype)
                 +(1|random_effect_term:Celltype)
                 ,
                 family=poisson,data=metadataExp,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

# standard errors of standard deviations (squre root of the variance parameters)
devfun = update(res.prop, devFunOnly=T)
pars = getME(res.prop, c("theta","fixef"))
hess = hessian(devfun, unlist(pars))
sdse.prop = data.frame(sd=unlist(pars), se=sqrt(diag(solve(hess))))

# posterior means and their standard deviations
res.prop.ranef = ranef(res.prop)
  
getCondVal=function(res.prop.ranef, id, ncells, nfactors=2, celltype){
        tmp = data.frame(res.prop.ranef)[data.frame(res.prop.ranef)[[1]]==id,]
        if(length(grep(":",tmp$grp))==0){
                cnam = matrix(as.character(tmp$term),ncells)[1,]
                rnam = matrix(as.character(tmp$grp),ncells)[,1]
        }else if(nfactors==2){
a=unique(matrix(unlist(strsplit(as.character(tmp$grp),":")),2)[1,])
b=unique(matrix(unlist(strsplit(as.character(tmp$grp),":")),2)[2,])
tmp=tmp[match(paste(a[rep(1:length(a),rep(length(b),length(a)))],b[rep(1:length(b),length(a))],sep=":"),tmp$grp),]
tmp$grp=paste(a[rep(1:length(a),rep(length(b),length(a)))],b[rep(1:length(b),length(a))],sep=":")
                cnam = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),2)[1,],ncells)[1,]
                rnam = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),2)[2,],ncells)[,1]
        }else if(nfactors==3){
a=unique(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[1,])
b=unique(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[2,])
C=unique(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[3,])
ab=paste(a[rep(1:length(a),rep(length(b),length(a)))],b[rep(1:length(b),length(a))],sep=":")
abc=paste(ab[rep(1:length(ab),rep(length(C),length(ab)))],C[rep(1:length(C),length(ab))],sep=":")
tmp=tmp[match(abc,tmp$grp),]
tmp$grp=abc
                cnam1 = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[1,],ncells)[1,]
                cnam2 = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[2,],ncells)[1,]
                cnam = paste(cnam1,cnam2,sep=":")
                rnam = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[3,],ncells)[,1]
        }
        condval = matrix(tmp$condval,ncells)
        condsd  = matrix(tmp$condsd, ncells)
    condval[is.na(condval)]=0
    condsd[is.na(condsd)]=1
        rownames(condval)=rownames(condsd)=rnam
    condval=condval[match(celltype,rnam),]
    condsd =condsd[match(celltype,rnam),]
        colnames(condval)=colnames(condsd)=cnam
        lfsr = pnorm(condval,0,condsd)
        lfsr[is.na(lfsr)]=1
        lfsr[lfsr>0.5]=1-lfsr[lfsr>0.5]
        list(condval=condval,  lfsr=lfsr)
}
# obtaining posterior mean and
# local false sign rate (lfsr)

getCondVal_modRik <- function(res.prop.ranef, id, ncells, nfactors=2, celltype, reference="D-1", referenceRow="none", fdr=F){
  tmp = data.frame(res.prop.ranef)[data.frame(res.prop.ranef)[[1]]==id,]
  if(length(grep(":",tmp$grp))==0){
    cnam = matrix(as.character(tmp$term),ncells)[1,]
    rnam = matrix(as.character(tmp$grp),ncells)[,1]
  }else if(nfactors==2){
    a=unique(matrix(unlist(strsplit(as.character(tmp$grp),":")),2)[1,])
    b=unique(matrix(unlist(strsplit(as.character(tmp$grp),":")),2)[2,])
    tmp=tmp[match(paste(a[rep(1:length(a),rep(length(b),length(a)))],b[rep(1:length(b),length(a))],sep=":"),tmp$grp),]
    tmp$grp=paste(a[rep(1:length(a),rep(length(b),length(a)))],b[rep(1:length(b),length(a))],sep=":")
    cnam = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),2)[1,],ncells)[1,]
    rnam = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),2)[2,],ncells)[,1]
  }else if(nfactors==3){
    a=unique(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[1,])
    b=unique(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[2,])
    C=unique(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[3,])
    ab=paste(a[rep(1:length(a),rep(length(b),length(a)))],b[rep(1:length(b),length(a))],sep=":")
    abc=paste(ab[rep(1:length(ab),rep(length(C),length(ab)))],C[rep(1:length(C),length(ab))],sep=":")
    tmp=tmp[match(abc,tmp$grp),]
    tmp$grp=abc
    cnam1 = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[1,],ncells)[1,]
    cnam2 = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[2,],ncells)[1,]
    cnam = paste(cnam1,cnam2,sep=":")
    rnam = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[3,],ncells)[,1]
  }
  condval = matrix(tmp$condval,ncells)
  condsd = matrix(tmp$condsd, ncells)
  condval[is.na(condval)]=0
  condsd[is.na(condsd)]=1
  rownames(condval)=rownames(condsd)=rnam
  condval=condval[match(celltype,rnam),]
  condsd =condsd[match(celltype,rnam),]
  colnames(condval)=colnames(condsd)=cnam
  if (referenceRow!="none") {
    condval = t(apply(condval,1,function(x) x-condval[referenceRow,]))
  }
  condval <- condval-condval[,reference] # Check this with Ni
  condsd <- sqrt(condsd^2+condsd[,reference]^2)
  lfsr = pnorm(condval,0,condsd)
  lfsr[is.na(lfsr)]=1
  lfsr[lfsr>0.5]=1-lfsr[lfsr>0.5]
  lfsr <- lfsr*2
  if (fdr) { lfsr[,colnames(lfsr)!=reference] <- matrix(p.adjust(lfsr[,colnames(lfsr)!=reference],method="fdr"),ncol = ncol(lfsr[,colnames(lfsr)!=reference])) } # Better to exclude reference from fdr correction..
  list(condval=condval, lfsr=lfsr, sd=condsd)
}

postmean = cbind(
    getCondVal(res.prop.ranef,"Celltype",ncells,celltype=colnames(Y),nfactors = 1)[[1]][,1,drop=F],
    NA,
    getCondVal(res.prop.ranef,"sex:Celltype",ncells,celltype=colnames(Y),nfactors = 2)[[1]],
    NA,
    getCondVal(res.prop.ranef,"corrected_ethnicity:Celltype",ncells,celltype=colnames(Y),nfactors = 2)[[1]],
    NA,
    getCondVal(res.prop.ranef,"centre:Celltype",ncells,celltype=colnames(Y),nfactors = 2)[[1]],
    NA,
    getCondVal_modRik(res.prop.ranef,"WHO_temp_severity_group:Celltype",ncells,celltype=colnames(Y),nfactors = 2,reference="mild_moderate")[[1]]
  )
  
  ltsr = cbind(
    getCondVal(res.prop.ranef,"Celltype",ncells,celltype=colnames(Y),nfactors = 1)[[2]][,1,drop=F],
    NA,
    getCondVal(res.prop.ranef,"sex:Celltype",ncells,celltype=colnames(Y),nfactors = 2)[[2]],
    NA,
    getCondVal(res.prop.ranef,"corrected_ethnicity:Celltype",ncells,celltype=colnames(Y),nfactors = 2)[[2]],
    NA,
    getCondVal(res.prop.ranef,"centre:Celltype",ncells,celltype=colnames(Y),nfactors = 2)[[2]],
    NA,
    getCondVal_modRik(res.prop.ranef,"WHO_temp_severity_group:Celltype",ncells,celltype=colnames(Y),nfactors = 2,reference="mild_moderate")[[2]]
  )

ranef_list_all_time[[paste0(annot,"_ranef")]] <- res.prop.ranef
ranef_list_all_time[[paste0(annot,"_postmean")]] <- postmean
ranef_list_all_time[[paste0(annot,"_ltsr")]] <- ltsr

postmean <- ranef_list_all_time[[paste0("cell_type_postmean")]]
ltsr <- ranef_list_all_time[[paste0("cell_type_ltsr")]]

maxL2FC <- 4
myLtsr <- as.data.frame(ltsr[,14,drop=FALSE])
myPostmean <- as.data.frame(postmean[,14,drop=FALSE])
myPostmean$annotation <- rownames(myPostmean)
myLtsr$annotation <- rownames(myLtsr)

#######################################
# want to output table near here

longPost <- as.data.frame(pivot_longer(myPostmean,names_to = "sampleInfo", values_to = "post",-annotation))
longLtsr <- as.data.frame(pivot_longer(myLtsr,names_to = "sampleInfo", values_to = "fdr",-annotation))

rownames(longPost) <- paste(longPost$annotation,longPost$sampleInfo)
longPost[paste(longLtsr$annotation,longLtsr$sampleInfo),"fdr"] <- longLtsr$fdr
longPost$post_log2 <- log2(exp(longPost$post))
write.table(longPost,"Rik_model_severity_results_level_2.tsv",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

##############################

longPost <- as.data.frame(pivot_longer(myPostmean,names_to = "sampleInfo", values_to = "post",-annotation))
longLtsr <- as.data.frame(pivot_longer(myLtsr,names_to = "sampleInfo", values_to = "fdr",-annotation))

rownames(longPost) <- paste(longPost$annotation,longPost$sampleInfo)
longPost[paste(longLtsr$annotation,longLtsr$sampleInfo),"fdr"] <- longLtsr$fdr

longPost$post <- log2(exp(longPost$post))

longPost$annotation <- factor(longPost$annotation,levels=rev(myPostmean$annotation))

longPost$post[longPost$post>maxL2FC] <- maxL2FC
longPost$post[longPost$post< -maxL2FC] <- -maxL2FC

colorLimit <- max(abs(longPost$post))
figureRatio <- length(unique(longPost$annotation))/length(unique(longPost$sample))

df<-longPost
df$label<-df$annotation
fdrMax<-10^-3
df$fdr[df$fdr<fdrMax] <- fdrMax

ggplot(df,aes(x=sampleInfo,y=label,col=post,size=-log10(fdr))) + 
    scale_color_gradientn(colors=c("blue","white","red"),breaks=c(-maxL2FC,-maxL2FC/2,0,maxL2FC/2,maxL2FC),limits=c(-maxL2FC,maxL2FC)) + 
    scale_radius(range = c(.1,6),limits = c(0,3),breaks=c(1,2,3)) + 
    geom_point(shape=16) + 
    xlab("WHO Severity Group")+
    labs(colour="Log2 Fold Change over\nMild-Moderate")+
    theme_classic()+
    theme(aspect.ratio = figureRatio,
         axis.text.x = element_text(size=8),
         axis.text.y=element_text(size=8),
         axis.title.y=element_blank(),
          axis.title.x.y=element_text(size=10),
         legend.title=element_text(size=10),
         legend.text=element_text(size=8))
ggsave("Riks_code_severity_adjust_size.pdf",height=8,width=6)
