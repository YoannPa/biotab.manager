#!/usr/bin/env Rscript
library("optparse")

opt_parser = OptionParser(description = "
##_Prepare Clinical Biotabs_####################################################
Version = '1.0'
Date = '2019-02-25'
Author = 'Yoann PAGEAUD'
Maintainer = 'Yoann PAGEAUD (yoann.pageaud@gmail.com)'
Dependencies = c('R version 3.5.3 (2019-03-11)',
'RStudio Version 1.1.463 – © 2009-2018','TCGAbiolinks','data.table','parallel')
Description = 'Subset Clinical Data based on (epi-)genomics data,create clinical
tables saved as CSV files.'
################################################################################
");
if (is.null(parse_args(opt_parser)$file) != TRUE){print_help(opt_parser);quit()}

##IMPORTS
Imports = c("TCGAbiolinks","data.table","parallel")
lapply(Imports, library, character.only = T)
source("subset_biotabs.R")

##PARAMETERS
TCGA.Somatic.MAF<-fread("TCGA/somatic_mutation_calls/mc3.v0.2.8.PUBLIC.maf.gz")
Methylation.data.dir<-"TCGA/downloaded/methylation"
Biotabs.dir<-"TCGA/Biotabs"
Clin.data.HM450K.dir<-"data/Clinical_Data_Subset_HM450K.RDS"
Merged.Clin.HM450K.dir<-"data/Merged_Clinical_Data_HM450K.RDS"
Clin.tbl.dir<-"data/Clinical_tables_HM450K"

#ANALYSIS
#List TGCGA project names
TCGA.projects<-list.files(Methylation.data.dir)
Projects.dir<-lapply(TCGA.projects, function(i){
  file.path(Methylation.data.dir,i)
})

#List data type for each project
Projects.data.type<-lapply(Projects.dir,list.files)

#Subset Clinical Biotabs following Data of Interest
Subset.TCGA.Biotabs<-lapply(seq_along(Projects.dir),function(i){
  cat(paste0("TCGA-",TCGA.projects[i],"\n"))
  #Load HM450K manifest
  if("hm450" %in% Projects.data.type[[i]]){
    if(file.exists(file.path(Projects.dir[[i]],"hm450","manifest.RDS"))){
      HM450K.Manifest<-readRDS(file.path(Projects.dir[[i]],"hm450",
                                         "manifest.RDS"))
    } else {
      warning(paste("No HM450k manifest for",TCGA.projects[i]))
    }
  }
  #Load WGBS manifest
  if("WGBS" %in% Projects.data.type[[i]]){
    if(file.exists(file.path(Projects.dir[[i]],"WGBS","manifest.RDS"))){
      WGBS.Manifest<-readRDS(file.path(Projects.dir[[i]],"WGBS",
                                       "manifest.RDS"))
    }
  }
  #Load Clinical Biotab
  if(file.exists(file.path(Biotabs.dir,paste0("TCGA-",TCGA.projects[i],
                                              ".RDS")))){
    Biotab<-readRDS(file.path(Biotabs.dir,paste0("TCGA-",TCGA.projects[i],
                                                 ".RDS")))
    Biotab<-lapply(Biotab, setDT)
  } else {
    warning(paste("No Biotab for",TCGA.projects[i]))
  }
  #Keep clinical.biotab overlapping HM450K data then MAF data
  if(exists("HM450K.Manifest")){
    # list.manifests<-list("HM450k"=HM450K.Manifest,"MAF"=TCGA.Somatic.MAF)
    # subset.HM450K.MAF<-multi.subset.biotab(Biotab.RDS = Biotab,
    #                                        list.manifests = list.manifests)
    subset.HM450K<-subset.clinical.biotab(Biotab.RDS = Biotab,
                                          Manifest = HM450K.Manifest,
                                          data.type = "HM450K")
    subset.HM450K.MAF<-subset.clinical.biotab(Biotab.RDS =
                                         subset.HM450K$Selected.biotab,
                                          Manifest = TCGA.Somatic.MAF,
                                          data.type = "MAF")
    #If WGBS Data Keep clinical.biotab overlapping WGBS
    if(exists("WGBS.Manifest")){
      subset.WGBS<-subset.clinical.biotab(Biotab.RDS =
                                            subset.HM450K.MAF$Selected.biotab,
                                          Manifest = WGBS.Manifest,
                                          data.type = "WGBS")
      Total.Selection.summary<-cbind(
        subset.HM450K$Selection.summary,
        "MAF" = subset.HM450K.MAF$Selection.summary$MAF,
        "WGBS"= subset.WGBS$Selection.summary$WGBS)
      
      res<-list("subset.HM450K" = subset.HM450K,
                "subset.HM450K.MAF" = subset.HM450K.MAF,
                "subset.WGBS" = subset.WGBS,
                "Total.Selection.summary" = Total.Selection.summary)
    } else {
      Total.Selection.summary<-subset.HM450K.MAF$Selection.summary
      res<-list("subset.HM450K" = subset.HM450K,
                "subset.HM450K.MAF" = subset.HM450K.MAF,
                "Total.Selection.summary" = Total.Selection.summary)
    }
  } else { #If No HM450K manifest
    res<-NULL
  }
  return(res)
})

#Give name to elements of the list
names(Subset.TCGA.Biotabs)<-paste0("TCGA-",TCGA.projects)
#Remove NULL elements
Subset.TCGA.Biotabs[sapply(Subset.TCGA.Biotabs,is.null)]<-NULL
#Save
saveRDS(Subset.TCGA.Biotabs, file = Clin.data.HM450K.dir)

##Focus on the HM450K Methylation data only first
#Get colnames of all tables
colnames.biotab<-lapply(seq_along(Subset.TCGA.Biotabs), function(i){
  lapply(Subset.TCGA.Biotabs[[i]]$subset.HM450K$Selected.biotab,colnames)
})
#Get column names (keys) in common between all tables and all tumor types
biotab.keys<-Reduce(intersect,lapply(seq_along(colnames.biotab),
                        function(i){Reduce(intersect, colnames.biotab[[i]])}))

#Merge all biotabs for each tumor type
TCGA.merged.Biotabs<-lapply(seq_along(Subset.TCGA.Biotabs), function(i){
  cat(paste0(names(Subset.TCGA.Biotabs)[i],"\n"))
  #Get HM450K methylation biotab
  biotab.hm450K<-Subset.TCGA.Biotabs[[i]]$subset.HM450K$Selected.biotab
  #Set keys in tables
  biotab.with.keys<-lapply(biotab.hm450K,setkeyv,cols = biotab.keys)
  #Rename colnames following the belonging table before merging
  biotab.renamed<-lapply(seq_along(biotab.with.keys), function(j){
    #Get colnames not being keys
    cols.not.keys<-match(biotab.keys,colnames(biotab.with.keys[[j]]))
    #Rename
    colnames(biotab.with.keys[[j]])[-c(cols.not.keys)]<-
      paste0(colnames(biotab.with.keys[[j]]),".",
             names(biotab.with.keys)[j])[-c(cols.not.keys)]
    biotab.with.keys[[j]]
  })
  #Rename Biotab
  names(biotab.renamed)<-names(biotab.with.keys)
  #Merge tables following table keys
  # (alow.cartesian to extend rows for the multiple drugs given to each patient)
  biotab.merged<-Reduce(function(x, y) merge(x,y,all = T,allow.cartesian=TRUE),
                        biotab.renamed)
  biotab.merged #Return merged biotab
})

##Clean empty columns
#Replace [Unknown],[Not Available],[Not Applicable],[Not Evaluated] by NA.
TCGA.merged.Biotabs<-lapply(seq_along(TCGA.merged.Biotabs), function(i){
  TCGA.merged.Biotabs[[i]][TCGA.merged.Biotabs[[i]]== "[Unknown]" |
                             TCGA.merged.Biotabs[[i]] ==
                             "[Not Available]" |
                             TCGA.merged.Biotabs[[i]] ==
                             "[Not Applicable]" |
                             TCGA.merged.Biotabs[[i]] ==
                             "[Not Evaluated]"]<-NA
  TCGA.merged.Biotabs[[i]]
})
names(TCGA.merged.Biotabs)<-names(Subset.TCGA.Biotabs)
#Remove columns containing only NAs
TCGA.merged.Biotabs<-lapply(TCGA.merged.Biotabs, function(i){
  i[,which(colSums(is.na(i)) == nrow(i)) := NULL]
})
#Save
saveRDS(TCGA.merged.Biotabs, file = Merged.Clin.HM450K.dir)

#Save separate clinical table for each tumor type
dir.create(Clin.tbl.dir)
lapply(seq_along(TCGA.merged.Biotabs), function(i){
  write.csv(TCGA.merged.Biotabs[[i]],
            file=file.path(Clin.tbl.dir,paste0(names(TCGA.merged.Biotabs)[i],
                                               ".csv")))
})
