#!/usr/bin/env Rscript
library("optparse")

opt_parser = OptionParser(description = "
##_Get TCGA Projects clinical biotab_###########################################
Version = '1.0'
Date = '2019-02-27'
Author = 'Yoann PAGEAUD'
Maintainer = 'Yoann PAGEAUD (yoann.pageaud@gmail.com)'
Dependencies = c('R version 3.5.3 (2019-03-11)',
'RStudio Version 1.1.463 – © 2009-2018')
Description = 'Retrieve all clinical data from TCGA projects and save several
logs generated during retrieval.'
################################################################################
");
if (is.null(parse_args(opt_parser)$file) != TRUE){print_help(opt_parser);quit()}

##IMPORTS
source("make_TCGAproj_biotab.R")

##PARAMETERS
Clinical.dir<-"TCGA"
Projects_list<-TCGAbiolinks:::getGDCprojects()$project_id
Biotab.dir<-file.path(Clinical.dir,"Biotabs")

##ANALYSIS
#Get projects biotabs
res<-get_TCGAproj_biotab(project.id = Projects_list,
                          clin.data.dir = Clinical.dir)

#Save Biotabs in separated rds files
save.list.RDS(res$Biotabs,dir.save = Biotab.dir)

#Save Names of unprocessed projects
save.unproc.proj(biotab.res = res, dir.save = Biotab.dir)

#Save Warnings
save.warnings(biotab.res = res,dir.save = Biotab.dir)
