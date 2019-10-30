##IMPORTS
Imports = c("TCGAbiolinks","simsalapar")
lapply(Imports, library, character.only = T)
source("src/fun_stringeval.R")

##FUNCTIONS

#' @description Function description.
#' 
#' @param proj.id A \code{character} specifying a TCGA project name.
#' @param dl.dir  A \code{character} specifying the path where downloaded
#'                clinical raw data should be saved. If NULL, a folder 'GDCdata'
#'                is created in the working directory (Default: dl.dir = NULL).
#' @param rm.data A \code{logical} specifying whether to remove the downloaded
#'                data or not once the biotab retrieved.
#' @value A \code{list} containing 6 elements:
#'        - patient.clinicals: a \code{data.frame} containing the project
#'          clinical data.
#'        - description.table: a \code{data.frame} containing short description
#'          of the biotab column names.
#'        - CDE_ID.table: a \code{data.frame} containing the Clinical Data
#'          Element IDs for each column of the biotab.
#'        - biotab: a \code{list} of the downloaded raw biotab tables.
#'        - warnings.query: a \code{list} of the possible warnings related to
#'          GDCquery().
#'        - warnings.DL: a \code{list} of the possible warnings related to
#'          GDCdownload().
#'        - warnings.preps: a \code{list} of the possible warnings related to
#'          GDCprepare().
#' @author Yoann Pageaud.
#' @export
#' @examples
#' @references

get.TCGA.biotab<-function(proj.id, dl.dir=NULL, rm.data=FALSE){
  if(grepl(pattern = "^TCGA-*", x = proj.id)){
    cat(paste0(proj.id,"\n"))
    #Create GDC query
    warn.query <- tryCatch.W.E(
      query <- GDCquery(project = proj.id, data.category = "Clinical",
                        data.type = "Clinical data", legacy = TRUE,
                        file.type = "txt"))$warning
    #Download clinical data
    if(!is.null(dl.dir)){ main_wd<-getwd() ; setwd(dl.dir) }
    warn.DL <-tryCatch.W.E(GDCdownload(query))$warning
    #Get clinical biotab 
    warn.prep <- tryCatch.W.E(Prep.query <- GDCprepare(query))$warning
    if(!is.null(dl.dir)){
      if(rm.data){
        unlink("GDCdata/", recursive = TRUE)
        if(file.exists("MANIFEST.txt")){ file.remove("MANIFEST.txt") }
      }
      setwd(main_wd)
    } else {
      if(rm.data){
        unlink("GDCdata/", recursive = TRUE)
        if(file.exists("MANIFEST.txt")){ file.remove("MANIFEST.txt") }
      }
    }
    #Find "clinical patient" table
    clname<-
      names(Prep.query)[grepl(pattern = "_patient_",x = names(Prep.query))]
    #Make column description table
    description.tbl<-data.frame(
      "column.names" = colnames(Prep.query[[clname]]),
      "column.descriptions" =c(as.matrix(Prep.query[[clname]][1,])))
    description.tbl<-description.tbl[sort(description.tbl$column.names),]
    #Make clinical data element table
    CDE_ID.tbl<-data.frame(
      "column.names" = colnames(Prep.query[[clname]]),
      "clinical.data.element.ID" = c(as.matrix(Prep.query[[clname]][2,]))
    )
    CDE_ID.tbl<-CDE_ID.tbl[sort(CDE_ID.tbl$column.names),]
    #Store clin.biotab
    patient_clinicals<-as.data.frame(Prep.query[[clname]][-c(1,2),])
    #If any upper case patient UUID convert to lower case 
    if(any(is.upper(str = patient_clinicals$bcr_patient_uuid))){
      patient_clinicals$bcr_patient_uuid[
        is.upper(str = patient_clinicals$bcr_patient_uuid)]<-
        tolower(patient_clinicals$bcr_patient_uuid[
          is.upper(str=patient_clinicals$bcr_patient_uuid)])
    }
  } else { stop("Unknown TCGA project.")}
  return(list(
    "patient.clinicals"=patient_clinicals, "description.table"=description.tbl,
    "CDE_ID.table" = CDE_ID.tbl, "biotab" = Prep.query,
    "warnings.query" = warn.query, "warnings.DL" = warn.DL,
    "warnings.preps" = warn.prep))
}

#' @description Returns a list of biotabs, one biotab per TCGA project.
#' 
#' @param proj.ids A \code{character} vector containing TCGA projects names.
#' @param dl.dir   A \code{character} specifying the path where downloaded
#'                 clinical raw data should be saved. If NULL, a folder
#'                 'GDCdata' is created in the working directory
#'                 (Default: dl.dir = NULL).
#' @param rm.data  A \code{logical} specifying whether to remove the downloaded
#'                 data or not once the biotab retrieved.               
#' @value A \code{list} containing results for each TCGA projects. A result is a
#'        \code{list} itself which contains 6 elements:
#'        - biotab: a \code{data.frame} containing the project clinical data.
#'        - description.table: a \code{data.frame} containing short description
#'          of the biotab column names.
#'        - CDE_ID.table: a \code{data.frame} containing the Clinical Data
#'          Element IDs for each column of the biotab.
#'        - warnings.query: a \code{list} of the possible warnings related to
#'          GDCquery().
#'        - warnings.DL: a \code{list} of the possible warnings related to
#'          GDCdownload().
#'        - warnings.preps: a \code{list} of the possible warnings related to
#'          GDCprepare().
#' @author Yoann Pageaud.
#' @export
#' @examples
#' @references

get.ls.TCGA.biotab<-function(proj.ids, dl.dir=NULL, rm.data=FALSE){ 
  #Keep only TCGA projects
  processable.proj<-proj.ids[grepl("^TCGA-*",proj.ids)]
  #Display projects that are not from TCGA
  if(length(proj.ids[!(proj.ids %in% processable.proj)]) != 0){
    unprocessable.proj<-proj.ids[!(proj.ids %in% processable.proj)]
    cat(paste(c("WARNING: Cannot process the following projects: ",
                unprocessable.proj,"\n")))
  } else { unprocessable.proj<-NULL }
  res.ls<-lapply(
    X=processable.proj, FUN=get.TCGA.biotab, dl.dir=dl.dir, rm.data=rm.data)
  names(res.ls)<-processable.proj
  if(!is.null(unprocessable.proj)){
    warning(paste(c(
      "The following projects have been discarded from the output: ",
      unprocessable.proj, "\n"))) 
  }
  return(res.ls)
}

##OBSOLETE #####################################################################

#' @description Function description.
#' 
#' @param param1 A \code{type} parameter description.
#' @value A \code{type} object returned description.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' @references

save.unproc.proj<-function(biotab.res,dir.save){
  if(file.exists(file.path(dir.save,"Unprocessed_Projects.txt"))){
    file.remove(file.path(dir.save,"Unprocessed_Projects.txt"))
  }
  write.res<-lapply(biotab.res[["unprocessed.projects"]], write,
                    file.path(dir.save,"Unprocessed_Projects.txt"),append = T)
}

#' @description Function description.
#' 
#' @param param1 A \code{type} parameter description.
#' @value A \code{type} object returned description.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' @references

save.warnings<-function(biotab.res,dir.save){
  if(file.exists(file.path(dir.save,"Warnings.txt"))){
    file.remove(file.path(dir.save,"Warnings.txt"))
  }
  warning.names<-names(biotab.res)[grepl("warnings",names(biotab.res))]
  write.res<-lapply(warning.names, function(i) {
    write(paste0("##",i),file = file.path(dir.save,"Warnings.txt"),append = T)
    if(length(biotab.res[[i]]) > 0){
      write.res<-lapply(seq_along(1:length(biotab.res[[i]])), function(j) {
        write(names(biotab.res[[i]])[j], file = file.path(dir.save,
                                                          "Warnings.txt"),
              append = T)
        write(paste0("\t",biotab.res[[i]][[names(biotab.res[[i]])[j]]]$message),
              file = file.path(dir.save,"Warnings.txt"), append = T)
      })
      write("\n",file = file.path(dir.save,"Warnings.txt"),append = T)
    } else {
      write("No Warnings",file = file.path(dir.save,"Warnings.txt"),append = T)
      write("\n",file = file.path(dir.save,"Warnings.txt"),append = T)
    }
  })
}
