##IMPORTS
Imports = c("TCGAbiolinks","simsalapar")
lapply(Imports, library, character.only = T)

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
      if(rm.data){ unlink("GDCdata/", recursive = TRUE) }
      setwd(main_wd)
    } else { if(rm.data){ unlink("GDCdata/", recursive = TRUE) } }
    #Make column description table
    description.tbl<-data.frame(
      "column.names" = colnames(Prep.query[[1]]),
      "column.descriptions" = c(as.matrix(Prep.query[[1]][1,])))
    description.tbl<-description.tbl[sort(description.tbl$column.names),]
    #Make clinical data element table
    CDE_ID.tbl<-data.frame(
      "column.names" = colnames(Prep.query[[1]]),
      "clinical.data.element.ID" = c(as.matrix(Prep.query[[1]][2,]))
    )
    CDE_ID.tbl<-CDE_ID.tbl[sort(CDE_ID.tbl$column.names),]
    #Store clin.biotab
    biotab<-as.data.frame(Prep.query[[1]][-c(1,2),])
  } else { stop("Unknown TCGA project.")}
  return(list("biotab" = biotab, "description.table" = description.tbl,
              "CDE_ID.table" = CDE_ID.tbl, "warnings.query" = warn.query,
              "warnings.DL" = warn.DL, "warnings.preps" = warn.prep))
}

#' @description Returns a list of biotabs, one biotab per TCGA project.
#' 
#' @param proj.ids A \code{character} vector containing TCGA projects names.
#' @param dl.dir   A \code{character} specifying the path where downloaded
#'                 clinical raw data should be saved.
#' @value A \code{type} object returned description.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' @references

#TODO: use get.TCGA.biotab in the function.
#TODO: replace for loop by lapply().
get.ls.TCGA.biotab<-function(proj.ids,dl.dir){ 
  #Keep only TCGA projects
  processable.proj<-proj.ids[grepl("^TCGA-*",proj.ids)]
  #Display projects that are not from TCGA
  if(length(proj.ids[!(proj.ids %in% processable.proj)]) != 0){
    unprocessable.proj<-proj.ids[!(proj.ids %in% processable.proj)]
    cat(paste(c("WARNING! Cannot process the following projects:",
                unprocessable.proj,"\n")))
  } else { unprocessable.proj<-NULL }
  #Change working dir
  main_wd<-getwd()
  setwd(dl.dir)
  #Set lists
  list.biotab<-list()
  list.warn.query<-list()
  list.warn.DL<-list()
  list.warn.prep<-list()
  #Download TCGA projects clinical data and create biotabs
  for (proj in processable.proj) {
    cat(paste0(proj,"\n"))
    #Create GDC query
    list.warn.query[[proj]] <- tryCatch.W.E(query <- GDCquery(project = proj,
                                                    data.category = "Clinical",
                                                    data.type = "Clinical data",
                                                    legacy = TRUE,
                                                    file.type = "txt"))$warning
    #Download Clin. data
    list.warn.DL[[proj]] <-tryCatch.W.E(GDCdownload(query))$warning
    #Get clinical biotab 
    list.warn.prep[[proj]] <- tryCatch.W.E(Prep.query <-
                                             GDCprepare(query))$warning
    #Store clin.biotab
    list.biotab[[proj]]<-Prep.query
  }
  #Move back to original Working dir 
  setwd(main_wd)
  #Create result with biotabs and warnings
  res<-list("Biotabs" = list.biotab,
            "unprocessed.projects" = unprocessable.proj,
            "warnings.query" = list.warn.query, "warnings.DL" = list.warn.DL,
            "warnings.preps" = list.warn.prep)
  return(res)
}

##save.list.RDS ################################################################

save.list.RDS<-function(list,dir.save){
  if(dir.exists(dir.save) == F){
    dir.create(dir.save)
  }
  save.ret<-lapply(seq_along(1:length(list)),
                   function(i) {saveRDS(list[[i]],
                                        file = file.path(dir.save,
                                                         paste0(names(list[i]),
                                                                ".RDS")))})
}

##save.unproc.proj #############################################################

save.unproc.proj<-function(biotab.res,dir.save){
  if(file.exists(file.path(dir.save,"Unprocessed_Projects.txt"))){
    file.remove(file.path(dir.save,"Unprocessed_Projects.txt"))
  }
  write.res<-lapply(biotab.res[["unprocessed.projects"]], write,
                    file.path(dir.save,"Unprocessed_Projects.txt"),append = T)
}

##save.warnings ################################################################

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
