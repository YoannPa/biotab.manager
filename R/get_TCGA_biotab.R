##IMPORTS
Imports = c("TCGAbiolinks", "simsalapar", "data.table")
lapply(Imports, library, character.only = T)
source("src/fun/fun_stringeval.R")

##FUNCTIONS

#' @description Returns a biotab for a given TCGA project ID.
#' 
#' @param proj.id A \code{character} specifying a TCGA project name.
#' @param dl.dir  A \code{character} specifying the path where downloaded
#'                clinical raw data should be saved. If NULL, a folder 'GDCdata'
#'                is created in the working directory (Default: dl.dir = NULL).
#' @param rm.data A \code{logical} specifying whether to remove the downloaded
#'                data or not once the biotab retrieved.
#' @return A \code{list} containing 6 elements:
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

get.TCGA.biotab <- function(proj.id, dl.dir = NULL, rm.data = FALSE){
  if(grepl(pattern = "^TCGA-*", x = proj.id)){
    cat(paste0(proj.id,"\n"))
    #Create GDC query
    warn.query <- tryCatch.W.E(
      query <- GDCquery(project = proj.id, data.category = "Clinical",
                        data.type = "Clinical data", legacy = TRUE,
                        file.type = "txt"))$warning
    #Download clinical data
    if(!is.null(dl.dir)){ main_wd <- getwd() ; setwd(dl.dir) }
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
    clname <- names(Prep.query)[
      grepl(pattern = "_patient_",x = names(Prep.query))]
    #Make column description table
    description.tbl <- data.frame(
      "column.names" = colnames(Prep.query[[clname]]),
      "column.descriptions" = c(as.matrix(Prep.query[[clname]][1,])))
    description.tbl <- description.tbl[sort(description.tbl$column.names), ]
    #Make clinical data element table
    CDE_ID.tbl <- data.frame(
      "column.names" = colnames(Prep.query[[clname]]),
      "clinical.data.element.ID" = c(as.matrix(Prep.query[[clname]][2, ]))
    )
    CDE_ID.tbl <- CDE_ID.tbl[sort(CDE_ID.tbl$column.names), ]
    #Store clin.biotab
    patient_clinicals <- as.data.frame(Prep.query[[clname]][-c(1,2), ])
    #If any upper case patient UUID convert to lower case 
    if(any(is.upper(str = patient_clinicals$bcr_patient_uuid))){
      patient_clinicals$bcr_patient_uuid[
        is.upper(str = patient_clinicals$bcr_patient_uuid)] <-
        tolower(patient_clinicals$bcr_patient_uuid[
          is.upper(str = patient_clinicals$bcr_patient_uuid)])
    }
  } else { stop("Unknown TCGA project.") }
  return(list(
    "patient.clinicals" = patient_clinicals,
    "description.table" = description.tbl, "CDE_ID.table" = CDE_ID.tbl,
    "biotab" = Prep.query, "warnings.query" = warn.query,
    "warnings.DL" = warn.DL, "warnings.preps" = warn.prep))
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
#' @return A \code{list} containing results for each TCGA projects. A result is
#'         a \code{list} itself which contains 6 elements:
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

get.ls.TCGA.biotab <- function(proj.ids, dl.dir = NULL, rm.data = FALSE){ 
  #Keep only TCGA projects
  processable.proj <- proj.ids[grepl("^TCGA-*", proj.ids)]
  #Display projects that are not from TCGA
  if(length(proj.ids[!(proj.ids %in% processable.proj)]) != 0){
    unprocessable.proj <- proj.ids[!(proj.ids %in% processable.proj)]
    cat(paste(c("WARNING: Cannot process the following projects: ",
                unprocessable.proj, "\n")))
  } else { unprocessable.proj <- NULL }
  res.ls <- lapply(X = processable.proj, FUN = get.TCGA.biotab, dl.dir = dl.dir,
                   rm.data = rm.data)
  names(res.ls) <- processable.proj
  if(!is.null(unprocessable.proj)){
    warning(paste(c(
      "The following projects have been discarded from the output: ",
      unprocessable.proj, "\n"))) 
  }
  return(res.ls)
}

#' Convert TCGA project ID into ICGC cohort and vice versa.
#' 
#' @param project A \code{character} vector of valid TCGA projects and/or
#'                ICGC cohorts.
#' @param to      A \code{character} specifying the output format of project
#'                names
#'                (Default: to = "TCGA"; Supported: to = c("TCGA", "ICGC")).
#' @return A \code{character} vector of the formated project names.
#' @author Yoann Pageaud.
#' @export

ICGC.to.TCGA.project <- function(project = NULL, to = "TCGA"){
  #Get TCGA project codes
  TCGA.codes <- TCGAbiolinks::getGDCprojects()$project_id
  TCGA.codes <- sort(TCGA.codes[grepl(pattern = "^TCGA-.+", x = TCGA.codes)])
  #Make table matching ICGC and TCGA codes.
  icgc.to.tcga.code <- data.frame("ICGC_code" = paste0(gsub(
    pattern = "^TCGA-", replacement = "", x = TCGA.codes), "-US"),
    "TCGA_code" = TCGA.codes)
  
  if(to == "TCGA"){
    #Convert detected ICGC codes into TCGA codes
    if(any(project %in% icgc.to.tcga.code$ICGC_code)){
      #Reorder table
      icgc.to.tcga.code <- icgc.to.tcga.code[order(match(
        icgc.to.tcga.code$ICGC_code, project)), ]
      #Subset table
      project[project %in% icgc.to.tcga.code$ICGC_code] <- icgc.to.tcga.code[
        icgc.to.tcga.code$ICGC_code %in% project, ]$TCGA_code
    }
    #Discard unsupported project codes
    if(any(grepl(pattern = "^TCGA-.+$", x = project) == FALSE)){
      warning("some invalid project codes have been discarded.")
      project <- project[grepl(pattern = "^TCGA-.+$", x = project)]
      if(length(project) == 0){
        stop("no valid project codes detected in 'project'.")
      }
    }
  } else if(to == "ICGC"){
    #Convert detected TCGA codes into ICGC codes
    if(any(project %in% icgc.to.tcga.code$TCGA_code)){
      #Reorder table
      icgc.to.tcga.code <- icgc.to.tcga.code[order(match(
        icgc.to.tcga.code$TCGA_code, project)), ]
      #Subset table
      project[project %in% icgc.to.tcga.code$TCGA_code] <- icgc.to.tcga.code[
        icgc.to.tcga.code$TCGA_code %in% project, ]$ICGC_code
    }
    #Discard unsupported project codes
    if(any(project %in% icgc.to.tcga.code$ICGC_code == FALSE)){
      warning("some invalid project codes have been discarded.")
      project <- project[project %in% icgc.to.tcga.code$ICGC_code]
      if(length(project) == 0){
        stop("no valid project codes detected in 'project'.")
      }
    }
  } else {
    stop("unsupported value for 'to'. Supported values are 'ICGC' or 'TCGA'.")
  }
  return(project)
}

#' Collects patients clinical data from specific TCGA projects. 
#' 
#' @param biotabs.dir     A \code{character} to specify the directory where all
#'                        biotabs have been downloaded using
#'                        get.ls.TCGA.biotab() or get.TCGA.biotab().
#' @param proj.id         A \code{character} vector of supported TCGA projects
#'                        and/or ICGC cohorts.
#' @param patient.barcode A \code{character} vector of valid TCGA patient IDs.
#' @param as.dt           A \code{logical} to specify whether the output table
#'                        should be a tibble (Default: as.dt = FALSE) or a
#'                        data.table (as.dt = TRUE).
#' @return A \code{list} with selected projects containing TCGA patients
#'         clinical data as a \code{tibble} or a \code{data.table}.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' @references

get.TCGA.clinical <- function(
  biotabs.dir, proj.id, patient.barcode, as.dt = FALSE){
  #Check biotabs.dir content
  proj.clinical.available <- gsub(
    pattern = ".RDS$", replacement = "", x = list.files(biotabs.dir))
  if(length(proj.clinical.available) == 0){
    stop("No biotab detected in the directory specified for 'biotabs.dir'")
  }
  #Convert all project IDs to TCGA
  proj.id <- ICGC.to.TCGA.project(project = proj.id)
  #Check if all project codes match available biotabs
  if(any(proj.id %in% proj.clinical.available == FALSE)){
    proj.not.avail <- proj.id[!proj.id %in% proj.clinical.available]
    stop(paste(
      "the following TCGA project(s) are not available in the directory of biotabs:",
      paste0(paste0(proj.not.avail, collapse = ", "),"."),
      "Please download the missing biotab(s) using get.ls.TCGA.biotab().",
      sep = "\n"))
  }
  #Checking patient barcodes
  if(any(grepl(pattern = "^TCGA-.{7}$", x = patient.barcode) == FALSE)){
    warning("some invalid patient barcodes will be discarded.")
    patient.barcode <- patient.barcode[
      grepl(pattern = "^TCGA-.{7}$", x = patient.barcode)]
    if(length(patient.barcode) == 0){
      stop("no valid patient barcodes detected in 'patient.barcode'.")
    }
  }
  ls.proj <- lapply(X = proj.id, FUN = function(p){
    load.path <- file.path(biotabs.dir, paste0(p,".RDS"))
    biotab <- readRDS(load.path)
    lapply(X = biotab, function(b){
      if(as.dt){
        data.table::as.data.table(b[
          b$bcr_patient_barcode %in% patient.barcode, ])  
      } else {
        b[b$bcr_patient_barcode %in% patient.barcode, ]
      }
    })
  })
  names(ls.proj) <- proj.id
  return(ls.proj)
}

#' Subset clinical data from collected patients clinical data, and return it
#' into a data.table.
#' 
#' @param TCGA.clinical A \code{list} of TCGA projects clinical data generated
#'                      using get.TCGA.clinical().
#' @param biotab.type   A \code{character} to specify the type of biotab to
#'                      extract data from
#'                      (Default: biotab.type = "clinical_patient"). Supported
#'                      biotab.type can match any table name without the project
#'                      suffix.
#' @param columns       A \code{character} vector specifying the column names
#'                      you want to retrieve clinical data from. Columns name
#'                      can appear in some projects only. If columns = 'all',
#'                      all existing columns throughout projects listed will be
#'                      retrieved.
#' @param to.ICGC       A \code{logical} to specify whether the TCGA project
#'                      should be converted into ICGC cohorts (to.ICGC = TRUE)
#'                      or not (Default: to.ICGC = FALSE).
#' @return A \code{type} object returned description.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' @references

TCGA.clinical.as.dt <- function(
  TCGA.clinical, biotab.type = "clinical_patient", columns = "all",
  to.ICGC = FALSE){
  #Retrieve all columns if "all"
  if(length(columns) == 1){
    if(columns == "all"){
      colnames.used <- lapply(X = TCGA.clinical, FUN = function(p){
        biotab.selected <- names(p)[
          grep(pattern = paste0("^", biotab.type, ".*$"), x = names(p))]
        colnames(p[[biotab.selected]])
      })
      columns <- unique(unlist(colnames.used))
    }  
  }
  #Subset data from data.type of the projects list
  ls.TCGA.clin <- lapply(X = TCGA.clinical, FUN = function(p){
    biotab.selected <- names(p)[
      grep(pattern = paste0("^", biotab.type, ".*$"), x = names(p))]
    col.avail <- columns[columns %in% colnames(p[[biotab.selected]])]
    col.not.avail <- columns[!columns %in% colnames(p[[biotab.selected]])]
    clin.res <- p[[biotab.selected]][, ..col.avail, ]
    dt.miss <- setNames(data.table(matrix(
      nrow = nrow(clin.res), ncol = length(col.not.avail))), col.not.avail)
    clin.res <- cbind(clin.res, dt.miss)
  })
  #Convert to TCGA project to ICGC cohort
  if(to.ICGC){
    names(ls.TCGA.clin) <- ICGC.to.TCGA.project(
      project = names(ls.TCGA.clin), to = "ICGC")
  }
  #Format result as a data.table
  dt.TCGA.clin <- rbindlist(
    l = ls.TCGA.clin, use.names = TRUE, idcol = "Cohorts")
  return(dt.TCGA.clin)
}