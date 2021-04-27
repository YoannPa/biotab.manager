
#' Subset clinical biotabs following a type of data and its manifest.
#'
#' @param Biotab.RDS A \code{list} of tibbles containing clinical data.
#' @param Manifest   A \code{data.frame} containing TCGA information about
#'                   samples and data to use for subsetting clinical data
#' @param data.type  A \code{character} to specify the type of data to use.
#' @return A \code{list} with biotabs as \code{tibbles} and a summary table as a
#'         \code{data.frame} of the subset.
#' @author Yoann Pageaud.
#' @export

subset.clinical.biotab <- function(Biotab.RDS, Manifest, data.type){
  #Get Clinical Tables size
  clinic.tab.size<-unlist(lapply(Biotab.RDS,nrow))

  #Select Clinical data for which there are HM450K methylation data
  if(data.type == "MAF"){
    MAF.patient.ID<-unique(substr(Manifest$Tumor_Sample_Barcode,
                                  start = 1, stop = 12))
    #Select Clinical data for which there are somatic mutation data
    Biotab.selection<-lapply(names(Biotab.RDS),
                              function(tab) {
                                Biotab.RDS[[tab]][
                                  Biotab.RDS[[tab]]$bcr_patient_barcode
                                  %in% MAF.patient.ID,]
                              })
  } else {
    Biotab.selection<-lapply(names(Biotab.RDS),
                             function(tab) {
                               Biotab.RDS[[tab]][
                                 Biotab.RDS[[tab]]$bcr_patient_barcode
                                 %in% Manifest$`Patient ID`,]
                             })
  }
  names(Biotab.selection)<-names(Biotab.RDS)
  #Get Clinical Tables size after selecting HM450K patients
  clinic.tab.selection.size<-unlist(lapply(Biotab.selection,nrow))
  #Track biotab selection process
  data.selection.df<-data.frame("Input.Biotab" = clinic.tab.size,
                                clinic.tab.selection.size)
  colnames(data.selection.df)[2]<-data.type
  return(list("Selected.biotab" = Biotab.selection,
              "Selection.summary" = data.selection.df))
}

#' Multi-step clinical biotabs subsetting following a type of data and its
#' manifest.
#'
#' @param Biotab.RDS     A \code{list} of tibbles containing clinical data.
#' @param list.manifests A \code{dataframe} vector containing multiple manifests
#'                       in a given order to successively subset clinical data
#'                       following the different manifests.
#' @return A \code{list} with biotabs as \code{tibbles} and a summary table as a
#'         \code{data.frame} of the different steps of subset.
#' @author Yoann Pageaud.
#' @export

multi.subset.biotab <- function(Biotab.RDS, list.manifests){
  init.res<-unlist(lapply(Biotab.RDS,nrow))

  subset.list<-lapply(seq_along(list.manifests), function(i){
    subset.res<-subset.clinical.biotab(Biotab.RDS = Biotab.RDS,
                                       Manifest = list.manifests[[i]],
                                       data.type = names(list.manifests)[i])
    #Change object value outside the scope of lapply() use <<-
    Biotab.RDS<<-subset.res$Selected.biotab
    return(subset.res$Selection.summary[[2]])
  })

  subset.df<-as.data.frame(do.call(cbind,subset.list))
  subset.df<-cbind(init.res,subset.df)
  colnames(subset.df)<-c("Input.Biotab",names(list.manifests))

  return(list("Selected.biotab" = Biotab.RDS, "Selection.summary" = subset.df))
}
