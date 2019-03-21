##IMPORTS
Imports = c("TCGAbiolinks","simsalapar")
lapply(Imports, library, character.only = T)

##FUNCTIONS

##get_TCGAproj_biotab ##########################################################

get_TCGAproj_biotab<-function(project.id,clin.data.dir){
  #Keep only TCGA projects
  processable.proj<-project.id[grepl("^TCGA-*",project.id)]
  
  #Display projects that are not from TCGA
  unprocessable.proj<-project.id[!(project.id %in% processable.proj)]
  cat(paste(c("WARNING! Cannot process the following projects:",
              unprocessable.proj,"\n")))
  #Change working dir
  main_wd<-getwd()
  setwd(clin.data.dir)
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
