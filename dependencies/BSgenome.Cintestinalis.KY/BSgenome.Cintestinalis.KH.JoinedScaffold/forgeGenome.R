library(BSgenome)
#kh <- readLines("JoinedScaffold")
#file.start <- grep('>',kh)
#file.end <- c(file.start[-1]-1,length(kh))
#file.names <- paste0("KH/",sub('>','',kh[file.start]),'.fa')
#mapply(
#  function(start,end,name) writeLines(kh[start:end],name),
#  file.start,
#  file.end,
#  file.names
#)
#forgeBSgenomeDataPkg("DESCRIPTION")
forgeBSgenomeDataPkg("de1")
#forgeMaskedBSgenomeDataPkg("de1")
