###
###

.pkgname <- "BSgenome.Cintestinalis.KY"

.seqnames <- NULL

.circ_seqs <- "GCA_009617815.1"

.mseqnames <- NULL

.onLoad <- function(libname, pkgname)
{
    if (pkgname != .pkgname)
        stop("package name (", pkgname, ") is not ",
             "the expected name (", .pkgname, ")")
    extdata_dirpath <- system.file("extdata", package=pkgname,
                                   lib.loc=libname, mustWork=TRUE)

    ## Make and export BSgenome object.
    bsgenome <- BSgenome(
        organism="Ciona intestinalis",
        common_name="Vase tunicate",
        genome="GCA_009617815.1",
        provider="NCBI",
        release_date="2019",
        source_url="https://www.ncbi.nlm.nih.gov/assembly/GCA_009617815.1",
        seqnames=.seqnames,
        circ_seqs=.circ_seqs,
        mseqnames=.mseqnames,
        seqs_pkgname=pkgname,
        seqs_dirpath=extdata_dirpath
    )

    ns <- asNamespace(pkgname)

    objname <- pkgname
    assign(objname, bsgenome, envir=ns)
    namespaceExport(ns, objname)

    old_objname <- "Cintestinalis"
    assign(old_objname, bsgenome, envir=ns)
    namespaceExport(ns, old_objname)
}

