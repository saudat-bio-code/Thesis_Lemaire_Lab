###
###

.pkgname <- "BSgenome.Crobusta.HT.KY"

.seqnames <- sub("\\.fa","",list.files("build",pattern="\\.fa"))

.circ_seqs <- "GCF_000224145.1"

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
        organism="Ciona robusta",
        common_name="Vase tunicate",
        genome="GCF_000224155.1",
        provider="NCBI",
        release_date="2019",
        source_url="https://www.ncbi.nlm.nih.gov/assembly/GCF_000224145.1",
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

    old_objname <- "Crobusta"
    assign(old_objname, bsgenome, envir=ns)
    namespaceExport(ns, old_objname)
}

