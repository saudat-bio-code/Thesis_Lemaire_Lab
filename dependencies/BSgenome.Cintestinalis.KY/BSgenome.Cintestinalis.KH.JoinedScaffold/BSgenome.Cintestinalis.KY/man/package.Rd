\name{BSgenome.Cintestinalis.KY}
\docType{package}

\alias{BSgenome.Cintestinalis.KY-package}
\alias{BSgenome.Cintestinalis.KY}
\alias{Cintestinalis}

\title{Ciona intestinalis HY}

\description{
  A BSgenome package for the Ciona intestinalis genome
}

\details{
  
}

\note{
  This BSgenome data package was made from the following source data files:
  \preformatted{
-- information not available --
  }

  See \code{?\link[BSgenome]{BSgenomeForge}} and the BSgenomeForge
  vignette (\code{vignette("BSgenomeForge")}) in the \pkg{BSgenome}
  software package for how to make a BSgenome data package.
}

\author{Yutaka Satou, Katsuhiko Mineta, Michio Ogasawara, Yasunori Sasakura, Eiichi Shoguchi, Keisuke Ueno, Lixy Yamada, Jun Matsumoto, Jessica Wasserscheid, Ken Dewar, Graham B Wiley, Simone L Macmil, Bruce A Roe, Robert W Zeller, Kenneth EM Hastings, Patrick Lemaire, Erika Lindquist, Toshinori Endo2, Kohji Hotta and Kazuo Inaba}

\seealso{
  \itemize{
    \item \link[BSgenome]{BSgenome} objects and the
          \code{\link[BSgenome]{available.genomes}} function
          in the \pkg{BSgenome} software package.
    \item \link[Biostrings]{DNAString} objects in the \pkg{Biostrings}
          package.
    \item The BSgenomeForge vignette (\code{vignette("BSgenomeForge")})
          in the \pkg{BSgenome} software package for how to make a BSgenome
          data package.
  }
}

\examples{
BSgenome.Cintestinalis.KY
genome <- BSgenome.Cintestinalis.KY
head(seqlengths(genome))


## ---------------------------------------------------------------------
## Genome-wide motif searching
## ---------------------------------------------------------------------
## See the GenomeSearching vignette in the BSgenome software
## package for some examples of genome-wide motif searching using
## Biostrings and the BSgenome data packages:
if (interactive())
    vignette("GenomeSearching", package="BSgenome")
}

\keyword{package}
\keyword{data}
