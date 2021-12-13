#wget -U firefox http://ghost.zool.kyoto-u.ac.jp/datas/JoinedScaffold.zip
#unzip JoinedScaffold.zip
#mkdir -p KH
Rscript forgeGenome.R
R CMD build BSgenome.Cintestinalis.KY --no-manual
R CMD check BSgenome.Cintestinalis.KY --no-manual
R CMD INSTALL BSgenome.Cintestinalis.KY
