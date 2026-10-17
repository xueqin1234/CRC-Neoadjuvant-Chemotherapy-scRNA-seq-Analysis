library(Seurat)
setwd("./test/00data")
stromal<-readRDS(file = "stromal_umap_new.rds")
table(Idents(stromal))

stromal@meta.data$cluster_annotation <- Idents(stromal)
setwd("./test/02stromal")
pdf("04maker_stromal.pdf",width=12, height=6)
#markers.to.plot <- c("MMP2","DCN","COL1A2","TAGLN","PDGFA","APOE","CCL8","ADAMDEC1","BMP5","FRZB","POSTN","CCDC80","OGN","FBLN1","EFEMP1","ACTA2","MMP1","MMP3","RGCC","RAMP3","LYVE1","PROX1","ACKR1","SYNPO2","CRYAB","CNN1","DES","RGS5","SOX10","S100B","PLP1")
#markers.to.plot <- c("DCN","COL1A2","FAP","MYLK","MMP2","MYH11","ACTG2","APOE","CCL8","ADAMDEC1","BMP5","FRZB","POSTN","CCDC80","OGN","FBLN1","EFEMP1","RGCC","COL4A1","RAMP3","CA4","CD320","EDN1","SELP","ACKR1","VCAM1","PROX1","LYVE1","RGS5","NDUFA4L2","PDGFRB","CRYAB","S100B","PLP1","CNN1","DES")
markers.to.plot <- c("DCN","PDPN","FAP","FHL2","MYLK","MYH11","ACTG2","APOE","CCL8","ADAMDEC1","BMP5","FRZB","POSTN","CCDC80","OGN","FBLN1","EFEMP1","RGCC","COL4A1","RAMP3","CA4","CD320","EDN1","SELP","ACKR1","VCAM1","PROX1","LYVE1","RGS5","NDUFA4L2","PDGFRB","CRYAB","S100B","PLP1","CNN1","DES")

DotPlot(stromal, features = markers.to.plot,assay='RNA',group.by = 'cluster_annotation')+theme(axis.text.x=element_text(face="italic",size=8)) + theme(axis.text.y = element_text(size=12)) +RotatedAxis() + theme(axis.text.x = element_text(vjust = 1)) + coord_fixed(ratio=3/3) #调整字体大小，横坐标倾斜率
dev.off()


