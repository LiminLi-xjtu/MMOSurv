if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
clinical<-GDCquery_clinic(project="TCGA-LAML", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-LAML.csv')

clinical<-GDCquery_clinic(project="TCGA-ACC", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-ACC.csv')

clinical<-GDCquery_clinic(project="TCGA-BLCA", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-BLCA.csv')

clinical<-GDCquery_clinic(project="TCGA-LGG", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-LGG.csv')

clinical<-GDCquery_clinic(project="TCGA-BRCA", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-BRCA.csv')

clinical<-GDCquery_clinic(project="TCGA-CESC", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-CESC.csv')

clinical<-GDCquery_clinic(project="TCGA-CHOL", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-CHOL.csv')

clinical<-GDCquery_clinic(project="TCGA-LCML", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-LCML.csv')

clinical<-GDCquery_clinic(project="TCGA-COAD", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-COAD.csv')

clinical<-GDCquery_clinic(project="TCGA-ESCA", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-ESCA.csv')

clinical<-GDCquery_clinic(project="TCGA-GBM", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-GBM.csv')

clinical<-GDCquery_clinic(project="TCGA-HNSC", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-HNSC.csv')

clinical<-GDCquery_clinic(project="TCGA-KICH", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-KICH.csv')

clinical<-GDCquery_clinic(project="TCGA-KIRC", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-KIRC.csv')

clinical<-GDCquery_clinic(project="TCGA-KIRP", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-KIRP.csv')

clinical<-GDCquery_clinic(project="TCGA-LIHC", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-LIHC.csv')

clinical<-GDCquery_clinic(project="TCGA-LUAD", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-LUAD.csv')

clinical<-GDCquery_clinic(project="TCGA-LUSC", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-LUSC.csv')

clinical<-GDCquery_clinic(project="TCGA-DLBC", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-DLBC.csv')

clinical<-GDCquery_clinic(project="TCGA-MESO", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-MESO.csv')

clinical<-GDCquery_clinic(project="TCGA-OV", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-OV.csv')

clinical<-GDCquery_clinic(project="TCGA-PAAD", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-PAAD.csv')

clinical<-GDCquery_clinic(project="TCGA-PCPG", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-PCPG.csv')

clinical<-GDCquery_clinic(project="TCGA-PRAD", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-PRAD.csv')

clinical<-GDCquery_clinic(project="TCGA-READ", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-READ.csv')

clinical<-GDCquery_clinic(project="TCGA-SARC", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-SARC.csv')

clinical<-GDCquery_clinic(project="TCGA-SKCM", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-SKCM.csv')

clinical<-GDCquery_clinic(project="TCGA-STAD", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-STAD.csv')

clinical<-GDCquery_clinic(project="TCGA-TGCT", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-TGCT.csv')

clinical<-GDCquery_clinic(project="TCGA-THYM", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-THYM.csv')

clinical<-GDCquery_clinic(project="TCGA-THCA", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-THCA.csv')

clinical<-GDCquery_clinic(project="TCGA-UCS", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-UCS.csv')

clinical<-GDCquery_clinic(project="TCGA-UCEC", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-UCEC.csv')

clinical<-GDCquery_clinic(project="TCGA-UVM", type="clinical")
write.csv(clinical, '../Origin_Data/Clinical/TCGA-UVM.csv')

