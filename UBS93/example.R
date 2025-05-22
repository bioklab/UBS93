####*UBS93 predicting####
source('UBS93.package.R')
load("CCLE_test_log2tpm.RData")
load("UBS93.data.RData")
CCLE_brca_subtype                   <- breast.cancer.predictor(expr.of.sample = CCLE_test_log2tpm,
                                                               expr.of.centroid = UBS93.data$UBS93.centroid,
                                                               marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                               HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                               ssgsea.cutoff = UBS93.data$ssgsea.cutoff)
CCLE_brca.subtype                   <- CCLE_brca_subtype$subtype
