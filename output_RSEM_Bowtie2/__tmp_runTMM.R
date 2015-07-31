library(edgeR)

rnaseqMatrix = read.table("Trinity_isoform.not_cross_norm.fpkm.tmp", header=T, row.names=1, com='', check.names=F)
rnaseqMatrix = round(rnaseqMatrix)
exp_study = DGEList(counts=rnaseqMatrix, group=factor(colnames(rnaseqMatrix)))
exp_study = calcNormFactors(exp_study)
exp_study$samples$eff.lib.size = exp_study$samples$lib.size * exp_study$samples$norm.factors
write.table(exp_study$samples, file="Trinity_isoform.not_cross_norm.fpkm.tmp.TMM_info.txt", quote=F, sep="\t", row.names=F)
