#Run QTLtools PCA
# QTLtools pca --bed PU1.norm_prop.txt.gz --center --scale --out PU1.pheno_pca
# QTLtools pca --vcf ../../../../../results/PU1/genotypes/Waszak_47_samples.merged.vcf.gz --center --scale --out PU1.geno_pca

#Import PCA results
phenotype_pca = readr::read_delim("processed/PU1/qtltools/input/cqn/PU1.pheno_pca.pca", delim = " ")
genotype_pca = readr::read_delim("processed/PU1/qtltools/input/cqn/PU1.geno_pca.pca", delim = " ")

#Filter genotype PCA to the same set of individuals
genotype_pca = genotype_pca[,colnames(phenotype_pca)]

#Merge them together
covariates = bind_rows(phenotype_pca[1:5,], genotype_pca[1:5,])
output_path = "processed/PU1/qtltools/input/cqn/PU1.covariates_prop.txt"
write.table(covariates, output_path, sep = " ", quote = FALSE, row.names = FALSE)

