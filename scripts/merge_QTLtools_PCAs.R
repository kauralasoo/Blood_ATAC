#Run QTLtools PCA
# QTLtools pca --bed CTCF.norm_prop.txt.gz --center --scale --out CTCF.pheno_pca
# QTLtools pca --vcf /gpfs/hpchome/a72094/rocket/datasets/CTCF/genotypes/vcf/GRCh38/CTCF_51_samples.GRCh38.final.vcf.gz --center --scale --out CTCF.geno_pca

#Import PCA results
phenotype_pca = readr::read_delim("processed/CTCF/qtltools/input/cqn/CTCF.pheno_pca.pca", delim = " ")
genotype_pca = readr::read_delim("processed/CTCF/qtltools/input/cqn/CTCF.geno_pca.pca", delim = " ")

#Filter genotype PCA to the same set of individuals
genotype_pca = genotype_pca[,colnames(phenotype_pca)]

#Merge them together
covariates = bind_rows(phenotype_pca[1:5,], genotype_pca[1:5,])
output_path = "processed/CTCF/qtltools/input/cqn/CTCF.covariates_prop.txt"
write.table(covariates, output_path, sep = " ", quote = FALSE, row.names = FALSE)

