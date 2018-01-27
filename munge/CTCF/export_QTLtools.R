ctcf_se = readRDS("results/SummarizedExperiment/CTCF_SummarizedExperiment.rds")

#Keep only autosomal chromosomes
valid_chr = c("1","10","11","12","13","14","15","16","17","18","19", "20", "21","22","3","4","5","6","7","8","9")
valid_rows = rowData(ctcf_se)[rowData(ctcf_se)$chr %in% valid_chr,]
filtered_se = ctcf_se[valid_rows$phenotype_id,]

#Extract data
event_metadata = rowData(filtered_se) %>% tbl_df2()

#Construct metadata
genepos = dplyr::mutate(event_metadata, transcript_id = phenotype_id, gene_id = phenotype_id) %>%
  constructQTLtoolsGenePos()
output_path = "processed/CTCF/qtltools/input/cqn/"

#Make qtltools matrix
matrix = prepareQTLtoolsMatrix(assays(filtered_se)$cqn, genepos)
saveFastqtlMatrices(list(CTCF = matrix), output_path, file_suffix = "norm_prop")
