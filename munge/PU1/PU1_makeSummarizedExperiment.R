library("rtracklayer")
library("dplyr")
library("cqn")
library("devtools")
library("SummarizedExperiment")
load_all("../../seqUtils/")

#Import peaks
peaks = import.gff("results/PU1/PU1.peaks.gtf") %>%
  as.data.frame() %>% tbl_df() %>%
  dplyr::transmute(phenotype_id = gene_id, chr = as.character(seqnames), start, end, length = width) %>%
  dplyr::mutate(strand = 1)

#Import GC-content
gc_content = readr::read_delim("results/PU1/PU1.peaks.GCcontent.txt", delim = "\t", col_types = "ccciiccccddiiiiiii")
gc_content = gc_content[,c(1,4,5,11)]
colnames(gc_content) = c("chr","start","end", "percentage_gc_content")

#Merge
peak_data = dplyr::left_join(peaks, gc_content, by = c("chr","start","end")) %>%
  as.data.frame()
rownames(peak_data) = peak_data$phenotype_id

#Import counts
counts = read.table("results/PU1/PU1.counts.merged.txt", header = TRUE)
count_matrix = dplyr::select(counts, -gene_id, -length) %>% as.matrix()
rownames(count_matrix) = counts$gene_id

#Normalise with CQN
cqn_matrix = calculateCQN(count_matrix, dplyr::rename(peak_data, gene_id = phenotype_id))

#Import the list of genotyped samples
genotyped = read.table("data/Waszak_2015/Waszak_2015_genotype_list.txt")

#Import 1000G metadata
kg_meta = readr::read_tsv("../../RNAseq_pipeline/metadata/GEUVADIS/1000_genomes_sample_metadata.tsv", col_names =
                            c("genotype_id", "sex", "biosample_id", "population_code", "population_name", 
                              "superpopulation_code", "superpopulation_name", "data_collections"), skip = 1) %>%
  dplyr::select(-data_collections)

#Import Waszak metadata
sample_meta = readr::read_tsv("data/Waszak_2015/Waszak_2015_clean_metadata.txt") %>%
  dplyr::filter(sample_id %in% colnames(count_matrix)) %>%
  dplyr::select(sample_id, genotype_id) %>%
  dplyr::left_join(kg_meta, by = "genotype_id") %>%
  as.data.frame()
rownames(sample_meta) = sample_meta$sample_id


#Construct a SummarizedExperiment object
se = SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = count_matrix[,rownames(sample_meta)], cqn = cqn_matrix[,rownames(sample_meta)]), 
  colData = sample_meta, 
  rowData = peak_data)
saveRDS(se, "results/SummarizedExperiment/PU1_SummarizedExperiment.rds")


