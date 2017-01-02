library("tidyr")
library("dplyr")

##### Blood progenitors dataset #####
#Import metadata
metadata = read.table("Blood_ATAC/data/sample_metadata.GEO.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#Import runinfo
runinfo = read.table("Blood_ATAC/data/runinfo/blood_progentiors.csv", sep =",",header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::select(Run, Experiment) %>% 
  dplyr::filter(Run != "Run") %>% 
  dplyr::rename(SRA_accession = Experiment, run_id = Run) %>%
  tbl_df()

#Parse columns
new_metadata = tidyr::separate(metadata, donor, c("d1","donor"), sep = "donorid: ") %>% 
  dplyr::select(-d1) %>%
  tidyr::separate(cell_description, c("d1","cell_description"), sep = "cell type: ") %>%
  dplyr::select(-d1) %>%
  dplyr::mutate(SRA_accession = gsub("SRA: https://www.ncbi.nlm.nih.gov/sra?term=", "", SRA_accession, fixed = TRUE)) %>%
  dplyr::group_by(cell_type, donor) %>%
  dplyr::mutate(sample_id = paste(cell_type, donor, seq(1:length(cell_type)), sep = "_")) %>%
  dplyr::ungroup() %>%
  dplyr::select(sample_id, everything()) %>%
  dplyr::left_join(runinfo, by = "SRA_accession")

#Save new metadata to disk
write.table(new_metadata, "Blood_ATAC/data/sample_metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#### Neutrophil dataset ####
#Import metadata
metadata = read.table("Blood_ATAC/data/neutrophils_sample_metadata.GEO.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

runinfo = read.table("Blood_ATAC/data/runinfo/neutrophils.csv", sep =",",header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::select(Run, Experiment) %>% 
  dplyr::filter(Run != "Run") %>% 
  dplyr::rename(SRA_accession = Experiment, run_id = Run) %>%
  tbl_df()

new_metadata = dplyr::left_join(metadata, runinfo, by = "SRA_accession")
write.table(new_metadata, "Blood_ATAC/data/neutrophils_sample_metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE)

