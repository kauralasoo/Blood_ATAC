#Process Blood ATAC data
#Dry run
snakemake -np processed/out.txt
#Cluster run
snakemake --cluster scripts/snakemake_submit.py -p processed/out.txt --jobs 100

#Process Macrophage ATAC data
#Dry run
snakemake -np macrophages/out.txt --configfile config_macrophages.yaml
#Cluster run
snakemake --cluster scripts/snakemake_submit.py -p macrophages/out.txt --jobs 100 --configfile config_macrophages.yaml



#Process PU.1 data
#Convert bams to fastq files
snakemake --cluster scripts/snakemake_submit_UT.py -np processed/PU1/out.txt -s bam_to_fastq_SE.snakefile --configfile configs/config_PU1.yaml --jobs 5

#Run the processing pipeline
snakemake --cluster scripts/snakemake_submit_UT.py -np processed/PU1/out.txt -s ChIP_pipeline_SE.snakefile --configfile configs/config_PU1.yaml --jobs 5
