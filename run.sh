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
