#BSUB -W 36:00
#BSUB -M 150000
#BSUB -n 75
#BSUB -J snakeATAC_kottyan_atopicDerm
#BSUB -o %J.out
#BSUB -e %J.err

module load anaconda3

source activate snakemake

cd /data/miraldiNB/Tareian/snakemake/snakeATAC

snakemake --cores 75 --use-conda --conda-frontend conda --rerun-incomplete --unlock
snakemake --cores 75 --use-conda --conda-frontend conda --rerun-incomplete 