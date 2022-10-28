

#make cellranger ref genome
cd /ref_genomes
wget https://ftp.ensembl.org/pub/release-108/gtf/mus_musculus/Mus_musculus.GRCm39.108.gtf.gz
gunzip Mus_musculus.GRCm39.108.gtf.gz

wget https://ftp.ensembl.org/pub/release-108/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

module load cellranger
cellranger mkgtf Mus_musculus.GRCm39.108.gtf Mus_musculus.GRCm39.108.filtered.gtf \
                 --attribute=gene_biotype:protein_coding \
                 --attribute=gene_biotype:lincRNA \
                 --attribute=gene_biotype:antisense \
                 --attribute=gene_biotype:IG_LV_gene \
                 --attribute=gene_biotype:IG_V_gene \
                 --attribute=gene_biotype:IG_V_pseudogene \
                 --attribute=gene_biotype:IG_D_gene \
                 --attribute=gene_biotype:IG_J_gene \
                 --attribute=gene_biotype:IG_J_pseudogene \
                 --attribute=gene_biotype:IG_C_gene \
                 --attribute=gene_biotype:IG_C_pseudogene \
                 --attribute=gene_biotype:TR_V_gene \
                 --attribute=gene_biotype:TR_V_pseudogene \
                 --attribute=gene_biotype:TR_D_gene \
                 --attribute=gene_biotype:TR_J_gene \
                 --attribute=gene_biotype:TR_J_pseudogene \
                 --attribute=gene_biotype:TR_C_gene

cellranger mkref --genome=mm39 \
                 --fasta=Mus_musculus.GRCm39.dna.primary_assembly.fa \
                 --genes=Mus_musculus.GRCm39.108.filtered.gtf 
                 
cellranger mkvdjref --genome=mm39_vdj \
                    --fasta=Mus_musculus.GRCm39.dna.primary_assembly.fa \
                 --genes=Mus_musculus.GRCm39.108.filtered.gtf 

#download data

#cellranger
module load cellranger
# if getting error, try cellranger/6.0.0
cd /data/cellranger
cellranger multi --id=T1 --csv=../configs/T1.csv
cellranger multi --id=T2 --csv=../configs/T2.csv
cellranger multi --id=T3 --csv=../configs/T3.csv
cellranger multi --id=T4 --csv=../configs/T4.csv
cellranger multi --id=T7 --csv=../configs/T7.csv
cellranger multi --id=T8 --csv=../configs/T8.csv
cellranger multi --id=T11 --csv=../configs/T11.csv
cellranger multi --id=T12 --csv=../configs/T12.csv
cellranger multi --id=T34 --csv=../configs/T34.csv
cellranger multi --id=T36 --csv=../configs/T36.csv

