#Exploring Population-Specific Cancer Variants on Chromosome 1 Using the 1000 Genomes Project

# download STEP 1: DOWNLOAD DATA & METADATA
# Download Chromosome 1 VCF (complete and correct file)
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

# Download index file for the VCF (so you don’t have to create it)
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi

# Download sample metadata
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20130502.phase3.analysis.sequence.index

# Download panel file
wget -c https://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel


#EXTRACT ONLY CANCER VARIANTS



bcftools view -R fixed_chr1_nochr.bed ALL.chr1*.vcf.gz -Oz -o chr1_cancer.vcf.gz

bcftools norm -m -any chr1_cancer.vcf.gz -Oz -o chr1_cancer.biallelic.vcf.gz

tabix -p vcf chr1_cancer.biallelic.vcf.gz


#CONVERT TO PLINK FORMAT + PCA
plink2 --vcf chr1_cancer.biallelic.vcf.gz --make-bed --out chr1_cancer --threads 4

plink2 --bfile chr1_cancer --pca --out pca_chr1

#EXTRACT GENOTYPES INTO TABLE
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' chr1_cancer.vcf.gz > genotypes.tsv

# good — valid REF and ALT for VEP
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' chr1_cancer.biallelic.vcf.gz > vep_input.tsv


#Downstream analysis in R

