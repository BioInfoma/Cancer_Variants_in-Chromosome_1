

# 1) create sample lists per superpopulation from panel file
# panel columns: SampleID  Population  SuperPopulation  Gender
awk '$3=="AFR" {print $1}' integrated_call_samples_v3.20130502.ALL.panel > samples_AFR.txt
awk '$3=="EUR" {print $1}' integrated_call_samples_v3.20130502.ALL.panel > samples_EUR.txt
awk '$3=="EAS" {print $1}' integrated_call_samples_v3.20130502.ALL.panel > samples_EAS.txt
awk '$3=="AMR" {print $1}' integrated_call_samples_v3.20130502.ALL.panel > samples_AMR.txt
awk '$3=="SAS" {print $1}' integrated_call_samples_v3.20130502.ALL.panel > samples_SAS.txt

# 2) compute allele frequencies per superpopulation using vcftools
# (vcftools must be installed; if not, use bcftools/htslib approach)
VCF=ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

vcftools --gzvcf $VCF --keep samples_AFR.txt --freq --out AFR_chr1
vcftools --gzvcf $VCF --keep samples_EUR.txt --freq --out EUR_chr1
vcftools --gzvcf $VCF --keep samples_EAS.txt --freq --out EAS_chr1
vcftools --gzvcf $VCF --keep samples_AMR.txt --freq --out AMR_chr1
vcftools --gzvcf $VCF --keep samples_SAS.txt --freq --out SAS_chr1

# This produces: AFR_chr1.frq, EUR_chr1.frq, EAS_chr1.frq, AMR_chr1.frq, SAS_chr1.frq
# Format: CHROM POS N_ALLELES N_CHR {ALLELE:MAC:MAF}...
