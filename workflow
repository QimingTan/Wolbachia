#01 Quality control filtering of raw data (Group member)
5,551 mosquito transcriptome samples collected from various regions across China
fastq/QC

#02 Species identification of mosquito samples (Group member)
Assemble contigs based on the mitochondrial COI gene, perform BLASTn alignment against the database, and identify mosquito species
A total of 10 mosquito genera were identified, with a small number of samples (200) unable to be assigned to a specific genus; 
The majority of these belong to the four major mosquito genera: Cluex,Aedes,Armigeres,Anopheles

#03 Remove host sequences
1. Download genome data of various mosquito species from NCBI and build blast indexes.

2.Use Bowtie2 to align sequences and select the unaligned sequences
bowtie2 -x $index -1 $raw_fq1 -2 $raw_fq2 | samtools fastq -f 13 -1 $clean_fq1 -2 $clean_fq2

#04 Identify Wolbachia-positive samples.
1.Identify single-copy orthologous genes of Wolbachia strain A based on BUSCO (362/364)
busco  -i WolA_genomic.fna  -l rickettsiales_odb10 -o output  -m genome
2.Align the nohost sequences to the reference genome
bowtie2 -x ${bacindex} -1 ${fq1} -2 ${fq2} --very-sensitive-local -p 4 --rg-id \"$sn\" --rg \"SM:$sn\" | samtools view -bS -F 4 | samtools sort -o ${sn}.bam

3.Select samples that contain more than one gene and have a coverage rate of over 50%
samtools coverage -H -q 10 ${sn}.bam
cut -f1 all_50.coverage | sort | uniq > unique_items.txt #3146 samples
awk '{print $1}' all_50.coverage | sort | uniq -c | sort -nr | awk '$1 > 100'|wc	#133 samples 

#05 Perform freebayes analysis based on single-copy genes to obtain VCF, and extract SNP sites from the VCF
freebayes -f Wolbachia_sco.fna -L wolbachia.list --standard-filters  > Wolbachia.vcf
vcftools --vcf $vcf --max-missing 0.4 --minQ 10 --remove-filtered-all --recode --recode-INFO-all --stdout | vcfsnps > snp.vcf #Filters out sites with a missing rate greater than 40%，
#06 Select missing sites and remove unreliable samples
vcftools --vcf snp.vcf --missing-indv
awk '$5 > 0.4' out.imiss | cut -f1 > lowDP.indv # Remove samples with more than 40% missing genotypes.
vcftools --vcf snp.vcf --remove lowDP.indv --recode --recode-INFO-all --stdout > freebayes_snp.vcf
bcftools stats Wolbachia_freebayes_p1_snp.vcf #samples:989, SNPs:1345

#06 Plink-PCA, Evaluate genetic variation within and between populations, and identify and correct population structure.
plink --allow-extra-chr --vcf freebayes_snp.vcf --recode 
plink --allow-extra-chr --file $name --noweb --make-bed 
plink --allow-extra-chr --threads 4 --bfile $name --pca 20 




