## Shell script to get Kelch genotypes from our and Jane's genomes
## Rifat has already split and cleaned up the individual genomes for his pfs project
## Script started 13 December 2015

vcfdir=/proj/julianog/users/RifatR
kelchdir=/proj/julianog/users/ChristianP/PvPf_Kelch
ref=/proj/julianog/refs/PvSAL1_v13.0/PlasmoDB-13.0_PvivaxSal1_Genome.fasta
gatk=/nas02/apps/biojars-1.0/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar

for name in `cat $vcfdir/names/all_goods_5x@80%.txt`
do

	echo $name

	## VCF TO FASTA
	java -Xmx2g -jar $gatk \
		-T FastaAlternateReferenceMaker \
		-R $ref \
		--variant $vcfdir/variants/indivs/$name.sans0.vcf \
		-o $kelchdir/global/$name.fa \
		-L Pv_Sal1_chr12:447729-449867

		## CREATING THE MULTIFASTA
		echo ">"$name >> $kelchdir/globalGenotypes.fa
		tail -n +2 $kelchdir/global/$name.fa >> $kelchdir/globalGenotypes.fa
			# The tail command removes the fasta header line

	done
done
