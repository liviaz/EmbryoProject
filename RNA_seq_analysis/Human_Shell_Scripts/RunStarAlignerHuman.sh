#!/bin/bash
# Type sudo -s before running this script (run as root since memory clearing operation at end requires it)

BASE_DIR="/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/HumanHiSeq"
STAR_PATH="/home/livia/Desktop/RNA_seq_analysis/STAR_2.3.0e.Linux_x86_64_static/STAR"
GENOME_DIR="/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/Homo_sapiens_UCSC_hg18/Sequence/STAR"
MAX_RAM=31800000000

# remove all previous STAR directories and make fresh ones

# 1. Generate Genome (only needs to be done once)

#cd $BASE_DIR
#rm -rf "./StarOut"
#rm -rf $GENOME_DIR
#mkdir "./StarOut"
#mkdir $GENOME_DIR
#cd "./StarOut"

#$STAR_PATH --runMode genomeGenerate --genomeDir $GENOME_DIR --genomeFastaFiles $HUMAN_GENOME_FASTA --sjdbGTFfile $REF_GTF_HUMAN --genomeChrBinNbits 12 --genomeSAindexNbases 12 --runThreadN 1 --limitGenomeGenerateRAM $MAX_RAM 


# declare -a sampleNum=(2 4 5 6 7 8 9)
# declare -a sampleNum=(10 11 12 13 15)
# declare -a sampleNum=(16 17 18 19 20)
declare -a sampleNum=(1 3 14 21 22)

# clear inactive memory and caches
sync && echo 3 | tee /proc/sys/vm/drop_caches 

# loop over samples
for i in ${sampleNum[@]}
do

	FILENAME="S"$i
	echo "Now aligning sample : "$FILENAME	

	cd $BASE_DIR"/StarOut"
	rm -rf "./"$FILENAME
	mkdir "./"$FILENAME 
	cd "./"$FILENAME

	# 2. Align samples

	$STAR_PATH --runMode alignReads --genomeDir $GENOME_DIR --readFilesIn $BASE_DIR"/RawData/"$FILENAME"_R1.fastq" $BASE_DIR"/RawData/"$FILENAME"_R2.fastq" --runThreadN 4 --genomeChrBinNbits 10 --genomeLoad NoSharedMemory --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitGenomeGenerateRAM $MAX_RAM 

	# 3. Convert SAM to BAM and sort
	cd $BASE_DIR
	mkdir -p "./HTSeq_Count_Out"
	mkdir -p "./"$FILENAME"_bam"

	echo "now sorting .bam file and running HTSeq-count for sample "$FILENAME
	
	# convert sam to bam
	samtools view -b -S "./StarOut/"$FILENAME"/Aligned.out.sam" > "./"$FILENAME"_bam/"$FILENAME"_pairedonly.bam"

	# sort bam file
	samtools sort -n "./"$FILENAME"_bam/"$FILENAME"_pairedonly.bam" "./"$FILENAME"_bam/"$FILENAME"_pairedonly_sorted"
	
	# then run HTSeq-count
	python -m HTSeq.scripts.count --stranded=no -f bam "./"$FILENAME"_bam/"$FILENAME"_pairedonly_sorted.bam" $REF_GTF_HUMAN > "./HTSeq_Count_Out/"$FILENAME"_paired_count.txt"

	# clear inactive memory and caches
	sync && echo 3 | tee /proc/sys/vm/drop_caches 

done




