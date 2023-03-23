#!/bin/bash

### SOFTWARE ### please change software paths to your own
export PATH="$PATH:$HOME/PROJECTS/SOFTWARE/sratoolkit.2.10.9-ubuntu64/bin"
export PATH="$PATH:/home/nbrait/PROJECTS/SOFTWARE/lofreq_star-2.1.2/bin"
export PATH="$PATH:/home/nbrait/TOOLS/bowtie2-2.3.5.1"
export LD_LIBRARY_PATH=/home/nbrait/TOOLS/htslib-1.9/
### PROJECT ###
R_DIR="/mnt/c/Users/nadja/Documents/LaptopAsus/PhD/Quaranjaproject/Aedes/diamond_blastx/Aedes/diamond_blast/EVE_ANALYSIS/R_folder/coverage_analysis"

### PATHS ###
READ_DIR="/mnt/c/Users/nadja/Documents/LaptopAsus/PhD/Quaranjaproject/COVERAGE/unmapped_reads" # unmapped reads from pre-processing workflow
PROC_DIR="/mnt/c/Users/nadja/Documents/LaptopAsus/PhD/Quaranjaproject/COVERAGE" # project directory
### INPUT ###
INPUT=$PROC_DIR/"Culicinae_accessions.txt"
###############################################################################################################################
cd $PROC_DIR
mkdir -p PROCESSED/{DEPTH_VARIANTS/PLOTS,CORRECTION} REF

# continuing workflow for every single SRR sample
while IFS="" read -u 8 SRR || [ -n "$SRR" ]
do
# making individual reference lists for each SRR sample
	cd $PROC_DIR/REF
	mkdir $SRR
	grep -A 1 $SRR $PROC_DIR/Culicinae_consensus_sequences.fasta --no-group-separator > ${SRR}/reference.fa

# Sample reference Indexing
	if [ ! -e $PROC_DIR/REF${SRR}/reference-index.1.bt2 ];then
        	bowtie2-build --threads 35 -f --quiet $PROC_DIR/REF/${SRR}/reference.fa $PROC_DIR/REF/${SRR}/reference-index
        fi
#MAPPING WITH BOWTIE2 for paired end reads
	cd $READ_DIR

	if [ ! -e $READ_DIR/${SRR}_1.unmapped.fastq ] || [ ! -e $READ_DIR/${SRR}_2.unmapped.fastq ] && [ ! -e $READ_DIR/${SRR}_single.unmapped.fastq ];then #either paired or single file 
                echo "$SRR: sample not there: $(date) " | tee -a $PROC_DIR/log_file
		n=$(awk -v x="$SRR" '$0~x {print NR}' $INPUT)
        	echo "--------------------------> $n analysed: $(date) " | tee -a $PROC_DIR/log_file
                continue
        fi
	echo "${SRR} mapping data:" | tee -a $PROC_DIR/map_log_file.txt
	if [ -e ${SRR}_1.unmapped.fastq ] && [ -e ${SRR}_2.unmapped.fastq ];then
		bowtie2 -p 15 --end-to-end -x $PROC_DIR/REF/${SRR}/reference-index \
		--fr -1 ${SRR}_1.unmapped.fastq -2 ${SRR}_2.unmapped.fastq -U \
		${SRR}_unpaired.unmapped.fastq 1> $PROC_DIR/${SRR}.sam 2>> $PROC_DIR/map_log_file.txt
#MAPPING WITH BOWTIE2 for single end reads
	else
		bowtie2 -p 15 --end-to-end -x $PROC_DIR/REF/${SRR}/reference-index ${SRR}_single.unmapped.fastq \
		1> $PROC_DIR/${SRR}.sam 2>> $PROC_DIR/map_log_file.txt
	fi
	cd $PROC_DIR
	if [ -s ${SRR}.sam ];then
		echo "$SRR: Mapped succesfully to target: $(date) " | tee -a $PROC_DIR/map_log_file.txt
	else
		echo "$SRR: no mapped sequences: $(date)" | tee -a $PROC_DIR/map_log_file.txt
	fi
# converting/sorting and coverage analysis
	if [ -f "${SRR}.sam" ];then
		if [ -e $READ_DIR/${SRR}_1.unmapped.fastq ] && [ -e $READ_DIR/${SRR}_2.unmapped.fastq ];then
			samtools view -bSh -q 2 -f 2 -o ${SRR}.bam ${SRR}.sam
	        else
			samtools view -bSh -q 2 -o ${SRR}.bam ${SRR}.sam
      		fi
		samtools sort ${SRR}.bam -o ${SRR}_sort.bam
		mkdir DIR_$SRR
		mv ${SRR}* DIR_${SRR}
		cd DIR_${SRR}
		bamtools split -in ${SRR}_sort.bam -reference
		samtools coverage ${SRR}_sort.REF_*.bam >> ${SRR}_coverage.txt
	fi
# processing reference sequences for correction later on
        cd $PROC_DIR/REF/$SRR
        if [ -f "reference.fa" ];then
               awk '{ sub("\r$", ""); print }' reference.fa > ref_linux.fa
               grep ">" ref_linux.fa | perl -pe 's|>SRR.*?_||' > ../${SRR}_refs.txt
               awk '/^>/ {OUT=substr($0,2) ".fa"}; OUT {print >OUT}' ref_linux.fa
               cd ..
# processing individual segments
		INPUT2="${SRR}_refs.txt" # list of all the segments per sample
		mkdir $PROC_DIR/PROCESSED/${SRR}
		mkdir $PROC_DIR/PROCESSED/DEPTH_VARIANTS/${SRR}
		while IFS= read -u 8 CHIM
		do
			cd $PROC_DIR/DIR_${SRR}
			if  [ -e ${SRR}_sort.REF_${SRR}_${CHIM}.bam ];then
				echo "worked"
				mv ${SRR}_sort.REF_${SRR}_${CHIM}.bam ${SRR}_${CHIM}_sort.bam
			fi
			samtools depth -a ${SRR}_${CHIM}_sort.bam -o ${SRR}_${CHIM}_cov.txt
			if [ -e ${SRR}_${CHIM}_sort.bam ];then
# processing steps for coverage analysis and variant calling
				TARGET=${SRR}_${CHIM}
				samtools index ${TARGET}_sort.bam
				bedtools genomecov -d -ibam ${TARGET}_sort.bam | grep -w "${TARGET}" > ${TARGET}.cov
				cp ${TARGET}.cov $PROC_DIR/PROCESSED/DEPTH_VARIANTS
				cd $PROC_DIR/REF/${SRR}
				samtools faidx ${TARGET}.fa
				lofreq call -f ${TARGET}.fa -o $PROC_DIR/DIR_${SRR}/${TARGET}.vcf $PROC_DIR/DIR_${SRR}/${TARGET}_sort.bam
				cd $PROC_DIR/DIR_${SRR}
				cp ${TARGET}.cov $R_DIR
				cp ${TARGET}.vcf $R_DIR
#  R-scripts for both coverage and variant calling - generation of ggplots
				cd $R_DIR
				Rscript $R_DIR/cov_analysis.r 2>> $PROC_DIR/PROCESSED/DEPTH_VARIANTS/R_log_file.txt
				rm ${TARGET}.vcf
				cp $PROC_DIR/DIR_${SRR}/{*.cov,*.vcf} $PROC_DIR/PROCESSED/CORRECTION
				cp $PROC_DIR/REF/${SRR}/*.fa $PROC_DIR/PROCESSED/CORRECTION
                                cd ${PROC_DIR}/PROCESSED/DEPTH_VARIANTS/PLOTS
				mv ${TARGET}.cov ${TARGET}_cov.png
				if [ -e ${TARGET} ];then
	                               	mv ${TARGET} ${TARGET}_vcf.png #if statement better because there isn't always a vcf plot
				fi
			fi
	cd $PROC_DIR/DIR_${SRR}
	ls *.cov | sed 's/.cov//' >> $PROC_DIR/PROCESSED/CORRECTION/listsample.txt
# R-script for correction of sequences
	Rscript $R_DIR/correction_analysis.r 2>> $PROC_DIR/PROCESSED/CORRECTION/R_log_file_cor.txt # cutoff is 3
	Rscript $R_DIR/correction1_analysis.r 2>> $PROC_DIR/PROCESSED/CORRECTION/R_log_file_cor.txt # cutoff is 1 - better for the EVEs
	rm $PROC_DIR/PROCESSED/CORRECTION/listsample.txt
	cat $PROC_DIR/PROCESSED/CORRECTION/*3_corrected.fasta > $PROC_DIR/PROCESSED/CORRECTION/all_corrected_3.fasta
	mkdir $PROC_DIR/PROCESSED/DEPTH_VARIANTS/DIR_${SRR}
	mv $PROC_DIR/PROCESSED/DEPTH_VARIANTS/${SRR}* $PROC_DIR/PROCESSED/DEPTH_VARIANTS/DIR_${SRR}
		done 8< "$INPUT2"
	fi
	rm $PROC_DIR/DIR_${SRR}/*sam
done 8<"$INPUT"
# Filtering EVE plots
		cd ${PROC_DIR}/PROCESSED/DEPTH_VARIANTS
		awk '{ sub("\r$", ""); print }' EVES.txt > EVES_ID.txt
		INPUT2="EVES_ID.txt"
		while IFS="" read -u 9 EVE || [ -n "$EVE" ]
		do
			cd ${PROC_DIR}/PROCESSED/DEPTH_VARIANTS/PLOTS
			mv ${EVE}* ${PROC_DIR}/PROCESSED/DEPTH_VARIANTS/PLOTS/EVES
			cat ${EVE}*.cov >> new_EVE_cov.txt
		done 9< "$INPUT2"
#progress:
	n=$(awk -v x="$SRR" '$0~x {print NR}' $INPUT)
	echo "--------------------------> $n samples analysed: $(date)"

done 8< "$INPUT"
