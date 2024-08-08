#! /bin/bash
# Parameters for slurm (don't remove the # in front of #SBATCH!)
#  Use partition debug:
#SBATCH --partition=shortterm
#  Use one node:
#  Request 10 cores (hard constraint):
#SBATCH -c 30
#  Request 10GB of memory (hard constraint):
#SBATCH --mem=100GB
#  Request one hour maximal execution time (hard constraint):
#SBATCH --time=1-0:0:0
#  Request 100 GB of local scratch disk (hard constraint):
#SBATCH --tmp=500G
#  Notify me at job start and end:
#SBATCH --mail-type=ALL
#  Send the notifications to:
#SBATCH --mail-user=
#  Find your job easier with a name:
#SBATCH --job-name=NCRF

#Initialize the module system:
source /etc/profile.d/modules.sh
# Allow aliases (required by some modules):
shopt -s expand_aliases
# Load your necessary modules (example):
module load filtlong
module load singularity-sylabs/v3.5.3
module load ncrf/v1.01.02
module load minimap2/v2.22
module load samtools/1.15
module load nanostat/1.5.0



# Define the output folder, home directory and path to sequencing data in .zip format: 


export JOB_OUTPUT_DIR="/home/projectname/" 

SEQ_DATA=/data/youdatafolder
mkdir $SCRATCH/inputdata
mkdir $CONFIG
mkdir $SCRATCH/unzip/
mkdir --parents $JOB_OUTPUT_DIR


# Prepare sequencing data and reference file:
cp -a $HOME_DIR/FGF14.fa $SCRATCH/inputdata/
unzip $SEQ_DATA/$1.zip -d $SCRATCH/unzip/
cp -a $SCRATCH/unzip/*/*/*/fastq_pass/* $SCRATCH/inputdata/



ls $SCRATCH/inputdata
RUN_PATH=$SCRATCH/inputdata
cd $RUN_PATH
for d in $(ls $RUN_PATH)
do
    SAMPLE=`basename $d`
	cd $d
	cd $RUN_PATH
	gunzip $d/*.fastq.gz
	ls $d
	cd $d
	cat *.fastq > all.fastq

	# Read length and quality of the reads
	NanoStat  --fastq all.fastq > $JOB_OUTPUT_DIR/$d/NanoStat_$d.txt
	
	# Alignment to the reference
	minimap2 -a -x map-ont $SCRATCH/inputdata/FGF14.fa all.fastq > all.sam
	samtools view -S -b all.sam > all.bam
	samtools sort all.bam -o all.sort.bam
	samtools view -h -F 4 all.sort.bam | perl -lane '$l = 0; $F[5] =~ s/(\d+)[MX=DN]/$l+=$1/eg; print if $l > 0 or /^@/' | samtools view -bS - > filt_$d.sorted.bam
	
	samtools sort filt_$d.sorted.bam -o filt_$d.sorted.bam
	samtools index filt_$d.sorted.bam
	
	
	# Define repeat motif and NCRF input folder
	MOTIF1=GAA
	CONFIG=$JOB_OUTPUT_DIR/$d/
	IDENTIFIER=$d
	mkdir $JOB_OUTPUT_DIR/$d/
	
	
	
	# Convert the fastq to fasta
	cat all.fastq | paste - - - - | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > $SCRATCH/fastq/$d.fasta
	cp -a $SCRATCH/fastq/$d.fasta $CONFIG

	# Detect the repeat length
	cat $SCRATCH/fastq/$d.fasta \
	| singularity exec $NCRF_CONTAINER NCRF $MOTIF1 --stats=events --positionalevents --maxnoise=80%  --minlength=108 \
	| singularity exec $NCRF_CONTAINER ncrf_sort.py --sortby=mratio  \
	| tee $CONFIG/"${MOTIF1}_raw_${IDENTIFIER}.summary" \
	| singularity exec $NCRF_CONTAINER ncrf_consensus_filter.py \
	| singularity exec $NCRF_CONTAINER ncrf_sort.py --sortby=mratio  \
	| tee $CONFIG/"${MOTIF1}_refined_${IDENTIFIER}.summary" \
	| singularity exec $NCRF_CONTAINER ncrf_summary.py \
	> $CONFIG/"${MOTIF1}_summary_${IDENTIFIER}.summary"

	echo "--->NCRF-MOTIF1=${MOTIF1}-Run concluded"

	# Part for the interruptions

	cp -a $JOB_OUTPUT_DIR/$d/GAA_refined_$d.summary $SCRATCH/inputdata/
	INPUT=$SCRATCH/inputdata/GAA_refined_$d.summary
	cd $SCRATCH/inputdata/


	# Extract important lines (for + positions)
	grep "^GAA+" -B 2 $INPUT > seq.txt
	echo "Extracted lines +"

	# Extract start and stop lines
	grep -B 1 "^GAA+" seq.txt | awk -F '\n' 'ln ~ /^$/ { ln = "matched"; print $1 } $1 ~ /^--$/ { ln = "" }' > start_stop.txt
	echo "start stop extraction +"

	# Extract start and Stop Position
	cat start_stop.txt | sed -e 's/^.* .* .*bp //g' -e 's/ .*$//g' > start_stop_a.txt
	cat start_stop_a.txt | sed -e 's/-/\t/g' > pos.tsv
	echo "start stop position extraction +"

	# Extract read name
	cat start_stop.txt |awk '{print $1}' > read.txt
	cat read.txt | sed -e 's/ /\t/g' > read.tsv
 
	# Extract Insertions and Mismatches
	grep "^GAA+" seq.txt > in_mm.txt
	cat in_mm.txt | sed -e 's/^GAA+                                   .*bp score=.* //g' > in_mm_a.txt
	cat in_mm_a.txt | perl -nle '$pos=-1;while (($off=index($_,"-",$pos))>=0) {print $off;$pos=$off+1;}; print "end of string"' > in_pos.txt

	cat in_mm_a.txt | perl -nle '$pos=-1;while (($off=index($_,"a",$pos))>=0) {print $off;$pos=$off+1;}; print "end of string"' > mm_a_pos.txt
	cat in_mm_a.txt | perl -nle '$pos=-1;while (($off=index($_,"g",$pos))>=0) {print $off;$pos=$off+1;}; print "end of string"' > mm_g_pos.txt
	cat in_mm_a.txt | perl -nle '$pos=-1;while (($off=index($_,"t",$pos))>=0) {print $off;$pos=$off+1;}; print "end of string"' > mm_t_pos.txt
	cat in_mm_a.txt | perl -nle '$pos=-1;while (($off=index($_,"c",$pos))>=0) {print $off;$pos=$off+1;}; print "end of string"' > mm_c_pos.txt

	cat in_pos.txt | sed -e 's/$/a/g' > in_pos_a.txt
	tr --delete '\n' < in_pos_a.txt > in_pos_b.txt
	cat in_pos_b.txt | sed -e 's/a/ /g' > in_pos_c.txt
	cat in_pos_c.txt | sed -e 's/end of string end of string/end of string 0 end of string/g' > in_pos_d.txt
	cat in_pos_d.txt | sed -e 's/end of string/\n/g' > in_pos_e.txt
	cat in_pos_e.txt | sed -e 's/^ //g' > in_pos_f.txt
	cat in_pos_f.txt | awk '!NF{$0="0"}1' > in_pos_g.txt
	cat in_pos_g.txt | sed -e 's/ /\t/g' > in_pos.tsv
	paste read.tsv pos.tsv in_pos.tsv > $JOB_OUTPUT_DIR/$d/ins_+_$d.tsv
	echo "inn and mm position extraction +"
	
	cat mm_a_pos.txt | sed -e 's/$/a/g' > mm_a_pos_a.txt
	tr --delete '\n' < mm_a_pos_a.txt > mm_a_pos_b.txt
	cat mm_a_pos_b.txt | sed -e 's/a/ /g' >mm_a_pos_c.txt
	cat mm_a_pos_c.txt | sed -e 's/end of string end of string/end of string 0 end of string/g' > mm_a_pos_d.txt
	cat mm_a_pos_d.txt | sed -e 's/end of string/\n/g' > mm_a_pos_e.txt
	cat mm_a_pos_e.txt | sed -e 's/^ //g' > mm_a_pos_f.txt
	cat mm_a_pos_f.txt | awk '!NF{$0="0"}1' > mm_a_pos_g.txt
	cat mm_a_pos_g.txt | sed -e 's/ /\t/g' > mm_a_pos.tsv
	paste read.tsv pos.tsv mm_a_pos.tsv > $JOB_OUTPUT_DIR/$d/mm_a_+_$d.tsv
	echo "inn and mm position extraction +"
	
	cat mm_g_pos.txt | sed -e 's/$/a/g' > mm_g_pos_a.txt
	tr --delete '\n' < mm_g_pos_a.txt > mm_g_pos_b.txt
	cat mm_g_pos_b.txt | sed -e 's/a/ /g' >mm_g_pos_c.txt
	cat mm_g_pos_c.txt | sed -e 's/end of string end of string/end of string 0 end of string/g' > mm_g_pos_d.txt
	cat mm_g_pos_d.txt | sed -e 's/end of string/\n/g' > mm_g_pos_e.txt
	cat mm_g_pos_e.txt | sed -e 's/^ //g' > mm_g_pos_f.txt
	cat mm_g_pos_f.txt | awk '!NF{$0="0"}1' > mm_g_pos_g.txt
	cat mm_g_pos_g.txt | sed -e 's/ /\t/g' > mm_g_pos.tsv
	paste read.tsv pos.tsv mm_g_pos.tsv > $JOB_OUTPUT_DIR/$d/mm_g_+_$d.tsv
	echo "inn and mm position extraction +"
	
	cat mm_c_pos.txt | sed -e 's/$/a/g' > mm_c_pos_a.txt
	tr --delete '\n' < mm_c_pos_a.txt > mm_c_pos_b.txt
	cat mm_c_pos_b.txt | sed -e 's/a/ /g' >mm_c_pos_c.txt
	cat mm_c_pos_c.txt | sed -e 's/end of string end of string/end of string 0 end of string/g' > mm_c_pos_d.txt
	cat mm_c_pos_d.txt | sed -e 's/end of string/\n/g' > mm_c_pos_e.txt
	cat mm_c_pos_e.txt | sed -e 's/^ //g' > mm_c_pos_f.txt
	cat mm_c_pos_f.txt | awk '!NF{$0="0"}1' > mm_c_pos_g.txt
	cat mm_c_pos_g.txt | sed -e 's/ /\t/g' > mm_c_pos.tsv
	paste read.tsv pos.tsv mm_c_pos.tsv > $JOB_OUTPUT_DIR/$d/mm_c_+_$d.tsv
	echo "inn and mm position extraction +"
	
	cat mm_t_pos.txt | sed -e 's/$/a/g' > mm_t_pos_a.txt
	tr --delete '\n' < mm_t_pos_a.txt > mm_t_pos_b.txt
	cat mm_t_pos_b.txt | sed -e 's/a/ /g' >mm_t_pos_c.txt
	cat mm_t_pos_c.txt | sed -e 's/end of string end of string/end of string 0 end of string/g' > mm_t_pos_d.txt
	cat mm_t_pos_d.txt | sed -e 's/end of string/\n/g' > mm_t_pos_e.txt
	cat mm_t_pos_e.txt | sed -e 's/^ //g' > mm_t_pos_f.txt
	cat mm_t_pos_f.txt | awk '!NF{$0="0"}1' > mm_t_pos_g.txt
	cat mm_t_pos_g.txt | sed -e 's/ /\t/g' > mm_t_pos.tsv
	paste read.tsv pos.tsv mm_t_pos.tsv > $JOB_OUTPUT_DIR/$d/mm_t_+_$d.tsv
	echo "inn and mm position extraction +"

	# Extract Deletions 
	cat start_stop.txt | sed -e 's/^.* .* .*bp .*-.* //g' > del_a.txt
	cat del_a.txt | perl -nle '$pos=-1;while (($off=index($_,"-",$pos))>=0) {print $off;$pos=$off+1;}; print "end of string"' > del_pos.txt

	cat del_pos.txt | sed -e 's/$/a/g' > del_pos_a.txt
	tr --delete '\n' < del_pos_a.txt > del_pos_b.txt
	cat del_pos_b.txt | sed -e 's/a/ /g' > del_pos_c.txt
	cat del_pos_c.txt | sed -e 's/end of string end of string/end of string 0 end of string/g' > del_pos_d.txt
	cat del_pos_d.txt | sed -e 's/end of string/\n/g' > del_pos_e.txt
	cat del_pos_e.txt | sed -e 's/^ //g' > del_pos_f.txt
	cat del_pos_f.txt | awk '!NF{$0="0"}1' > del_pos_g.txt
	cat del_pos_g.txt | sed -e 's/ /\t/g' > del_pos.tsv
	paste read.tsv pos.tsv del_pos.tsv > $JOB_OUTPUT_DIR/$d/del_+_$d.tsv
	echo "del position extraction +"


	# Extract important lines (for - positions)
	grep "^GAA-" -B 2 $INPUT > seq.txt
	echo "Extracted lines -"

	# Extract start and stop lines
	grep -B 1 "^GAA-" seq.txt | awk -F '\n' 'ln ~ /^$/ { ln = "matched"; print $1 } $1 ~ /^--$/ { ln = "" }' > start_stop.txt
	echo "start stop extraction -"

	# Extract start and Stop Position
	cat start_stop.txt | sed -e 's/^.* .* .*bp //g' -e 's/ .*$//g' > start_stop_a.txt
	cat start_stop_a.txt | sed -e 's/-/\t/g' > pos.tsv
	echo "start stop position extraction -"

	# Extract read name
	cat start_stop.txt |awk '{print $1}' > read.txt
	cat read.txt | sed -e 's/ /\t/g' > read.tsv

	# Extract Insertions and Mismatches
	grep "^GAA-" seq.txt > in_mm.txt
	cat in_mm.txt | sed -e 's/^AGG-                                   .*bp score=.* //g' > in_mm_a.txt
	cat in_mm_a.txt | perl -nle '$pos=-1;while (($off=index($_,"-",$pos))>=0) {print $off;$pos=$off+1;}; print "end of string"' > in_pos.txt

	cat in_mm_a.txt | perl -nle '$pos=-1;while (($off=index($_,"a",$pos))>=0) {print $off;$pos=$off+1;}; print "end of string"' > mm_a_pos.txt
	cat in_mm_a.txt | perl -nle '$pos=-1;while (($off=index($_,"g",$pos))>=0) {print $off;$pos=$off+1;}; print "end of string"' > mm_g_pos.txt
	cat in_mm_a.txt | perl -nle '$pos=-1;while (($off=index($_,"t",$pos))>=0) {print $off;$pos=$off+1;}; print "end of string"' > mm_t_pos.txt
	cat in_mm_a.txt | perl -nle '$pos=-1;while (($off=index($_,"c",$pos))>=0) {print $off;$pos=$off+1;}; print "end of string"' > mm_c_pos.txt

	cat in_pos.txt | sed -e 's/$/a/g' > in_pos_a.txt
	tr --delete '\n' < in_pos_a.txt > in_pos_b.txt
	cat in_pos_b.txt | sed -e 's/a/ /g' > in_pos_c.txt
	cat in_pos_c.txt | sed -e 's/end of string end of string/end of string 0 end of string/g' > in_pos_d.txt
	cat in_pos_d.txt | sed -e 's/end of string/\n/g' > in_pos_e.txt
	cat in_pos_e.txt | sed -e 's/^ //g' > in_pos_f.txt
	cat in_pos_f.txt | awk '!NF{$0="0"}1' > in_pos_g.txt
	cat in_pos_g.txt | sed -e 's/ /\t/g' > in_pos.tsv
	paste read.tsv pos.tsv in_pos.tsv > $JOB_OUTPUT_DIR/$d/ins_-_$d.tsv
	echo "inn and mm position extraction -"
	
	cat mm_a_pos.txt | sed -e 's/$/a/g' > mm_a_pos_a.txt
	tr --delete '\n' < mm_a_pos_a.txt > mm_a_pos_b.txt
	cat mm_a_pos_b.txt | sed -e 's/a/ /g' >mm_a_pos_c.txt
	cat mm_a_pos_c.txt | sed -e 's/end of string end of string/end of string 0 end of string/g' > mm_a_pos_d.txt
	cat mm_a_pos_d.txt | sed -e 's/end of string/\n/g' > mm_a_pos_e.txt
	cat mm_a_pos_e.txt | sed -e 's/^ //g' > mm_a_pos_f.txt
	cat mm_a_pos_f.txt | awk '!NF{$0="0"}1' > mm_a_pos_g.txt
	cat mm_a_pos_g.txt | sed -e 's/ /\t/g' > mm_a_pos.tsv
	paste read.tsv pos.tsv mm_a_pos.tsv > $JOB_OUTPUT_DIR/$d/mm_a_-_$d.tsv
	echo "inn and mm position extraction +"
	
	cat mm_g_pos.txt | sed -e 's/$/a/g' > mm_g_pos_a.txt
	tr --delete '\n' < mm_g_pos_a.txt > mm_g_pos_b.txt
	cat mm_g_pos_b.txt | sed -e 's/a/ /g' >mm_g_pos_c.txt
	cat mm_g_pos_c.txt | sed -e 's/end of string end of string/end of string 0 end of string/g' > mm_g_pos_d.txt
	cat mm_g_pos_d.txt | sed -e 's/end of string/\n/g' > mm_g_pos_e.txt
	cat mm_g_pos_e.txt | sed -e 's/^ //g' > mm_g_pos_f.txt
	cat mm_g_pos_f.txt | awk '!NF{$0="0"}1' > mm_g_pos_g.txt
	cat mm_g_pos_g.txt | sed -e 's/ /\t/g' > mm_g_pos.tsv
	paste read.tsv pos.tsv mm_g_pos.tsv > $JOB_OUTPUT_DIR/$d/mm_g_-_$d.tsv
	echo "inn and mm position extraction +"
	
	cat mm_c_pos.txt | sed -e 's/$/a/g' > mm_c_pos_a.txt
	tr --delete '\n' < mm_c_pos_a.txt > mm_c_pos_b.txt
	cat mm_c_pos_b.txt | sed -e 's/a/ /g' >mm_c_pos_c.txt
	cat mm_c_pos_c.txt | sed -e 's/end of string end of string/end of string 0 end of string/g' > mm_c_pos_d.txt
	cat mm_c_pos_d.txt | sed -e 's/end of string/\n/g' > mm_c_pos_e.txt
	cat mm_c_pos_e.txt | sed -e 's/^ //g' > mm_c_pos_f.txt
	cat mm_c_pos_f.txt | awk '!NF{$0="0"}1' > mm_c_pos_g.txt
	cat mm_c_pos_g.txt | sed -e 's/ /\t/g' > mm_c_pos.tsv
	paste read.tsv pos.tsv mm_c_pos.tsv > $JOB_OUTPUT_DIR/$d/mm_c_-_$d.tsv
	echo "inn and mm position extraction +"
	
	cat mm_t_pos.txt | sed -e 's/$/a/g' > mm_t_pos_a.txt
	tr --delete '\n' < mm_t_pos_a.txt > mm_t_pos_b.txt
	cat mm_t_pos_b.txt | sed -e 's/a/ /g' >mm_t_pos_c.txt
	cat mm_t_pos_c.txt | sed -e 's/end of string end of string/end of string 0 end of string/g' > mm_t_pos_d.txt
	cat mm_t_pos_d.txt | sed -e 's/end of string/\n/g' > mm_t_pos_e.txt
	cat mm_t_pos_e.txt | sed -e 's/^ //g' > mm_t_pos_f.txt
	cat mm_t_pos_f.txt | awk '!NF{$0="0"}1' > mm_t_pos_g.txt
	cat mm_t_pos_g.txt | sed -e 's/ /\t/g' > mm_t_pos.tsv
	paste read.tsv pos.tsv mm_t_pos.tsv > $JOB_OUTPUT_DIR/$d/mm_t_-_$d.tsv
	echo "inn and mm position extraction +"

	# Extract Deletions 
	cat start_stop.txt | sed -e 's/^.* .* .*bp .*-.* //g' > del_a.txt
	cat del_a.txt | perl -nle '$pos=-1;while (($off=index($_,"-",$pos))>=0) {print $off;$pos=$off+1;}; print "end of string"' > del_pos.txt

	cat del_pos.txt | sed -e 's/$/a/g' > del_pos_a.txt
	tr --delete '\n' < del_pos_a.txt > del_pos_b.txt
	cat del_pos_b.txt | sed -e 's/a/ /g' > del_pos_c.txt
	cat del_pos_c.txt | sed -e 's/end of string end of string/end of string 0 end of string/g' > del_pos_d.txt
	cat del_pos_d.txt | sed -e 's/end of string/\n/g' > del_pos_e.txt
	cat del_pos_e.txt | sed -e 's/^ //g' > del_pos_f.txt
	cat del_pos_f.txt | awk '!NF{$0="0"}1' > del_pos_g.txt
	cat del_pos_g.txt | sed -e 's/ /\t/g' > del_pos.tsv
	paste read.tsv pos.tsv del_pos.tsv > $JOB_OUTPUT_DIR/$d/del_-_$d.tsv
	echo "del position extraction -"
done



