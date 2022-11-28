cd "$1"
ref="$2"
threads="$3"
threads="${threads:-4}"
for R1 in $(ls *R1_001.fastq.gz)
do
	smpl="${R1/_R1_001.fastq.gz/}"
	# R1="${smpl}_R1_001.fastq.gz"
	R2="${smpl}_R2_001.fastq.gz"

	bwa mem -t $threads "$ref" "$R1" "$R2" > "${smpl}.sam"

	samtools view -bS "${smpl}.sam" | samtools sort -@ $threads -o "${smpl}.bam"
	samtools index "${smpl}.bam"

	samtools mpileup -q 0 -Q 0 -B -d 10000000 -A -f "$ref" \
	     "${smpl}.bam" | java -jar /Users/semiquant/bioinfomatics/VarScan.jar mpileup2cns --variants 1 --output-vcf 1 --min-coverage 1 --min-avg-qual 0 \
	    --min-var-freq 0 --strand-filter 0 --p-value 1 --min-reads2 1 > "${smpl}.vcf"

	for ref_name in $(grep "^>" "$ref" | cut -c2-)
	do
		process.py 20 "${smpl}.bam" "$ref_name" "$ref"
	done

done