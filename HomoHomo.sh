# HOMOlogy HOMOpolymer repair
# started 05.01.2015. AI Manolov
# requirements: bowtie2, samtools, bedtools, perl, BioPerl 

source config.txt

# A r g u m e n t s
args="$#"
echo "args = $args"

if [ $args -lt 2 ]
then
	echo "HomoHomo.sh <fasta of assembly> <fastq of reads> [bam] [vcf]"
	exit -1
fi

fastaFile=$1
fastqFile=$2

fastaName=$(basename "$1")
ext="${fastaName##*.}"
fastaName="${fastaName%.*}"

fastqName=$(basename "$2")
ext="${fastqName##*.}"
fastqName="${fastqName%.*}"

projName=${fastaName}_${fastqName}

bamFile=""
vcfFile=""

if [ $args -ge 3 ]
then
	bamFile=$3	
	if [ $args -eq 4 ]
	then 
		vcfFile=$4
	fi
fi

echo "fasta: $fastaName fastq: $fastqName  bam: $bamFile  vcf: $vcfFile"
echo "PROJECT NAME: $projName"

# F o l d e r s
mkdir -p tmp
mkdir -p tmp/reads
mkdir -p tmp/assembly
mkdir -p tmp/bam
mkdir -p tmp/vcf
mkdir -p tmp/pieces
mkdir -p corrected

# M a i n

# bam

if [[ -e "tmp/bam/${projName}.bam" ]]  &&  [[ "$reuse" == 1 ]]
then
	bamFile="tmp/bam/${projName}.bam"	
fi

if [ "$bamFile" = "" ]
then
	echo "no BAM file. Mapping..."
        mkdir -p tmp/bam/bwt
	bamFile="tmp/bam/$projName.bam"
	bowtie2-build $fastaFile tmp/bam/bwt/$fastaName
        bowtie2 tmp/bam/bwt/$fastaName $fastqFile -p $threads --no-unal -S tmp/bam/$projName.sam
	samtools view -Sb tmp/bam/$projName.sam  > tmp/bam/${projName}_.bam 
        samtools sort tmp/bam/${projName}_.bam tmp/bam/$projName
	rm tmp/bam/${projName}_.bam
        samtools index $bamFile 
	[ ! -f $bamFile ] && { echo "mapping failed..."; exit 1; }
fi
echo "BAM file: $bamFile"

# vcf
if [[ -e "tmp/vcf/$projName.vcf" ]] && [[ "$reuse" == 1 ]]
then
	vcfFile="tmp/vcf/$projName.vcf"	
fi

if [ "$vcfFile" == "" ]
then
	echo "no VCF file. Calculating..."
	vcfFile="tmp/vcf/$projName.vcf"
	samtools mpileup -uf $fastaFile $bamFile | bcftools view -gvc - > $vcfFile
	[ ! -f $vcfFile ] && { echo "vcf calculation failed..."; exit 2; }
fi
echo "VCF file: $vcfFile"

# pieces
perl makeIndelFasta.pl $vcfFile $fastaFile $projName $indent $lowFreq 

fasta2blast="tmp/pieces/indel_pieces_$projName.fasta"
fastaFromBed -name -bed "tmp/pieces/indel_pieces_$projName.bed" \
		-fi $fastaFile -fo $fasta2blast;

[ ! -f $fasta2blast ] && { echo "no pieces for blast..."; exit 3; }

# blast
echo "BLASTING..."
blastOut="tmp/pieces/indel_pieces_$projName.xml"
blastn -db "nt" -query $fasta2blast -num_threads $threads -outfmt 5 \
		-max_target_seqs 1 -task 'blastn' -evalue $evalue \
		> $blastOut

[ ! -f $blastOut ] && { echo "no blast results..."; exit 4; }

# repair
perl repair.pl $blastOut tmp/pieces/indel_pieces_$projName.vcf $fastaFile ${projName} $indent 

