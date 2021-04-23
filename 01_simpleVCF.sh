##################################################
######            Details                   ######
##################################################

# author:   Steven Yates
# contact:  steven.yates@usys.ethz.ch
# Year:     2021
# Citation: TBA

##################################################
######            Description               ######
##################################################

# A script to download a genome and install a 
# bowtie2 reference using the Euler computing system.


##################################################
######              Usage                   ######
##################################################

# This script will make a file for sorted BAM files in a given folder. 
# It requires the following inputs:

# -b a directory with the sorted BAM files
# -g the FASTA file with your genome (nicely formatted)
# -o an output VCF file
# -d the minimum depth  

# sh SAYVariantDetection/01_simpleVCF.sh -b BAMsorted -g GENOME/genome.fasta -o output.vcf -d 60

##################################################
######              Script                  ######
##################################################


while getopts b:g:o:d: flag
do
    case "${flag}" in
        b) BAM=${OPTARG};;
 	g) GEN=${OPTARG};;
 	o) OUT=${OPTARG};;
 	d) DEP=${OPTARG};;

    esac
done

echo "looking for BAM files in directory:";
echo $BAM;

BAMSort=$(find $BAM -name '*bamsort')

echo "found the following bam files:";
echo $BAMsort;

echo "genome index is:";
echo $GEN; 

echo "Loading modules:";
module load gdc
module load samtools/1.7

echo "aligning to genome:"
samtools mpileup -ug --output-tags DP -f $GEN $BAMSort | bcftools call -mv | vcfutils.pl varFilter -d $DEP > $OUT

echo "Done :)";
