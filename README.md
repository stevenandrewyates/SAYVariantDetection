# SAYVariantDetection

A GitHub repository for getting information from short-read genome alignments.

# Description

This guide will help you detect variants using SAMtools on the  Euler computing system at ETH Zurich (http://scicomp.ethz.ch/wiki/Euler).

 It assumes that you have two things; a directory full of sorted *BAM* files and a genome reference
 
# Simple variant detetction

I like SAMtools, so that's what we're going to use. To detect variants you need some sorted *BAM* files, in a directory. To save typing all of these we will make a variable `$BAMSort` that contains them all by looking (`find`) for them in a directory: in this case called *BAMsorted*. This can be made using the GitHub repository [here](https://github.com/stevenandrewyates/SAYReadMappingDNA). Below is the line to make the `$BAMSort` variable.

```
BAMSort=$(find BAMsorted -name '*bamsort')
```

Next we need to load the *SAMtools* module.

```
module load gdc
module load samtools/1.7
```

next we can begin variant calling using

```
samtools mpileup -ug --output-tags DP -f GENOME/genome.fasta $BAMSort | bcftools call -mv | vcfutils.pl varFilter -d 60 > output.vcf
```

In the above we make a `samtools mpileup` file and specify to: report the depth (`--output-tags DP`) and the genome sequence in fasta format (`-f GENOME/genome.fasta`). This is then piped (`|`) to `bcftools call` Finally the output is filtered (`vcfutils.pl varFilter`) for minimum depth of of 60 in this case. Finally we write (`>`) the data to an `output.vcf`.

All of this can be done with the following:

```
module load git
git clone https://github.com/stevenandrewyates/SAYVariantDetection
sh SAYVariantDetection/01_simpleVCF.sh -b BAMsorted -g GENOME/genome.fasta -o output.vcf -d 60
```

Where the following arguments are passed:

- -b a directory with the sorted BAM files
- -g the FASTA file with your genome (nicely formatted)
- -o an output VCF file
- -d the minimum depth  

# Example 
For heuristic reasons the full code from zero to *vcf* can be found below:

```
cd $SCRATCH
mkdir CASS
cd CASS

module load git
git clone https://github.com/stevenandrewyates/SAYEulerDataManagement
sh SAYEulerDataManagement/01_DownloadSRAtoolkit.sh
sh SAYEulerDataManagement/02_DownloadCassavaSix.sh -f CASS -n 100000


git clone https://github.com/stevenandrewyates/SAYReadMappingDNA
sh SAYReadMappingDNA/01_DownloadGenome.sh -f ftp://ftp.ensemblgenomes.org/pub/plants/release-50/fasta/manihot_esculenta/dna/Manihot_esculenta.Manihot_esculenta_v6.dna.toplevel.fa.gz
bash SAYReadMappingDNA/03_ReadMapBsub.sh -N 3 -R 2 -M 2 -F FASTQ -G GENOME
bash SAYReadMappingDNA/04_ReadMapBsubBAM.sh -N 1 -R 1 -M 1 -S SAM -G genome.fasta
bash SAYReadMappingDNA/05_ReadMapBsubBAMsorted.sh -N 1 -R 1 -M 1 -S BAM -G genome.fasta


git clone https://github.com/stevenandrewyates/SAYVariantDetection
sh SAYVariantDetection/01_simpleVCF.sh -b BAMsorted -g GENOME/genome.fasta -o output.vcf -d 60
```

# Future

This will only get you so far, because it may take longer than Euler will give you or more RAM than needed. Soon I will add a version to run this using `bsub` where you can increase the time and RAM.


```
bsub -W 240 -n 4 -R "rusage[mem=1024]" samtools mpileup -ug --output-tags DP -f GENOME/genome.fasta $BAMSort | bcftools call -mv | vcfutils.pl varFilter -d 60 > output.vcf
```

Although the line above will specify four hours run time (`-W 240`), with four cores (`-n 4`), each with 1 GB of RAM (`"rusage[mem=1024]"` so 4GB RAM total).
