# Working with RNAseq data on the HPC

## Preface

Most of these scripts were adapted from ones I got from Jeff, or were built with the help of AI. They were created mostly to handle the RNAseq data from the 12 treatments and control groups in my nutrient deficiency project, hence the looping and nutrient names. Looking back, most of these have a lot of room for improvement in terms of efficiency and syntax (especially when it comes to looping - I would do arrays instead!). But, hopefully they are a good jumping off point.

## Paired-end read QC

For QC, I used fastQC. I would also recommend fastp (does trimming in addition to QC).
This is the HPC script I used to run fastQC:
```
#!/bin/bash
# Set the directory where your .fq.gz files are
input_dir="/grphome/fslg_atriplex/nobackup/archive/Ashley/nutrient_study/Ashley_nutrient/all_treatments"
# Make sure output folders exist
mkdir -p scripts
mkdir -p fastqc_output
# Find all .fq.gz files
find -L "$input_dir" -type f -name "*.fq.gz" | while read file
do
    # Get just the filename
    filename=(basename"file")
    name="${filename%.fq.gz}"
    # Create a subfolder for each sample's script
    mkdir -p scripts/${name}
    # Create a SLURM batch script
    cat > scripts/name/{name}.sh <<EOF
#!/bin/bash
#SBATCH --time=72:00:00   # walltime
#SBATCH --ntasks=24   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4096M   # memory per CPU core
#SBATCH --job-name=fastqc_${name}
#SBATCH --output=scripts/name/{name}.out
#SBATCH --error=scripts/name/{name}.err
module purge
module load conda-pws
module load conda/fastqc
fastqc -t 24 ${file} -o fastqc_output  
EOF
sbatch ./scripts/name/{name}.sh
done
```
## Read trimming

Trim any lingering adapter sequences using Trimmomatic. This script is set for paired-end reads:
```
#!/bin/sh
#this is a shell script that you launch using:$ 'sh trimmomatic_paired-end.sh' 
#it will launch separate jobs for all pairs of  illumina files it finds in the file where it was launched - it does this by looking for a forward file of the pair, labeled _1.fq.gz (this could be changed) -the reverse file is _2.fq.gz.
#this is for PE (paried-end illumina) files.
#it is checking for all of the adapters in the following kits: NexteraPE-PE.fa  TruSeq2-PE.fa  TruSeq2-SE.fa  TruSeq3-PE-2.fa  TruSeq3-PE.fa  TruSeq3-SE.fa
#make sure you make a trim_scripts and a Trimmed_Reads folders in the folder you launch the script from.
# Set up directories
mkdir -p trim_scripts
mkdir -p Trimmed_Reads
# Loop over all forward read files matching *_1.fq.gz
for forward_file in *_1.fq.gz
do
    # Get the sample name by removing _1.fq.gz
    name="${forward_file%_1.fq.gz}"
    # Derive the matching reverse read filename
    reverse_file="${name}_2.fq.gz"
cat > ./trim_scripts/${name}.sh <<EOF
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4096M
#SBATCH -J "trimmomatic"
#SBATCH --mail-user=ashmarch@byu.edu
#SBATCH --mail-type=FAIL
#SBATCH --qos=pws
module purge
module load conda-pws
module load trimmomatic
trimmomatic PE -threads 4 -summary ${name}_stats_trim.txt \
${forward_file} ${reverse_file} \
-baseout ./Trimmed_Reads/${name}.fq.gz \
ILLUMINACLIP:/nobackup/archive/grp/fslg_pws_module/software/.conda/envs/trimmomatic/adapters/all_adaptor.fasta:2:40:15 \
LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:75
EOF
sbatch ./trim_scripts/${name}.sh
done
```

## Read mapping
To map reads back to your reference, use HISAT2. These settings allow for read multimapping for proper feature counting in a polyploid. These scripts I ran separately in multiple steps, but could easily be combined.

### Step 1. Build an index
```
#!/bin/bash
#SBATCH --time=72:00:00   # walltime
#SBATCH --ntasks=32   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=16G   # memory per CPU core
#SBATCH -J "hisat2-build"   # job name
#SBATCH --mail-user=ashmarch@byu.edu   # email address
#SBATCH --mail-type=FAIL
#SBATCH --qos=pws

module purge
module load conda-pws
module load hisat2_v2.2.1
reference=Cquinoa_QQ74_v2_final.fasta #Don't give the full path.  Run this where the reference is.
hisat2-build -p 32 ${reference} ${reference}
```

### Step 2. Map reads
```
#!/bin/sh

reference=/grphome/fslg_atriplex/nobackup/archive/Ashley/nutrient_study/hisat2/Cquinoa_QQ74_v2_final.fasta   #full path to reference and index.
reads=/grphome/fslg_atriplex/nobackup/archive/Ashley/nutrient_study/Ashley_nutrient/all_treatments/Trimmed_Reads/Mo  #full path to reads directory
output=/grphome/fslg_atriplex/nobackup/archive/Ashley/nutrient_study/hisat2/Mo  #full path to output directory

mkdir -p ${output}/hisat_scripts

for forward_file in ${reads}/*_1P.fq.gz; do
    name=(basename "forward_file" _1P.fq.gz)
    cat > output/hisat_scripts/{name}.sh <<EOF
#!/bin/bash

#SBATCH --time=08:00:00   # walltime
#SBATCH --ntasks=8   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4G   # memory per CKU core
#SBATCH -J "hisat2"   # job name
#SBATCH --mail-user=ashmarch@byu.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --qos=pws

module load conda-pws
module load hisat2_v2.2.1

hisat2 -p 8 -x ${reference} -1 reads/{name}_1P.fq.gz -2 reads/{name}_2P.fq.gz -k 10 -S output/{name}.sam

EOF

sbatch output/hisat_scripts/{name}.sh

done
```
### Step 3. Convert sam to bam and sort
```
#!/bin/sh

output=/grphome/fslg_atriplex/nobackup/archive/Ashley/nutrient_study/hisat2/Mo
mkdir -p ${output}/sam_scripts

for sam in ${output}/*.sam
do
    name=(basename"sam" .sam)
    cat > output/sam_scripts/{name}.sh <<EOF
#!/bin/bash

#SBATCH --time=72:00:00   # walltime
#SBATCH --ntasks=32   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=16G   # memory per CPU core
#SBATCH -J "sam2bam"   # job name
#SBATCH --mail-user=ashmarch@byu.edu   # email address
#SBATCH --mail-type=FAIL
#SBATCH --qos=pws

module purge
module load conda-pws
module load samtools/1.19.2

cd ${output}
samtools view -Sb --threads 32 ${name}.sam > ${name}.bam

samtools sort --threads 32 -o ${name}.sorted.bam ${name}.bam
samtools index -@ 32 ${name}.sorted.bam
EOF

sbatch output/sam_scripts/{name}.sh

done
```

## Feature counting
Because multimapping was allowed during read mapping, with these parameters features are able to be counted fractionally with featureCounts. This are the scripts I made to run this, and once again, it could all be combined if desired:

### Step 1: Convert gff to gtf (I used longest isoform gff)
```
#!/bin/bash

#SBATCH --time=01:00:00   # walltime
#SBATCH --ntasks=1        # number of cores
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH -J "gff2gtf"      # job name
#SBATCH --mail-user=ashmarch@byu.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --qos=pws

module purge
module load conda-pws
module load gffread_v0.12.7

# Replace these paths with yours
gff_input=/grphome/fslg_atriplex/nobackup/archive/Ashley/nutrient_study/featureCounts/quinoa_braker_functional_longest_isoform.gff3
gtf_output=/grphome/fslg_atriplex/nobackup/archive/Ashley/nutrient_study/featureCounts/quinoa_braker_functional_longest_isoform.gtf

gffread $gff_input -T -o $gtf_output
```
### Step 2: Run featureCounts
```
#!/bin/sh

output=/grphome/fslg_atriplex/nobackup/archive/Ashley/nutrient_study/hisat2/Mo #change as needed
gtf=quinoa_braker_functional_longest_isoform.gtf #change as needed

mkdir -p featurecounts_scripts

for bam in ${output}/*.coord.sorted.bam
do
    name=(basename"bam" .coord.sorted.bam)
    script=featurecounts_scripts/${name}_featurecounts.sh
    cat > $script <<EOF
#!/bin/bash

#SBATCH --time=08:00:00
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8G
#SBATCH -J "featureCounts"
#SBATCH --mail-user=ashmarch@byu.edu
#SBATCH --mail-type=FAIL
#SBATCH --qos=pws

module purge
module load conda-pws
module load subread_v2.0.8

featureCounts -T 8 -p --countReadPairs -M --fraction -t exon -g gene_id -a $gtf -o ${name}_counts.txt $bam
EOF
        sbatch $script
done
```

## Preparing files for DEG analysis in R
To make the necessary files for DESeq2 in R, I had to combine the counts.txt files that were created by featureCounts into one big table. Also, I needed another table with information about the samples. Here is how to prepare those files:

### Step 1: Make a merged count matrix
```
#!/bin/bash
output="L_combined_counts_matrix.tsv"
count_files=(*_L_*_counts.txt)
# Sort and clean all files
for f in "${count_files[@]}"; do
    grep -v '^#' "$f" | awk 'NR>1 {print $1"\t"$7}' | sort > "${f}.clean"
done
# Initialize with first file
cp "${count_files[0]}.clean" merged.tsv
# Merge all remaining
for f in "${count_files[@]:1}"; do
    join -t $'\t' merged.tsv "${f}.clean" > tmp && mv tmp merged.tsv
done
# Build header
header="Geneid"
for f in "${count_files[@]}"; do
    sample=(grep-v'#'"f" | head -1 | awk '{print $7}')
    sample=(basename"sample".sorted.bam)
    header="$header\t$sample"
done
# Add header
echo -e "header"|cat-merged.tsv>"output"
# Cleanup
rm *.clean merged.tsv
echo "✅ Safely merged matrix written to $output"
```

### Step 2: Create a table with this format (change accordingly)
```
sample  condition  replicate
B_L_1  B  1
B_L_2  B  2
...
```


