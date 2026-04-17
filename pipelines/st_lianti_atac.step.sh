#!/usr/bin/env bash

# GitHub-ready entrypoint copied from the development workspace.
# Machine-specific tools and references are resolved from config/pipeline.env.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
export SPATIAL_EPIGENOME_REPO_ROOT="${repo_root}"
# shellcheck source=../lib/common.sh
source "${repo_root}/lib/common.sh"
load_pipeline_config "${repo_root}"
SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-$(threads_or_default 4)}"


## software
cutadapt="${CUTADAPT:-cutadapt}"
bowtie2="${BOWTIE2:-bowtie2}"

thread=12
wd=`pwd`
ddate=$(date -R | awk '{print $2$3}')
BST_src="${BST_SRC_S1000:-/path/to/BSTMatrix_v2.3.r}"
umi_tools="${UMI_TOOLS:-umi_tools}"
spatial_lianti="${SPATIAL_LIANTI:-spatial_lianti}"
Rscript="${RSCRIPT:-Rscript}"
plot_abortive="${src}/plot_abortive_distribution.R"

if [[ ! -d ./logs ]];then   
    mkdir ./logs
fi 


function barcode_calling(){
    read1=$1
    read2=$2
    outdir=$3

    job="barcode_calling"
    if [[ ! -d $outdir ]];then   
        mkdir -p $outdir
    fi 

    sample=$(basename -s "_R1.fastq.gz" $read1)
    job=${job}_${sample}
    sbatch <<RUN
#!/usr/bin/env bash
#SBATCH -J $job
#SBATCH -o logs/${job}.${ddate}.log
#SBATCH -e logs/${job}.${ddate}.log
#SBATCH --mem=120g
#SBATCH --cpus-per-task=$thread
#SBATCH --partition=compute_fat
#SBATCH -t 7-00:00:00
    export PERL5LIB="/usr/lib64/perl5/vendor_perl"

    $BST_src/fastq2BcUmi -i1 $read1 \
        -i2 $read2 \
        -d $outdir \
        -p $sample \
        -n $thread

RUN
}


function cutadapt(){
    read1=$1
    read2=$2
    outdir=$3

    job="cutadapt"
    if [[ ! -d $outdir ]];then   
        mkdir -p $outdir
    fi 

    sample=$(basename -s "_R1.fastq.gz" $read1 | sed s'/SQ.*-BMK//')
    job=${job}_${sample}
    sbatch <<RUN
#!/usr/bin/env bash
#SBATCH -J $job
#SBATCH -o logs/${job}.${ddate}.log
#SBATCH -e logs/${job}.${ddate}.log
#SBATCH --mem=120g
#SBATCH --cpus-per-task=$thread
#SBATCH --partition=compute_fat
#SBATCH -t 7-00:00:00
    $cutadapt -j $thread \
        -G "^NAGATGTGTATAAGAGACAG" \
        -G "^NNAGATGTGTATAAGAGACAG" \
        -G "^NNNAGATGTGTATAAGAGACAG" \
        -G "^NNNNAGATGTGTATAAGAGACAG" \
        --no-indels \
        -e 2 \
        -o $outdir/${sample}_R1.fastq.gz \
        -p $outdir/${sample}_R2.fastq.gz \
        --untrimmed-output $outdir/${sample}_woME_R1.fastq.gz \
        --untrimmed-paired-output $outdir/${sample}_woME_R2.fastq.gz \
        $read1 $read2
RUN
    done 
}

function stat_abortive(){
    read1=$1
    read2=$2
    outdir=$3

    job="stat_abortive"
    
    if [[ ! -d $outdir ]];then   
        mkdir -p $outdir
    fi 

    sample=$(basename -s "_R1.fastq.gz" $read1 | sed s'/SQ.*-BMK//')
    job=${job}_${sample}
    sbatch <<RUN
#!/usr/bin/env bash
#SBATCH -J $job
#SBATCH -o logs/${job}.${ddate}.log
#SBATCH -e logs/${job}.${ddate}.log
#SBATCH --mem=120g
#SBATCH --cpus-per-task=$thread
#SBATCH --partition=compute_fat
#SBATCH -t 7-00:00:00
    $cutadapt -j $thread \
        -A "A{9};min_overlap=9" \
        --no-indels \
        -e 2 \
        -o $outdir/${sample}_polyA_trimmed.R1.fastq.gz \
        -p $outdir/${sample}_polyA_trimmed.R2.fastq.gz \
        $read1 $read2
    
    seqkit fx2tab $outdir/${sample}_polyA_trimmed.R2.fastq.gz -l -i -Q -n > $outdir/len.stat.txt

    $Rscript $plot_abortive_distribution $outdir/len.stat.txt  $outdir/$sample


RUN
    done 

}

function barcode2fq(){
    fastq=$1
    bc_read_map=$2
    umi_cor=$3
    outdir=$4

    if [[ ! -d $outdir ]];then   
        mkdir -p $outdir
    fi 

    sampleid=$(basename -s "_R2.fastq.gz" $fastq)
    job="bc2fq"
    sbatch <<RUN
#!/usr/bin/env bash
#SBATCH -J $job
#SBATCH -o logs/${job}.${ddate}.log
#SBATCH -e logs/${job}.${ddate}.log
#SBATCH --mem=120g
#SBATCH --cpus-per-task=$thread
#SBATCH --partition=compute_fat
#SBATCH -t 7-00:00:00

    
    $spatial_lianti tag2fq \
        --fastq $fastq \
        --bc_read_map $bc_read_map \
        --umi_cor $umi_cor \
        -o $outdir \
        --n_lines 1663753936 \
        --sampleName $sampleid \
        --threads $thread \
        --out_unmapped \
        --verbose   
RUN
}