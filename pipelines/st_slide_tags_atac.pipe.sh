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



set -e
set -o pipefail

helpdoc() {
    cat <<EOF
Description:
    [main::] Spatial ATAC-seq based on IVT data of slide-tags processing pipeline created by Maoxu Wang, Xie Lab, Peking University.
    Last update: 2025.04.12
Usage:
    $run [options]
        -h print user manual information
Options:
    --atac_read1        -- atac read1 fastq.gz file
    --atac_read2        -- atac read2 fastq.gz file
    --spatial_r1        -- spatial read1 fastq.gz file
    --spatial_r2        -- spatial read2 fastq.gz file
    --hdmi_read         -- hdmi read fastq.gz file
    --sampleid          -- sampleid name specified 
    --species           -- e.g. hg38 or mm10
    --n_line            -- number of lines to process (default: all lines)
    --outdir            -- output directory
    --option            -- steps to run
    --debug             -- debug mode
    --help              -- print manual
EOF
}
ARGS=$(getopt -o o:i:h -l atac_read1:,atac_read2:,spatial_r1:,spatial_r2:,chip_id:,DAPI:,HE:,hdmi_read:,sampleid:,species:,n_line:,outdir:,min_mismatch:,option:,debug:,help -- "$@")
if [ $? != 0 ]; then
    echo "terminating..." >&2
    exit 1
fi

eval set -- "$ARGS"

while true; do
    case "$1" in
    --atac_read1)
        atac_read1=$2
        shift 2
        ;;
    --atac_read2)
        atac_read2=$2
        shift 2
        ;;
    --spatial_r1)
        spatial_r1=$2
        shift 2
        ;;
    --spatial_r2)
        spatial_r2=$2
        shift 2
        ;;
    --hdmi_read)
        hdmi_read=$2
        shift 2
        ;;
    --chip_id)
        chip_id=$2
        shift 2
        ;;
    --DAPI)
        DAPI=$2
        shift 2
        ;;
    --HE)
        HE=$2
        shift 2
        ;;
    --sampleid)
        sampleid=$2
        shift 2
        ;; 
    --species)
        species=$2
        shift 2
        ;;   
    --n_line)
        n_line=$2
        shift 2
        ;;
    -o | --outdir)
        outdir=$2
        shift 2
        ;;
    --option)
        option=$2
        shift 2
        ;;
    --debug)
        debug=$2
        shift 2
        ;;
    --min_mismatch)
        min_mismatch=$2
        shift 2
        ;;
    -h | --help)
        helpdoc
        exit 1
        ;;
    --)
        shift
        echo "Parameters recevied."
        break
        ;;
    *)
        echo "unknown parameter:" $1
        helpdoc
        exit 1
        ;;
    esac
done

### cell barcode on atac read1
if [[ ! -n $outdir ]]; then outdir="./result"; fi
if [[ ! -n $n_line ]]; then n_line="-1"; fi
if [[ ! -n $option ]]; then option="pipe"; fi

if [[ ! -n $debug ]]; then debug=False; fi
if [[ ! -n $min_mismatch ]]; then min_mismatch=2; fi

if [[ ! -n $sampleid ]]; then
    echo "Please specify sampleid name"
    exit
fi

if [[ ! -n $species ]]; then
    echo "Please specify species"
    exit
fi

if [[ $SLURM_CPUS_PER_TASK > 12 ]];then 
    sort_thread=12
else
    sort_thread=$SLURM_CPUS_PER_TASK
fi 

trim_galore="${TRIM_GALORE:-trim_galore}"
cutadapt="${CUTADAPT:-cutadapt}"
bedtools="${BEDTOOLS:-bedtools}"
bowtie2="${BOWTIE2:-bowtie2}"
samtools="${SAMTOOLS:-samtools}"
picard="${PICARD:-picard}"
Rscript="${RSCRIPT:-Rscript}"
macs2="${MACS2:-macs2}"
bamCoverage="${BAMCOVERAGE:-bamCoverage}"
umi_tools="${UMI_TOOLS:-umi_tools}"
spatial_IVT="${SPATIAL_IVT:-spatial_IVT}"
seqkit="${SEQKIT:-seqkit}"
python="${PYTHON:-python}"
python_cv2="${PYTHON_CV2:-python}"
bigWigToBedGraph="${BIGWIG_TO_BEDGRAPH:-bigWigToBedGraph}"
bedGraphToBigWig="${BEDGRAPH_TO_BIGWIG:-bedGraphToBigWig}"
src="${SPATIAL_EPIGENOME_SCRIPT_DIR:-${repo_root}/scripts}"

slide_tag_src="${src}/slide_tag_atac"
spatial_src="${src}"

createSeuratObj="$slide_tag_src/createSeuratObj.R"

hg38_idx="${HG38_BOWTIE2_INDEX:-/path/to/bowtie2/hg38}"
mm10_idx="${MM10_BOWTIE2_INDEX:-/path/to/bowtie2/mm10}"
chrom_size_var="${species^^}_CHROM_SIZES"
chromSize="${!chrom_size_var}"
hg38_autosome="${HG38_AUTOSOME_BED:-/path/to/hg38.autosome.txt}"
mm10_autosome="${MM10_AUTOSOME_BED:-/path/to/mm10.autosome.txt}"
rmsk_hg38="${HG38_RMSK_BED:-/path/to/hg38.rmsk.bed}"
rmsk_mm10="${MM10_RMSK_BED:-/path/to/mm10.rmsk.bed}"
star_hg38="${HG38_STAR_INDEX:-/path/to/star/hg38}"
star_mm10="${MM10_STAR_INDEX:-/path/to/star/mm10}"
signac_chromVar="$src/DBiT/signac_chromVar.R"
create_gene_activity="${src}/create_gene_activity.py"
gtf_hg38="${HG38_GTF:-/path/to/hg38/genes.gtf}"
gtf_mm10="${MM10_GTF:-/path/to/mm10/genes.gtf}"

STAR="${STAR:-STAR}"

declare -A gtf_idxes
gtf_idxes["hg38"]=$gtf_hg38
gtf_idxes["mm10"]=$gtf_mm10
gtf=${gtf_idxes[$species]}

declare -A star_idxes
star_idxes["hg38"]=$star_hg38
star_idxes["mm10"]=$star_mm10
star_idx_path=${star_idxes[$species]}

read_len=20

declare -A Ref_fa
Ref_fa["hg38"]=$hg38_idx
Ref_fa["mm10"]=$mm10_idx
bowti2_idx=${Ref_fa[$species]}

declare -A effective_size
effective_size["hg38"]=2913022398
effective_size["mm10"]=2652783500
effective_genome_size=${effective_size[$species]}


if [[ $SLURM_CPUS_PER_TASK -ge 8 ]]; then
    thread_trim=8
else
    thread_trim=$SLURM_CPUS_PER_TASK
fi 

if [[ $SLURM_CPUS_PER_TASK -ge 12 ]]; then
    thread_sort=12
else
    thread_sort=$SLURM_CPUS_PER_TASK
fi 

function run_seekSpace(){
    seekspacetools run \
        --fq1 $atac_read1 \
        --fq2 $atac_read2 \
        --spatialfq1 $spatial_r1 \
        --spatialfq2 $spatial_r2 \
        --hdmifq $hdmi_read \
        --samplename $sampleid \
        --outdir $outdir/$sampleid \
        --star_path $STAR \
        --genomeDir $star_idx_path \
        --gtf $gtf \
        --chemistry DDVS \
        --core $SLURM_CPUS_PER_TASK \
        --include-introns \
        --forceCell 80000 \
        --min_umi 10 \
        --chip_id $chip_id \
        --DAPI $DAPI \
        --HE $HE
}


function stat_abortive(){
    
    if [[ ! -d $outdir/stat_abortive ]];then 
        mkdir -p $outdir/stat_abortive
    fi 
    thread_trim=8


    echo "Downsampling to 100K reads..."

    $cutadapt -j $thread_trim \
        -a "A{7};max_error_rate=0.2;min_overlap=3" \
        -a "AGATCGGAAGAGC" \
        -a "G{9};max_error_rate=0.2;min_overlap=3" \
        --trim-n \
        --no-indels \
        --json $outdir/stat_abortive/$sampleid.trimming_report.txt \
        -e 2 \
        -o $outdir/stat_abortive/${sampleid}_polyA_trimmed.R1.fastq.gz \
        $outdir/stat_abortive/$sampleid.100K.R1.fastq.gz

    $seqkit fx2tab $outdir/stat_abortive/${sampleid}_polyA_trimmed.R1.fastq.gz \
        -l -i -Q -n > $outdir/stat_abortive/$sampleid.len.stat.txt
    awk '{print $2}' $outdir/stat_abortive/$sampleid.len.stat.txt > $outdir/stat_abortive/$sampleid.len.stat.tmp
    mv $outdir/stat_abortive/$sampleid.len.stat.tmp $outdir/stat_abortive/$sampleid.len.stat.txt

    $Rscript $spatial_src/plot_abortive_distribution.R $outdir/stat_abortive/$sampleid.len.stat.txt $outdir/stat_abortive/$sampleid

}

function trim_fq(){
    if [[ $SLURM_CPUS_PER_TASK -ge 12 ]]; then
        thread_trim=12
    else
        thread_trim=$SLURM_CPUS_PER_TASK
    fi 

    if [[ ! -d $outdir/trim ]];then 
        mkdir -p $outdir/trim
    fi 

    ## specify atac_reads by seekSpace tool
    ### Cellbarcode_UMI are added to read name
    atac_read1="${outdir}/${sampleid}/Analysis/scRNA-seq_Analysis/step1/${sampleid}_1.fq.gz"
    atac_read2="${outdir}/${sampleid}/Analysis/scRNA-seq_Analysis/step1/${sampleid}_2.fq.gz"

    if [[ ! -s $atac_read1 ]]; then
        echo $atac_read1 "wrong Read1 specified"
        exit
    fi

    if [[ ! -s $atac_read2 ]]; then
        echo $atac_read2 "wrong Read2 specified"
        exit
    fi

    ## MEthod2: using first 50bp 
    # #input_fq="$outdir/barcode/${sampleid}.R2.barcoded.fastq.gz"

    zcat $atac_read2 | seqkit subseq -r 1:50 | pigz > $outdir/trim/${sampleid}.R2.barcoded_first50.fastq.gz
    input_fq="$outdir/trim/${sampleid}.R2.barcoded_first50.fastq.gz"

    if [[ ! -s $outdir/trim/${sampleid}_clean.R2.fastq.gz ]];then 
        $cutadapt -j $thread_trim \
            -q 20 \
            -a "A{9};max_error_rate=0.1;min_overlap=7" \
            -a "AGATCGGAAGAGC" \
            -a "G{9};max_error_rate=0.1;min_overlap=7" \
            --trim-n \
            --no-indels \
            --minimum-length $read_len \
            --json $outdir/trim/$sampleid.trimming_report.txt \
            -e 2 \
            -o $outdir/trim/${sampleid}_clean_barcoded.R2.fastq.gz \
            $input_fq
    else
        echo "sequence is trimmed already."
    fi 
}

function align(){

    if [[ ! -d $outdir/align ]];then 
        mkdir -p $outdir/align
    fi 
    input_fq="$outdir/trim/${sampleid}_clean_barcoded.R2.fastq.gz"

    # outdir/trim/${sampleid}_clean_barcoded.R2.fastq.gz
    if [[ ! -s $outdir/align/${sampleid}.sorted.bam ]];then 
        $bowtie2 -x $bowti2_idx \
            -U $input_fq \
            --no-unal \
            -p $SLURM_CPUS_PER_TASK \
            | awk -v OFS='\t' '$0~/^@/ { print; next } { split($1,bc,"_"); total=length(bc); print $0 "\tCB:Z:" bc[1] "\tUB:Z:" bc[2]}' \
            | $samtools sort -@ $thread_sort -m "12G" \
            > $outdir/align/${sampleid}.sorted.bam
    else
        echo "Align is finished."
    fi 

    $samtools index $outdir/align/${sampleid}.sorted.bam
    $samtools stats $outdir/align/${sampleid}.sorted.bam > $outdir/align/${sampleid}.bam.stat

}


function snapATAC2(){
    if [[ ! -d $outdir/fragments ]];then 
        mkdir -p $outdir/fragments
    fi 

    if [[ ! -d $outdir/spatial/${sampleid} ]];then 
        mkdir -p $outdir/spatial/${sampleid}
    fi 

    ## make tissue barcode 
    awk '{FS=",";print $1}' \
        ${outdir}/${sampleid}/Analysis/Spatial_Positioning/${sampleid}_cell_locations.csv \
        | grep -v Cell_Barcode > $outdir/spatial/${sampleid}/ROI_barcodes.txt

    if [[ $paired == "TRUE" ]];then
        $python $src/DBiT/snapatac2_preprocess.py \
            --input_bam $outdir/align/${sampleid}.sorted.bam \
            --output_dir $outdir/fragments \
            --sampleid $sampleid \
            --species $species \
            --max_iter 1 \
            --is_paired \
            --n_features 50000 \
            --min_fragments_n 1 \
            --tissue_barcode $outdir/spatial/${sampleid}/ROI_barcodes.txt
    else
        $python $src/DBiT/snapatac2_preprocess.py \
            --input_bam $outdir/align/${sampleid}.sorted.bam \
            --output_dir $outdir/fragments \
            --sampleid $sampleid \
            --species $species \
            --max_iter 1 \
            --n_features 50000 \
            --min_fragments_n 1 \
            --tissue_barcode $outdir/spatial/${sampleid}/ROI_barcodes.txt
    fi
    pigz -d $outdir/fragments/${sampleid}_CB.fragments.tsv.gz 
    sort -k1,1V -k2,2n \
        -T . $outdir/fragments/${sampleid}_CB.fragments.tsv | awk '{print $1, $2, $3, $4, $5}' OFS="\t" > $outdir/fragments/${sampleid}.fragments.tsv
    bgzip $outdir/fragments/${sampleid}.fragments.tsv

    tabix -p bed $outdir/fragments/${sampleid}.fragments.tsv.gz
    pigz $outdir/fragments/${sampleid}_CB.fragments.tsv

    $Rscript $signac_chromVar \
        $outdir/fragments/${sampleid}.peak_matrix.h5 \
        $outdir/fragments/${sampleid}.fragments.tsv.gz \
        $species

    ls $outdir/fragments/*metadata.txt | while read meta_file;do 
        prefix_sample=$(basename -s ".metadata.txt" $meta_file)
        $Rscript $src/atac_QC.R 100 $meta_file $outdir/fragments/$prefix_sample.qc.pdf
    done 

}


function spatial_position(){
    if [[ ! -d $outdir/spatial ]];then 
        mkdir -p $outdir/spatial
    fi 

    $python $create_gene_activity \
        --input_h5ad $outdir/fragments/${sampleid}_CB.h5ad \
        --species $species

    $Rscript $createSeuratObj \
        --input_h5ad $outdir/fragments/${sampleid}_CB.gene_activity.h5ad \
        --metadata_file $outdir/fragments/${sampleid}_CB.metadata.txt \
        --img $outdir/${sampleid}/Outs/${sampleid}_aligned_DAPI.png \
        --position_file $outdir/${sampleid}/Outs/${sampleid}_filtered_feature_bc_matrix/cell_locations.tsv.gz \
        --outdir $outdir/spatial

}


function bulk(){
    ### calculate as bulk samples
    if [[ ! -d $outdir/bulk ]];then 
        mkdir -p $outdir/bulk
    fi  

   $picard MarkDuplicates \
        -I $outdir/align/${sampleid}.sorted.bam \
        -O $outdir/bulk/${sampleid}.rmdup.bam \
        -M $outdir/bulk/${sampleid}.dedup.txt \
        -REMOVE_DUPLICATES true \
        -READ_NAME_REGEX null \
        --TMP_DIR ${outdir}/bulk
    $samtools index -@ $SLURM_CPUS_PER_TASK $outdir/bulk/${sampleid}.rmdup.bam

    ## filter the chrM reads
    $samtools idxstats -@ $SLURM_CPUS_PER_TASK $outdir/bulk/${sampleid}.rmdup.bam \
        | cut -f 1 | grep -v MT \
        | xargs $samtools view -b $outdir/bulk/${sampleid}.rmdup.bam > $outdir/bulk/${sampleid}.clean.bam
    
    $samtools index -@ $SLURM_CPUS_PER_TASK $outdir/bulk/${sampleid}.clean.bam

    effect_genome_size=2652783500
    
    $bamCoverage --numberOfProcessors $SLURM_CPUS_PER_TASK \
        --binSize 10 \
        --normalizeUsing RPGC \
        --effectiveGenomeSize $effective_genome_size \
        --bam $outdir/bulk/${sampleid}.clean.bam \
        -o $outdir/bulk/${sampleid}.clean.bw

}

function stat(){

    if [[ ! -d $outdir/stat ]];then 
        mkdir -p $outdir/stat 
    fi 
    
    valid_barcode_reads=$(grep valid ${outdir}/${sampleid}/Analysis/scRNA-seq_Analysis/${sampleid}_summary.json | awk '{print $NF}' | sed 's/,//')
    total_raw_reads=$(grep total ${outdir}/${sampleid}/Analysis/scRNA-seq_Analysis/${sampleid}_summary.json | awk '{print $NF}' | sed 's/,//')

    valid_barcode_percent=$(echo "scale=2; $valid_barcode_reads / $total_raw_reads * 100" | bc | xargs)
    reads_passing_filter=$(grep 'output":' $outdir/trim/${sampleid}.trimming_report.txt | awk '{print $2}' | head -n 1 | sed 's/,//g')
    reads_passing_filter_rate=$(echo "scale=2; $reads_passing_filter / $valid_barcode_reads * 100" | bc | xargs)
    
    unique_align_reads=$($samtools view -q 20 -c -@ 12 $outdir/align/${sampleid}.sorted.bam | xargs)

    align_rate=$(echo "scale=2; $unique_align_reads / $reads_passing_filter * 100" | bc | xargs)

    if [[ ! -s ${outdir}/stat/stat.txt ]]; then
        echo -e "SampleID\tSpecies\tTotal_reads\tBarcode_valid_percent\tReads_high_qualtity_percent\tUnique_aligned_reads\tUnique_align_rate" > ${outdir}/stat/stat.txt
    fi

    echo -e "$sampleid\t${species}\t" \
        "${total_raw_reads}\t" \
        "${valid_barcode_percent}%\t" \
        "${reads_passing_filter_rate}%\t" \
        "${unique_align_reads}%\t" \
        "${align_rate}\t" \
        >>${outdir}/stat/stat.txt

}




set -x 


case "$option" in
pipe)
    # #stat_abortive
    run_seekSpace
    trim_fq
    align 
    snapATAC2
    spatial_position
    stat
    bulk
    ;;
align)
    trim_fq
    align
    ;;
QC)
    stat_abortive
    ;;
atac)
    snapATAC2
    ;;
dedup)
    dedup
    ;;
bulk)
    bulk
    # trim_Nextera
    ;;
dedup_UMI)
    dedup_UMI
    ;;
spatial)
    spatial_position
    ;;
no_IVT)
    paired=FALSE
    barcode_calling
    trim_Nextera
    paired=TRUE
    snapATAC2
    bulk
    spatial_position
    ;;
*)
    echo "unknown parameter:" $1
    helpdoc
    exit 1
    ;;
esac