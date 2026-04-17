# Spatial Epigenome LIANTI-style Pipelines

This repository contains shell-centered pipelines and helper scripts for spatial epigenomics method development. The core experimental idea is to use LIANTI-like linear DNA amplification to improve the low resolution and low detection rate that often limit spatial epigenome assays.

The code was reorganized from a development workspace into a GitHub-ready layout. Raw sequencing data, BAM files, h5ad objects, and run outputs are intentionally excluded by `.gitignore`.

## Repository Layout

```text
.
├── pipelines/                 # Main shell entrypoints
├── scripts/                   # Python/R helper scripts called by pipelines
│   ├── DBiT/                  # DBiT-specific helpers
│   ├── RongFan_sp_atac/       # RongFan/DBiT variant helper
│   ├── slide_tag_atac/        # Slide-tags helper scripts
│   └── visium_hd_atac/        # Visium HD helper scripts
├── tools/                     # Standalone QC/conversion utilities
├── resources/                 # Small barcode/probe resources copied with the repo
├── config/                    # Environment and zUMIs templates
├── lib/                       # Shared shell helpers
└── docs/                      # Script index and maintenance notes
```

## Main Entrypoints

| Entrypoint | Main use case |
| --- | --- |
| `pipelines/st_lianti_atac.pipe.sh` | LIANTI-style spatial ATAC pipeline for BMK/S1000/S3000-style chips. |
| `pipelines/st_DBiT_T7_atac.pipe.sh` | DBiT T7/IVT spatial ATAC pipeline. |
| `pipelines/st_visium_hd_atac.pipe.sh` | 10x Visium HD spatial ATAC/IVT analysis. |
| `pipelines/st_slide_tags_atac.pipe.sh` | Slide-tags spatial ATAC workflow. |
| `pipelines/rongfan_sp_atac.pipe.sh` | RongFan spatial ATAC / DBiT variant workflow. |
| `pipelines/dbit_multiome_rna.pipe.sh` | DBiT multiome RNA helper workflow using zUMIs. |
| `pipelines/st_lianti_atac.step.sh` | Older SLURM step-by-step helper for LIANTI development runs. |

## Configure the Environment

Create a local config file and edit all machine-specific paths:

```bash
cp config/pipeline.env.example config/pipeline.env
vim config/pipeline.env
```

Then source it before running pipelines:

```bash
source config/pipeline.env
```

The shell entrypoints also auto-load `config/pipeline.env` when it exists. You can use a different config file by setting:

```bash
export SPATIAL_EPIGENOME_CONFIG=/path/to/my.pipeline.env
```

Important variables include:

- Tool executables: `CUTADAPT`, `BOWTIE2`, `SAMTOOLS`, `BEDTOOLS`, `PICARD`, `MACS2`, `RSCRIPT`, `PYTHON`.
- References: `HG38_BOWTIE2_INDEX`, `MM10_BOWTIE2_INDEX`, `HG38_CHROM_SIZES`, `MM10_CHROM_SIZES`, `HG38_GFF3`, `MM10_GFF3`.
- External method tools: `SPATIAL_IVT`, `MERGE_SUPER_SPOT`, `BST_SRC_S1000`, `BST_SRC_S3000`, `SEEKSPACE_TOOLS`, `ZUMIS`.

## Minimal Run Examples

LIANTI-style spatial ATAC:

```bash
bash pipelines/st_lianti_atac.pipe.sh \
  --read1 data/sample_R1.fastq.gz \
  --read2 data/sample_R2.fastq.gz \
  --sampleid sample01 \
  --species mm10 \
  --decode_file resources/DBiT/spatial_barcodes.txt \
  --outdir results/sample01 \
  --option pipe
```

DBiT T7 spatial ATAC:

```bash
bash pipelines/st_DBiT_T7_atac.pipe.sh \
  --read1 data/sample_R1.fastq.gz \
  --read2 data/sample_R2.fastq.gz \
  --sampleid dbit01 \
  --species mm10 \
  --decode_file resources/DBiT/spatial_barcodes.txt \
  --outdir results/dbit01 \
  --option pipe
```

Visium HD spatial ATAC:

```bash
bash pipelines/st_visium_hd_atac.pipe.sh \
  --read1 data/sample_R1.fastq.gz \
  --read2 data/sample_R2.fastq.gz \
  --sampleid visium01 \
  --species mm10 \
  --spaceranger_bam /path/to/spaceranger/possorted_genome_bam.bam \
  --outdir results/visium01 \
  --option pipe
```

Slide-tags spatial ATAC:

```bash
bash pipelines/st_slide_tags_atac.pipe.sh \
  --atac_read1 data/atac_R1.fastq.gz \
  --atac_read2 data/atac_R2.fastq.gz \
  --spatial_r1 data/spatial_R1.fastq.gz \
  --spatial_r2 data/spatial_R2.fastq.gz \
  --hdmi_read data/hdmi.fastq.gz \
  --sampleid slide01 \
  --species mm10 \
  --outdir results/slide01 \
  --option pipe
```

## Common `--option` Values

Most entrypoints support an `--option` flag to run a full pipeline or a single stage:

- `pipe`: run the default full workflow for that platform.
- `QC`: run abortive-extension or QC-related summaries when available.
- `align`: run trimming/alignment stages.
- `atac`: run ATAC object creation and downstream analysis stages.
- `bulk`: make pseudobulk BAM/bigWig/peaks where supported.
- `spatial`: build spatial coordinate outputs or spatial Seurat/Signac objects.
- `stat`: collect summary metrics where supported.

Check the bottom `case "$option"` block in each entrypoint for the exact stages enabled by that script.

## Expected Inputs

Typical inputs are:

- Paired FASTQ files with `R1` and `R2` in the filenames.
- A species label: `mm10` or `hg38`.
- A platform-specific barcode decode file, usually from `resources/DBiT/` or a chip vendor.
- Optional image or spatial files for HE/DAPI/spaceranger/SeekSpace workflows.
- Reference indexes and annotations configured in `config/pipeline.env`.

## Main Outputs

Output directories are pipeline-specific, but common subdirectories are:

- `barcode/`: barcode assignment statistics and barcoded FASTQs.
- `trim/`: cutadapt/trim_galore outputs.
- `align/`: aligned BAM files and indexes.
- `fragments/`: SnapATAC2 fragments, metadata, and h5ad objects.
- `spatial/`: spatial barcodes, Seurat/Signac objects, and spatial summaries.
- `bulk/`: pseudobulk BAM, peaks, and bigWig tracks.
- `stat/`: QC tables such as FRiP, fragment length, and spot-swap estimates.

## Utility Scripts

```bash
# Estimate spot-swap rate from 6-column fragments.
bash tools/cal_swap_score.sh results/sample/fragments/sample.fragments.tsv.gz

# Remove likely spot-swap fragments with local/remote filtering.
bash tools/remove_spot_swap.sh \
  results/sample/fragments/sample.fragments.tsv.gz \
  results/sample/fragments/sample.cleaned.fragments.tsv.gz \
  16 64 0.1

# Convert fragments to bigWig coverage.
bash tools/convertFrag2bw.sh \
  results/sample/fragments/sample.fragments.tsv.gz \
  results/sample/bulk/sample \
  /path/to/mm10.chrom.sizes
```

## Usage Tips

1. Start with `--option QC` or a small FASTQ subset before running `--option pipe` on a full dataset.
2. Keep one `config/pipeline.env` per cluster or conda environment. Do not commit it.
3. Use absolute paths in `config/pipeline.env`, but keep committed scripts path-free.
4. For SLURM jobs, set `SLURM_CPUS_PER_TASK` explicitly; otherwise scripts default to 4 threads.
5. Confirm barcode orientation and decode-file chemistry before interpreting spatial dropouts.
6. Inspect `stat/`, cutadapt reports, and barcode statistics before downstream clustering.
7. For spot-swap analysis, confirm fragment files are 6-column SE-style fragments with barcode coordinates encoded in the barcode name.
8. For public releases, keep raw FASTQ/BAM/h5ad files outside the repository.

## Development Notes

- New Python functions should use type hints and logging.
- Shell scripts should keep platform-specific paths in `config/pipeline.env`.
- R scripts should prefer CLI arguments or environment variables over hard-coded paths.
- Add the smallest meaningful validation run when changing logic.

## Quick Validation

Run syntax checks after editing:

```bash
bash -n pipelines/st_lianti_atac.pipe.sh
bash -n pipelines/st_DBiT_T7_atac.pipe.sh
bash -n pipelines/st_visium_hd_atac.pipe.sh
bash -n pipelines/st_slide_tags_atac.pipe.sh
python -m py_compile scripts/create_gene_activity.py scripts/cal_length.py
```
