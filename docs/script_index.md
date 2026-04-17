# Script Index

This index records how the GitHub-ready repository is organized around shell entrypoints.

## Pipeline Entrypoints

| Script | Calls or depends on | Notes |
| --- | --- | --- |
| `pipelines/st_lianti_atac.pipe.sh` | `scripts/snapatac2_preprocess.py`, `scripts/create_gene_activity.py`, `scripts/createSeuratObj.R`, `scripts/atac_QC.R`, `scripts/drawQC.R` | Main LIANTI-style spatial ATAC workflow. Requires BST/ROI tools configured in `config/pipeline.env`. |
| `pipelines/st_DBiT_T7_atac.pipe.sh` | `scripts/DBiT/extract_barcodes.py`, `scripts/DBiT/snapatac2_preprocess.py`, `scripts/DBiT/gene_activity.py`, `scripts/DBiT/createSeuratObj.R`, `scripts/DBiT/signac_chromVar.R` | DBiT T7/IVT spatial ATAC workflow. |
| `pipelines/st_visium_hd_atac.pipe.sh` | `scripts/visium_hd_atac/extract_barcodes.py`, `scripts/visium_hd_atac/snapatac2_preprocess.py`, `scripts/visium_hd_atac/createSeuratObj.R`, `scripts/visium_hd_atac/run_signac.R`, `scripts/cal_length.py` | Visium HD spatial ATAC/IVT workflow. |
| `pipelines/st_slide_tags_atac.pipe.sh` | `scripts/slide_tag_atac/createSeuratObj.R`, `scripts/DBiT/snapatac2_preprocess.py`, `scripts/DBiT/signac_chromVar.R`, `scripts/create_gene_activity.py` | Slide-tags spatial ATAC workflow. |
| `pipelines/rongfan_sp_atac.pipe.sh` | `scripts/RongFan_sp_atac/BC_process.py`, `scripts/DBiT/extract_barcodes.py`, `scripts/DBiT/createSeuratObj.R` | RongFan/DBiT variant workflow. |
| `pipelines/dbit_multiome_rna.pipe.sh` | `scripts/DBiT/extract_barcodes.py`, `scripts/DBiT/gene_activity.py`, `scripts/DBiT/svg2png.py`, `scripts/DBiT/makeMetadata.py`, `config/DBiT_multiome_RNA.yaml` | DBiT multiome RNA helper around zUMIs. |
| `pipelines/st_lianti_atac.step.sh` | `scripts/plot_abortive_distribution.R` | Older SLURM step helper retained for reproducibility. |

## Standalone Tools

| Script | Purpose |
| --- | --- |
| `tools/cal_ME.sh` | Split ME/Tn5-adaptor-positive and untrimmed reads with cutadapt. |
| `tools/cal_swap_score.sh` | Estimate spot-swap rate from 6-column fragments. |
| `tools/remove_spot_swap.sh` | Apply a two-pass hybrid local/remote filter to remove likely swapped fragments. |
| `tools/convertFrag2bw.sh` | Convert fragments to bedGraph and bigWig coverage tracks. |

## Resources

| Path | Purpose |
| --- | --- |
| `resources/DBiT/spatial_barcodes*.txt` | DBiT barcode coordinate resources used as decode files. |
| `resources/DBiT/DBiT.whitelist.txt` | Barcode whitelist used by DBiT RNA/zUMIs workflows. |
| `resources/visium_hd_atac/Visium_Mouse_Transcriptome_Probe_Set_v2.0_mm10-2020-A.csv` | Probe set used by Visium HD helper scripts. |

## Maintenance Checklist

1. Keep new platform-specific paths in `config/pipeline.env.example`, not in code.
2. Add usage examples to `README.md` when adding a new `--option` or entrypoint.
3. Run `bash -n` for changed shell scripts and `python -m py_compile` for changed Python scripts.
4. Keep raw FASTQ/BAM/h5ad/output files outside Git.
