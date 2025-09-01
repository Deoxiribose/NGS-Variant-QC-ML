# NGS Variant Calling & Deep Learning QC

This project demonstrates an end-to-end workflow for variant calling and quality control (QC) using machine learning and deep learning.  
It includes a lightweight QuickStart demo that runs entirely in Python, as well as a realistic pipeline using the human mitochondrial genome (chrM) for alignment and variant calling.

---

## QuickStart (Python-only)

The QuickStart path runs on synthetic VCF data and does not require external bioinformatics tools.  
It shows how variant features such as QUAL, DP, and MQ can be extracted and used in machine learning models.

- Random Forest (scikit-learn) as a baseline model  
- Keras MLP (TensorFlow) as a deep learning tabular model  
- Evaluation with ROC AUC, precision–recall, and calibration

To run the QuickStart:

```bash
pip install -r requirements.txt
jupyter notebook NGS_Variant_QC_complete_v2.ipynb
```

---

## Realistic Pipeline (chrM)

For realism, the repository also documents a small but real NGS pipeline using the mitochondrial genome (chrM).  
This keeps runtimes manageable while following familiar steps:

1. Download chrM reference (hg38)  
2. Simulate paired-end reads with `wgsim`  
3. Align reads with BWA-MEM  
4. Sort and index BAM files with Samtools  
5. Call variants with bcftools  
6. Parse the resulting VCF into features for machine learning  

These commands are included in the notebook and require Linux or WSL with `bwa`, `samtools`, `bcftools`, `tabix`, `wgsim`, and `aria2` installed.

```bash
mkdir -p data && cd data
aria2c -x 8 -s 8 -o chrM.fa.gz https://hgdownload.soe.ucsc.edu/goldenpath/hg38/chromosomes/chrM.fa.gz
gunzip -f chrM.fa.gz
samtools faidx chrM.fa
bwa index chrM.fa

# Simulate a small paired-end dataset
wgsim -N 5000 -1 100 -2 100 chrM.fa sample_R1.fastq sample_R2.fastq
gzip -f sample_R1.fastq sample_R2.fastq

# Align and produce a sorted, indexed BAM
bwa mem -t 4 chrM.fa sample_R1.fastq.gz sample_R2.fastq.gz | samtools sort -o sample.chrM.sorted.bam
samtools index sample.chrM.sorted.bam

# Call variants and index the VCF
bcftools mpileup -f chrM.fa sample.chrM.sorted.bam -Ou | bcftools call -mv -Oz -o sample.chrM.vcf.gz
tabix -p vcf sample.chrM.vcf.gz
```

---

## Snakemake Workflow (Optional)

A minimal Snakemake workflow is included for reproducibility. It encodes the same steps as above.

Run it:

```bash
# Optional: prepare chrM reference and reads via the helper script
bash setup_chrM.sh

# Then execute the workflow
snakemake -j 4
```

Snakemake will generate the sorted BAM and compressed VCF (with tabix index) under `data/`.

---

## Repository Contents

- `NGS_Variant_QC_complete_v2.ipynb` — main Jupyter notebook (QuickStart and chrM pipeline)  
- `README.md` — project overview  
- `requirements.txt` — Python dependencies for QuickStart  
- `Snakefile` — Snakemake workflow for the chrM pipeline (optional)  
- `setup_chrM.sh` — helper script for chrM reference and read simulation  

---

## Results

- The QuickStart workflow generates a labeled dataset and trains both Random Forest and Keras MLP models.  
- Metrics such as ROC AUC and average precision illustrate model performance on synthetic variants.  
- This serves as a reproducible example of ML/DL applied to variant quality prediction and can be extended to larger genomes.

---

## Requirements

- Python 3.9 or newer  
- Python packages: `pandas`, `numpy`, `matplotlib`, `scikit-learn`, `tensorflow==2.16.1`

Install dependencies:

```bash
pip install -r requirements.txt
```

Optional system tools (for the chrM pipeline):  
`bwa`, `samtools`, `bcftools`, `tabix`, `wgsim`, `aria2`

---

## License

MIT License. You are free to use, modify, and share this project.
