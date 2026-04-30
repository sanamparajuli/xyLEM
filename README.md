# Gene Family Identification Pipeline

A **Nextflow DSL2** pipeline for genome-wide identification of plant gene family
members using a multi-evidence approach.

```
Proteome (FASTA)
      │
      ▼
┌──────────────────────────┐
│  1. HMMER                │  hmmsearch – profile HMM scan
│     + isoform dedup      │  one representative per gene (best bitscore)
└────────────┬─────────────┘
             │ unique gene candidates
             ▼
┌──────────────────────────┐
│  2. InterProScan         │  domain confirmation (Pfam, PANTHER, Gene3D …)
└────────────┬─────────────┘
             │ confirmed candidates
             ▼
┌──────────────────────────┐
│  3. BLAST RBH            │  orthology ASSIGNMENT only — no sequences removed
│     (annotation step)    │  labels: RBH_ORTHOLOG / PUTATIVE_HOMOLOG / NO_BLAST_HIT
└────────────┬─────────────┘
             │ all confirmed candidates (+ orthology annotations)
             ▼
┌──────────────────────────┐
│  4. MUSCLE v5            │  multiple sequence alignment
│     + optional outgroup  │  outgroup sequences tracked in sequence_manifest.tsv
└────────────┬─────────────┘
             ▼
┌──────────────────────────┐
│  5. TrimAl               │  alignment trimming
└────────────┬─────────────┘
             ▼
┌──────────────────────────┐
│  6. IQ-TREE3             │  maximum-likelihood phylogeny + UFBoot
│     (unrooted by default)│  root manually in FigTree / iTOL / ggtree
└──────────────────────────┘
```

---

## Requirements

| Tool          | Version tested | Install via               |
|---------------|---------------|---------------------------|
| Nextflow      | ≥ 23.04       | `curl -s get.nextflow.io` |
| HMMER         | 3.4           | bioconda                  |
| InterProScan  | 5.67-99.0     | bioconda / manual         |
| BLAST+        | 2.15.0        | bioconda                  |
| MUSCLE        | 5.1           | bioconda                  |
| TrimAl        | 1.4.1         | bioconda                  |
| IQ-TREE3      | 3.0.0         | bioconda                  |
| Python        | ≥ 3.10        | conda / system            |
| Biopython     | ≥ 1.83        | bioconda                  |

---

## Quick Start

### 1. Install Nextflow
```bash
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/
```

### 2. Prepare inputs

```
data/
├── target_proteome.fasta      # Your plant proteome
├── profiles/
│   └── gene_family.hmm        # HMM profile(s) for your gene family
├── reference_proteins.fasta   # Known members (e.g. from Arabidopsis TAIR)
└── outgroup.fasta             # Optional: outgroup sequences for phylogenetics
```

### 3. Run (conda)
```bash
nextflow run main.nf \
    -profile conda \
    --proteome        data/target_proteome.fasta \
    --hmm_profiles    "data/profiles/*.hmm" \
    --reference_db    data/reference_proteins.fasta \
    --target_domains  "PF03106" \
    --outdir          results/
```

### 4. Run with outgroup
```bash
nextflow run main.nf \
    -profile conda \
    --proteome        data/target_proteome.fasta \
    --hmm_profiles    "data/profiles/*.hmm" \
    --reference_db    data/reference_proteins.fasta \
    --target_domains  "PF03106" \
    --outgroup        data/outgroup.fasta \
    --outdir          results/
```
The outgroup sequences are included in the alignment and labelled in
`04_alignment/sequence_manifest.tsv`. The tree is left **unrooted** by default —
root it at the outgroup in your preferred viewer (see [Tree visualisation](#tree-visualisation)).

### 5. Run on SLURM
```bash
nextflow run main.nf -profile slurm,conda \
    --proteome ... [other params]
```

---

## Parameters

| Parameter            | Default        | Description                                         |
|----------------------|---------------|-----------------------------------------------------|
| `--proteome`         | required      | Target plant proteome FASTA                         |
| `--hmm_profiles`     | required      | HMM profile file(s), glob OK                        |
| `--reference_db`     | required      | Reference known gene family members FASTA           |
| `--target_domains`   | required      | Pfam/IPR accessions (comma-separated)               |
| `--outgroup`         | `null`        | FASTA of outgroup sequence(s) — optional            |
| `--interpro_db`      | system default| InterProScan local database path                    |
| `--hmmer_evalue`     | `1e-5`        | E-value cutoff for HMMER                            |
| `--hmmer_coverage`   | `0.5`         | Min fraction of HMM profile covered                 |
| `--blast_evalue`     | `1e-5`        | E-value cutoff for BLAST                            |
| `--blast_identity`   | `30.0`        | Min % identity for BLAST hits                       |
| `--blast_coverage`   | `50.0`        | Min query coverage % for BLAST hits                 |
| `--iqtree_model`     | `TEST`        | Substitution model (`TEST` = auto-select)           |
| `--iqtree_bootstrap` | `1000`        | UFBoot replicates                                   |
| `--trimal_method`    | `automated1`  | TrimAl strategy: `automated1`, `gappyout`, `strict` |
| `--threads`          | `8`           | Default CPU threads per process                     |
| `--outdir`           | `results`     | Output directory                                    |

---

## Pipeline design notes

### Isoform deduplication (Step 1)
Many plant proteomes contain multiple isoforms per gene (e.g. `AT1G01010.1`,
`AT1G01010.2`). The HMMER filter step automatically detects and strips isoform
suffixes (`.N`, `.N.N`, `-N`) to group hits by gene, then retains only the
isoform with the **highest bitscore** (tie-break: lowest E-value, then highest
HMM coverage). The selected representative and its gene ID are recorded in
`01_hmmer/hmmer_filtered.tsv`.

### RBH as orthology annotation (Step 3)
BLAST reciprocal best hits are used to **annotate** sequences, not filter them.
Every sequence confirmed by InterProScan proceeds to alignment regardless of its
BLAST result. Each sequence receives one of three labels in
`03_blast/orthology_table.tsv`:

| Label              | Meaning                                             |
|--------------------|-----------------------------------------------------|
| `RBH_ORTHOLOG`     | Strict reciprocal best hit with a reference member  |
| `PUTATIVE_HOMOLOG` | Significant forward hit to reference, but no RBH    |
| `NO_BLAST_HIT`     | No significant BLAST hit (possible lineage-specific)|

### Outgroup handling (Steps 4–6)
If `--outgroup` is provided, outgroup sequences are merged into the alignment
and flagged in `04_alignment/sequence_manifest.tsv` (`is_outgroup = true`).
The phylogenetic tree is produced **unrooted** — this is standard practice and
gives you maximum flexibility to root the tree in the tool of your choice.

---

## Outputs

```
results/
├── 01_hmmer/
│   ├── hmmer_results.domtblout
│   ├── hmmer_filtered.tsv          ← includes gene_id + selected isoform
│   └── candidate_ids.txt
├── 02_interproscan/
│   ├── interproscan_results.tsv
│   ├── confirmed_ids.txt
│   └── domain_confirmation_summary.tsv
├── 03_blast/
│   ├── blast_forward.tsv
│   ├── blast_reverse.tsv
│   ├── orthology_table.tsv         ← RBH_ORTHOLOG / PUTATIVE_HOMOLOG / NO_BLAST_HIT
│   └── rbh_pairs.tsv
├── 04_alignment/
│   ├── combined_aligned.fasta
│   └── sequence_manifest.tsv       ← source (candidate/reference/outgroup) per seq
├── 05_trimming/
│   ├── trimmed_alignment.fasta
│   └── trimmed_alignment.html
├── 06_phylogeny/
│   ├── phylogeny.treefile           ← Main result: unrooted ML tree
│   ├── phylogeny.contree            ← Consensus tree
│   └── phylogeny.iqtree             ← Full IQ-TREE log
├── sequences/
│   ├── hmmer_candidates.fasta
│   ├── confirmed_candidates.fasta
│   └── final_candidates.fasta
├── pipeline_summary.txt
└── final_gene_ids.txt
```

---

## Finding HMM profiles and domain accessions

Download from **InterPro** (https://www.ebi.ac.uk/interpro) by searching your
gene family. Common examples:

| Gene family    | Pfam accession |
|----------------|---------------|
| WRKY TF        | `PF03106`     |
| MYB TF         | `PF00249`     |
| HSP70          | `PF00012`     |
| Kinase domain  | `PF00069`     |
| LRR domain     | `PF00560`     |

Build your own HMM from a curated alignment:
```bash
hmmbuild my_family.hmm reference_alignment.fasta
```

---

## Resume a failed run

```bash
nextflow run main.nf -resume [same params]
```

---

## Test run

```bash
nextflow run main.nf -profile test,conda
```

---

## Citation

If you use this pipeline, please cite the underlying tools:

- **HMMER**: Eddy SR (2011) *PLoS Comput Biol* 7(10):e1002195
- **InterProScan**: Jones P et al. (2014) *Bioinformatics* 30(9):1236–1240
- **BLAST+**: Camacho C et al. (2009) *BMC Bioinformatics* 10:421
- **MUSCLE**: Edgar RC (2022) *Nature Methods* 19:714–717
- **TrimAl**: Capella-Gutierrez S et al. (2009) *Bioinformatics* 25(15):1972–1973
- **IQ-TREE3**: Wong TK et al. (2025)
- **Nextflow**: Di Tommaso P et al. (2017) *Nat Biotechnol* 35:316–319
