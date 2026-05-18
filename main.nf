#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ============================================================
// Gene Family Identification Pipeline (Plant Genome/Proteome)
// ============================================================
// Steps:
//   1. HMMER  – HMM-based candidate identification
//              + isoform deduplication (best score per gene)
//   2. InterProScan – domain confirmation
//   3. BLAST RBH    – orthology assignment (no filtering)
//   4. MUSCLE       – multiple sequence alignment (+ optional outgroup)
//   5. TrimAl       – alignment trimming
//   6. IQ-TREE3     – maximum-likelihood phylogeny (+ optional outgroup rooting)
// ============================================================

log.info """
╔══════════════════════════════════════════════════════════╗
║        Gene Family Identification Pipeline               ║
║        Plant Genome / Proteome Analysis                  ║
╚══════════════════════════════════════════════════════════╝
  proteome      : ${params.proteome}
  hmm_profiles  : ${params.hmm_profiles}
  reference_db  : ${params.reference_db}
  outgroup      : ${params.outgroup ?: "not set"}
  outgroup_taxon: ${params.outgroup_taxon ?: "not set (unrooted tree)"}
  target_domains: ${params.target_domains}
  outdir        : ${params.outdir}
  threads       : ${params.threads}
""".stripIndent()

// ── Module imports ─────────────────────────────────────────────────────────
// FIX: imports now reference root-level .nf files (no modules/ subdirectory)
// FIX: EXTRACT_SEQUENCES and MAKEBLASTDB aliased because DSL2 forbids calling
//      the same process more than once without aliasing.
include { HMMER_SEARCH                              } from './hmmer'
include { FILTER_HMMER                              } from './hmmer'
include { INTERPROSCAN                              } from './interproscan'
include { CONFIRM_CANDIDATES                        } from './interproscan'
include { MAKEBLASTDB as MAKEBLASTDB_SUBJECT        } from './blast'
include { MAKEBLASTDB as MAKEBLASTDB_REF            } from './blast'
include { BLAST_FORWARD                             } from './blast'
include { BLAST_REVERSE                             } from './blast'
include { RBH_ORTHOLOGY                             } from './blast'
include { EXTRACT_SEQUENCES as EXTRACT_HMMER_SEQS   } from './utils'
include { EXTRACT_SEQUENCES as EXTRACT_CONFIRM_SEQS } from './utils'
include { EXTRACT_SEQUENCES as EXTRACT_FINAL_SEQS   } from './utils'
include { MUSCLE_ALIGN                              } from './muscle'
include { TRIMAL_TRIM                               } from './trimal'
include { IQTREE3                                   } from './iqtree'
include { SUMMARY_REPORT                            } from './utils'

// ── Workflow ───────────────────────────────────────────────────────────────
workflow {

    // ── Input channels ──────────────────────────────────────────────────
    proteome_ch      = Channel.fromPath(params.proteome,     checkIfExists: true)
    hmm_profiles_ch  = Channel.fromPath(params.hmm_profiles, checkIfExists: true)
    reference_db_ch  = Channel.fromPath(params.reference_db, checkIfExists: true)

    // Optional outgroup FASTA (empty file passed if not set)
    outgroup_ch = params.outgroup
        ? Channel.fromPath(params.outgroup, checkIfExists: true)
        : Channel.fromPath("$projectDir/assets/empty.fasta", checkIfExists: false)
              .ifEmpty { file("$projectDir/assets/empty.fasta") }

    // ── Step 1: HMMER + isoform deduplication ───────────────────────────
    HMMER_SEARCH(proteome_ch, hmm_profiles_ch.collect())
    FILTER_HMMER(HMMER_SEARCH.out.domtblout)

    // ── Step 2: Extract HMMER candidate sequences ────────────────────────
    EXTRACT_HMMER_SEQS(
        proteome_ch,
        FILTER_HMMER.out.candidate_ids,
        "hmmer_candidates"
    )

    // ── Step 3: InterProScan – domain confirmation ───────────────────────
    INTERPROSCAN(EXTRACT_HMMER_SEQS.out.sequences)
    CONFIRM_CANDIDATES(
        INTERPROSCAN.out.tsv,
        FILTER_HMMER.out.candidate_ids,
        params.target_domains
    )

    // Extract confirmed sequences
    EXTRACT_CONFIRM_SEQS(
        proteome_ch,
        CONFIRM_CANDIDATES.out.confirmed_ids,
        "confirmed_candidates"
    )
    confirmed_seqs_ch = EXTRACT_CONFIRM_SEQS.out.sequences

    // ── Step 4: BLAST RBH – orthology ASSIGNMENT only ───────────────────
    // All confirmed sequences proceed downstream regardless of RBH status.
    // RBH result is annotation only — no sequences removed.
    MAKEBLASTDB_SUBJECT(proteome_ch,     "subject_db")
    MAKEBLASTDB_REF(reference_db_ch,     "reference_db")

    BLAST_FORWARD(
        confirmed_seqs_ch,
        MAKEBLASTDB_REF.out.db
    )
    BLAST_REVERSE(
        reference_db_ch,
        MAKEBLASTDB_SUBJECT.out.db
    )
    RBH_ORTHOLOGY(
        BLAST_FORWARD.out.results,
        BLAST_REVERSE.out.results,
        CONFIRM_CANDIDATES.out.confirmed_ids
    )

    // Use ALL confirmed sequences downstream (RBH_ORTHOLOGY does not filter)
    EXTRACT_FINAL_SEQS(
        proteome_ch,
        RBH_ORTHOLOGY.out.all_ids,
        "final_candidates"
    )
    final_seqs_ch = EXTRACT_FINAL_SEQS.out.sequences

    // ── Step 5: MUSCLE – MSA with optional outgroup ─────────────────────
    combined_ch = final_seqs_ch.combine(reference_db_ch)
    MUSCLE_ALIGN(combined_ch, outgroup_ch)

    // ── Step 6: TrimAl ──────────────────────────────────────────────────
    TRIMAL_TRIM(MUSCLE_ALIGN.out.alignment)

    // ── Step 7: IQ-TREE3 – ML tree with optional outgroup rooting ───────
    // --outgroup_taxon must match a sequence ID present in the alignment.
    IQTREE3(TRIMAL_TRIM.out.trimmed_aln)

    // ── Summary ─────────────────────────────────────────────────────────
    SUMMARY_REPORT(
        FILTER_HMMER.out.candidate_ids,
        CONFIRM_CANDIDATES.out.confirmed_ids,
        RBH_ORTHOLOGY.out.all_ids,
        IQTREE3.out.treefile
    )
}

workflow.onComplete {
    log.info ( workflow.success
        ? "\n Pipeline complete! Results: ${params.outdir}"
        : "\n Pipeline failed. Check logs above." )
}
