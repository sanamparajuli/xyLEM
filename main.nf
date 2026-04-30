#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ============================================================
// Gene Family Identification Pipeline (Plant Genome/Proteome)
// ============================================================
// Steps:
//   1. HMMER        – HMM-based candidate identification
//                     + isoform deduplication (best score per gene)
//   2. InterProScan – domain confirmation
//   3. BLAST RBH    – orthology assignment (no filtering)
//   4. MUSCLE       – multiple sequence alignment (+ optional outgroup)
//   5. TrimAl       – alignment trimming
//   6. IQ-TREE3     – maximum-likelihood phylogeny (+ optional outgroup rooting)
// ============================================================

log.info """
╔══════════════════════════════════════════════════════════╗
║         Gene Family Identification Pipeline              ║
║         Plant Genome / Proteome Analysis                 ║
╚══════════════════════════════════════════════════════════╝
  proteome        : ${params.proteome}
  hmm_profiles    : ${params.hmm_profiles}
  reference_db    : ${params.reference_db}
  outgroup        : ${params.outgroup    ?: "not set"}
  outgroup_taxon  : ${params.outgroup_taxon ?: "not set (unrooted tree)"}
  target_domains  : ${params.target_domains}
  outdir          : ${params.outdir}
  threads         : ${params.threads}
""".stripIndent()

// ── Module imports ────────────────────────────────────────────
include { HMMER_SEARCH        } from './modules/hmmer'
include { FILTER_HMMER        } from './modules/hmmer'
include { INTERPROSCAN        } from './modules/interproscan'
include { CONFIRM_CANDIDATES  } from './modules/interproscan'
include { MAKEBLASTDB         } from './modules/blast'
include { BLAST_FORWARD       } from './modules/blast'
include { BLAST_REVERSE       } from './modules/blast'
include { RBH_ORTHOLOGY       } from './modules/blast'
include { EXTRACT_SEQUENCES   } from './modules/utils'
include { MUSCLE_ALIGN        } from './modules/muscle'
include { TRIMAL_TRIM         } from './modules/trimal'
include { IQTREE3             } from './modules/iqtree'
include { SUMMARY_REPORT      } from './modules/utils'

// ── Workflow ──────────────────────────────────────────────────
workflow {

    // ── Input channels ────────────────────────────────────────
    proteome_ch      = Channel.fromPath(params.proteome,     checkIfExists: true)
    hmm_profiles_ch  = Channel.fromPath(params.hmm_profiles, checkIfExists: true)
    reference_db_ch  = Channel.fromPath(params.reference_db, checkIfExists: true)

    // Optional outgroup FASTA (empty file passed if not set)
    outgroup_ch = params.outgroup
        ? Channel.fromPath(params.outgroup, checkIfExists: true)
        : Channel.fromPath("$projectDir/assets/empty.fasta", checkIfExists: false)
              .ifEmpty { file("$projectDir/assets/empty.fasta") }

    // ── Step 1: HMMER + isoform deduplication ────────────────
    HMMER_SEARCH(proteome_ch, hmm_profiles_ch.collect())
    FILTER_HMMER(HMMER_SEARCH.out.domtblout)
    // FILTER_HMMER now also deduplicates isoforms — one representative per gene

    // ── Step 2: Extract HMMER candidate sequences ─────────────
    EXTRACT_SEQUENCES(
        proteome_ch,
        FILTER_HMMER.out.candidate_ids,
        "hmmer_candidates"
    )

    // ── Step 3: InterProScan – domain confirmation ────────────
    INTERPROSCAN(EXTRACT_SEQUENCES.out.sequences)
    CONFIRM_CANDIDATES(
        INTERPROSCAN.out.tsv,
        FILTER_HMMER.out.candidate_ids,
        params.target_domains
    )

    // Extract confirmed sequences (all pass — no filtering at RBH step)
    EXTRACT_SEQUENCES(
        proteome_ch,
        CONFIRM_CANDIDATES.out.confirmed_ids,
        "confirmed_candidates"
    )
    confirmed_seqs_ch = EXTRACT_SEQUENCES.out.sequences

    // ── Step 4: BLAST RBH – orthology ASSIGNMENT only ─────────
    // All confirmed sequences proceed downstream regardless of RBH status.
    // RBH result is annotation attached to each sequence in the output table.
    MAKEBLASTDB(proteome_ch,     "subject_db")
    MAKEBLASTDB(reference_db_ch, "reference_db")

    BLAST_FORWARD(
        confirmed_seqs_ch,
        MAKEBLASTDB.out.db.filter { it[0] == "reference_db" }
    )
    BLAST_REVERSE(
        reference_db_ch,
        MAKEBLASTDB.out.db.filter { it[0] == "subject_db" }
    )

    RBH_ORTHOLOGY(
        BLAST_FORWARD.out.results,
        BLAST_REVERSE.out.results,
        CONFIRM_CANDIDATES.out.confirmed_ids   // all confirmed IDs pass through
    )

    // Use ALL confirmed sequences downstream (RBH_ORTHOLOGY does not filter)
    EXTRACT_SEQUENCES(
        proteome_ch,
        RBH_ORTHOLOGY.out.all_ids,
        "final_candidates"
    )
    final_seqs_ch = EXTRACT_SEQUENCES.out.sequences

    // ── Step 5: MUSCLE – MSA with optional outgroup ───────────
    combined_ch = final_seqs_ch.combine(reference_db_ch)

    MUSCLE_ALIGN(combined_ch, outgroup_ch)

    // ── Step 6: TrimAl ────────────────────────────────────────
    TRIMAL_TRIM(MUSCLE_ALIGN.out.alignment)

    // ── Step 7: IQ-TREE3 – ML tree with optional outgroup ─────
    // --outgroup_taxon must match a sequence ID present in the alignment.
    // If using an outgroup FASTA, set --outgroup_taxon to that sequence's ID.
    IQTREE3(TRIMAL_TRIM.out.trimmed_aln)

    // ── Summary ───────────────────────────────────────────────
    SUMMARY_REPORT(
        FILTER_HMMER.out.candidate_ids,
        CONFIRM_CANDIDATES.out.confirmed_ids,
        RBH_ORTHOLOGY.out.all_ids,
        IQTREE3.out.treefile
    )
}

workflow.onComplete {
    log.info ( workflow.success
        ? "\n✅  Pipeline complete! Results: ${params.outdir}"
        : "\n❌  Pipeline failed. Check logs above." )
}
