// modules/utils.nf
// Utility processes: sequence extraction and summary report

process EXTRACT_SEQUENCES {
    tag "extract_${label}"
    label 'low_cpu'
    publishDir "${params.outdir}/sequences", mode: 'copy'

    input:
    path  proteome
    path  id_list
    val   label

    output:
    path "${label}.fasta",         emit: sequences
    path "${label}_not_found.txt", emit: not_found

    script:
    """
    #!/usr/bin/env python3
    import sys

    with open("${id_list}") as fh:
        target_ids = set(line.strip() for line in fh if line.strip())

    found     = set()
    not_found = set(target_ids)
    writing   = False

    with open("${proteome}") as fh, \
         open("${label}.fasta", "w") as out:
        for line in fh:
            if line.startswith(">"):
                seq_id  = line[1:].split()[0]
                writing = seq_id in target_ids
                if writing:
                    found.add(seq_id)
                    not_found.discard(seq_id)
                    out.write(line)
            elif writing:
                out.write(line)

    with open("${label}_not_found.txt", "w") as out:
        for sid in sorted(not_found):
            out.write(sid + "\n")

    if not_found:
        print(f"WARNING: {len(not_found)} IDs not found in proteome", file=sys.stderr)
        for sid in sorted(not_found):
            print(f"  Missing: {sid}", file=sys.stderr)

    print(f"Extracted {len(found)}/{len(target_ids)} sequences for '{label}'", file=sys.stderr)

    if len(found) == 0:
        print("ERROR: No sequences extracted. Check ID format.", file=sys.stderr)
        sys.exit(1)
    """
}


process SUMMARY_REPORT {
    tag "summary_report"
    label 'low_cpu'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path hmmer_ids       // IDs after HMMER + isoform dedup
    path confirmed_ids   // IDs after InterProScan domain confirmation
    path final_ids       // IDs from RBH_ORTHOLOGY (= all confirmed, no filtering)
    path treefile        // IQ-TREE3 treefile

    output:
    path "pipeline_summary.txt",  emit: summary
    path "final_gene_ids.txt",    emit: final_ids

    script:
    """
    #!/usr/bin/env python3
    import os
    from datetime import datetime

    def count_lines(filepath):
        try:
            with open(filepath) as fh:
                return sum(1 for line in fh if line.strip())
        except Exception:
            return 0

    def read_ids(filepath):
        try:
            with open(filepath) as fh:
                return set(line.strip() for line in fh if line.strip())
        except Exception:
            return set()

    n_hmmer     = count_lines("${hmmer_ids}")
    n_confirmed = count_lines("${confirmed_ids}")
    n_final     = count_lines("${final_ids}")

    final_set = read_ids("${final_ids}")
    tree_ok   = os.path.exists("${treefile}") and os.path.getsize("${treefile}") > 10

    with open("final_gene_ids.txt", "w") as out:
        for sid in sorted(final_set):
            out.write(sid + "\n")

    def pct(a, b):
        return f"{a/b*100:.1f}%" if b > 0 else "N/A"

    with open("pipeline_summary.txt", "w") as out:
        out.write("=" * 60 + "\n")
        out.write("  Gene Family Identification Pipeline — Summary\n")
        out.write("=" * 60 + "\n")
        out.write(f"  Date : {datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n")
        out.write("  Step                              Sequences\n")
        out.write("  " + "-" * 50 + "\n")
        out.write(f"  1. HMMER (post isoform dedup)   : {n_hmmer:>6}\n")
        out.write(f"  2. InterProScan confirmed        : {n_confirmed:>6}  ({pct(n_confirmed, n_hmmer)} retained)\n")
        out.write(f"  3. Orthology assigned (no filter): {n_final:>6}  (= all confirmed)\n")
        out.write("\n")
        out.write(f"  Final gene family members        : {n_final}\n")
        out.write(f"  Phylogenetic tree produced       : {'Yes' if tree_ok else 'No'}\n\n")
        out.write("  Note: RBH orthology assignments are in\n")
        out.write("  03_blast/orthology_table.tsv\n\n")
        out.write("  Output files:\n")
        out.write("    sequences/final_candidates.fasta\n")
        out.write("    04_alignment/sequence_manifest.tsv  (outgroup flags)\n")
        out.write("    06_phylogeny/phylogeny.treefile\n")
        out.write("    final_gene_ids.txt\n")
        out.write("=" * 60 + "\n")

    print("Pipeline complete. See pipeline_summary.txt")
    """
}
