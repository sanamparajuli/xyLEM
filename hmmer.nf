// modules/hmmer.nf
// ── HMMER: profile-based sequence search ─────────────────────

process HMMER_SEARCH {
    tag "hmmsearch"
    label 'medium_cpu'
    publishDir "${params.outdir}/01_hmmer", mode: 'copy'

    input:
    path proteome        // Target plant proteome (FASTA)
    path hmm_profiles    // HMM profile file(s) – collected into one dir

    output:
    path "hmmer_results.domtblout", emit: domtblout
    path "hmmer_results.tblout",    emit: tblout
    path "hmmer_results.out",       emit: stdout

    script:
    """
    # If multiple profiles supplied, press and concatenate them
    if [ \$(ls *.hmm 2>/dev/null | wc -l) -gt 1 ]; then
        cat *.hmm > combined_profiles.hmm
    else
        cp *.hmm combined_profiles.hmm
    fi

    # Press the HMM database for faster searching
    hmmpress combined_profiles.hmm

    # Run hmmsearch
    hmmsearch \\
        --cpu ${task.cpus} \\
        -E ${params.hmmer_evalue} \\
        --domE ${params.hmmer_evalue} \\
        --domtblout hmmer_results.domtblout \\
        --tblout    hmmer_results.tblout \\
        combined_profiles.hmm \\
        ${proteome} \\
        > hmmer_results.out
    """
}

process FILTER_HMMER {
    tag "filter_hmmer"
    label 'low_cpu'
    publishDir "${params.outdir}/01_hmmer", mode: 'copy'

    input:
    path domtblout

    output:
    path "candidate_ids.txt",          emit: candidate_ids
    path "hmmer_filtered.tsv",         emit: filtered_table
    path "hmmer_filter_summary.txt",   emit: summary

    script:
    """
    #!/usr/bin/env python3
    import sys
    import re

    evalue_thresh   = float("${params.hmmer_evalue}")
    coverage_thresh = float("${params.hmmer_coverage}")

    # ── Step 1: parse all hits passing E-value + coverage thresholds ──
    all_hits = []
    with open("${domtblout}") as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.split()
            if len(cols) < 23:
                continue

            seq_id     = cols[0]
            hmm_name   = cols[3]
            hmm_len    = int(cols[5])
            seq_evalue = float(cols[6])
            dom_evalue = float(cols[11])
            bitscore   = float(cols[7])
            hmm_from   = int(cols[15])
            hmm_to     = int(cols[16])

            hmm_coverage = (hmm_to - hmm_from + 1) / hmm_len

            if seq_evalue <= evalue_thresh and hmm_coverage >= coverage_thresh:
                all_hits.append({
                    'seq_id':     seq_id,
                    'hmm_name':   hmm_name,
                    'seq_evalue': seq_evalue,
                    'dom_evalue': dom_evalue,
                    'bitscore':   bitscore,
                    'coverage':   hmm_coverage,
                })

    # ── Step 2: isoform deduplication ─────────────────────────────────
    # Gene ID is inferred by stripping common isoform suffixes:
    #   AT1G01010.1, AT1G01010.2  → gene = AT1G01010
    #   Os01g0100100-01            → gene = Os01g0100100
    #   Solyc01g005000.1.1         → gene = Solyc01g005000
    #   gene_id_transcript_id      → try splitting on last dot or dash+digits
    def get_gene_id(seq_id):
        # Try stripping trailing .N or .N.N (common in plant proteomes)
        m = re.match(r'^(.+?)(\\.\\d+)+\$', seq_id)
        if m:
            return m.group(1)
        # Try stripping trailing -NN (e.g. rice RAP-DB isoforms)
        m = re.match(r'^(.+?)-\\d+\$', seq_id)
        if m:
            return m.group(1)
        # No recognizable suffix — treat the whole ID as the gene ID
        return seq_id

    # Group hits by gene ID, keep the one with best bitscore
    # (tie-break: lowest e-value, then highest HMM coverage)
    gene_best = {}
    for hit in all_hits:
        gene_id = get_gene_id(hit['seq_id'])
        if gene_id not in gene_best:
            gene_best[gene_id] = hit
        else:
            prev = gene_best[gene_id]
            # Primary sort: bitscore descending
            if hit['bitscore'] > prev['bitscore']:
                gene_best[gene_id] = hit
            elif hit['bitscore'] == prev['bitscore']:
                # Secondary: evalue ascending
                if hit['seq_evalue'] < prev['seq_evalue']:
                    gene_best[gene_id] = hit
                elif hit['seq_evalue'] == prev['seq_evalue']:
                    # Tertiary: HMM coverage descending
                    if hit['coverage'] > prev['coverage']:
                        gene_best[gene_id] = hit

    # ── Step 3: write outputs ──────────────────────────────────────────
    candidates = {hit['seq_id']: hit for hit in gene_best.values()}

    with open("candidate_ids.txt", "w") as out:
        for sid in sorted(candidates):
            out.write(sid + "\\n")

    with open("hmmer_filtered.tsv", "w") as out:
        out.write("seq_id\\tgene_id\\thmm_name\\tseq_evalue\\tdom_evalue\\tbitscore\\thmm_coverage\\n")
        for gene_id, hit in sorted(gene_best.items()):
            out.write(
                f"{hit['seq_id']}\\t{gene_id}\\t{hit['hmm_name']}\\t"
                f"{hit['seq_evalue']:.2e}\\t{hit['dom_evalue']:.2e}\\t"
                f"{hit['bitscore']:.1f}\\t{hit['coverage']:.3f}\\n"
            )

    n_before = len(all_hits)
    n_after  = len(candidates)
    n_removed = n_before - n_after

    with open("hmmer_filter_summary.txt", "w") as out:
        out.write("HMMER Filter + Isoform Deduplication Summary\\n")
        out.write("=============================================\\n")
        out.write(f"E-value threshold        : {evalue_thresh}\\n")
        out.write(f"Coverage threshold       : {coverage_thresh}\\n")
        out.write(f"Raw hits passing filters : {n_before}\\n")
        out.write(f"Isoforms removed         : {n_removed}\\n")
        out.write(f"Unique genes retained    : {n_after}\\n")
        out.write("\\nSelection rule: best bitscore per gene ID;\\n")
        out.write("tie-break: lowest e-value, then highest HMM coverage.\\n")

    print(f"HMMER: {n_before} hits → {n_after} unique genes "
          f"({n_removed} isoforms removed)", file=sys.stderr)
    """
}
