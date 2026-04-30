// modules/interproscan.nf
// InterProScan: domain-based confirmation of HMMER candidates

process INTERPROSCAN {
    tag "interproscan"
    label 'high_cpu'
    publishDir "${params.outdir}/02_interproscan", mode: 'copy'

    input:
    path candidates_fasta

    output:
    path "interproscan_results.tsv",  emit: tsv
    path "interproscan_results.xml",  emit: xml
    path "interproscan_results.gff3", emit: gff3

    script:
    """
    # InterProScan does not accept stop-codon asterisks
    sed 's/\*//g' ${candidates_fasta} > candidates_clean.fasta

    interproscan.sh \
        --input        candidates_clean.fasta \
        --output-dir   . \
        --formats      TSV,XML,GFF3 \
        --cpu          ${task.cpus} \
        --applications Pfam,PANTHER,Gene3D,SUPERFAMILY,PRINTS,ProSiteProfiles \
        --goterms \
        --iprlookup \
        --pathways \
        ${params.interpro_db ? "--data-dir ${params.interpro_db}" : ""} \
        --outfile-base interproscan_results

    # Normalize output filenames (IPS naming varies by version)
    for ext in tsv xml gff3; do
        if ! [ -f "interproscan_results.\${ext}" ]; then
            mv interproscan_results*."\${ext}" "interproscan_results.\${ext}" 2>/dev/null || true
        fi
    done
    """
}

process CONFIRM_CANDIDATES {
    tag "confirm_domains"
    label 'low_cpu'
    publishDir "${params.outdir}/02_interproscan", mode: 'copy'

    input:
    path iprscan_tsv
    path hmmer_ids
    val  target_domains

    output:
    path "confirmed_ids.txt",               emit: confirmed_ids
    path "rejected_ids.txt",                emit: rejected_ids
    path "domain_confirmation_summary.tsv", emit: summary

    script:
    """
    #!/usr/bin/env python3
    import sys
    from collections import defaultdict

    target_set = set(d.strip() for d in "${target_domains}".split(",") if d.strip())

    with open("${hmmer_ids}") as fh:
        hmmer_ids = set(line.strip() for line in fh if line.strip())

    # InterProScan TSV columns:
    # 0:seq_id  1:md5  2:length  3:analysis  4:acc  5:name
    # 6:start   7:stop 8:score   9:status   10:date
    # 11:ipr_acc 12:ipr_desc 13:go_terms 14:pathways
    seq_domains = defaultdict(set)
    with open("${iprscan_tsv}") as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            seq_id   = cols[0]
            pfam_acc = cols[4]
            ipr_acc  = cols[11] if len(cols) > 11 else ""
            if pfam_acc:
                seq_domains[seq_id].add(pfam_acc)
            if ipr_acc:
                seq_domains[seq_id].add(ipr_acc)

    confirmed = []
    rejected  = []
    rows      = []

    for sid in sorted(hmmer_ids):
        found      = seq_domains.get(sid, set()) & target_set
        has_domain = len(found) > 0
        status     = "CONFIRMED" if has_domain else "REJECTED"
        (confirmed if has_domain else rejected).append(sid)
        rows.append((sid, ";".join(sorted(found)) or "none", status))

    with open("confirmed_ids.txt", "w") as out:
        out.write("\n".join(confirmed) + ("\n" if confirmed else ""))

    with open("rejected_ids.txt", "w") as out:
        out.write("\n".join(rejected) + ("\n" if rejected else ""))

    with open("domain_confirmation_summary.tsv", "w") as out:
        out.write("seq_id\tmatched_domains\tstatus\n")
        for r in rows:
            out.write("\t".join(r) + "\n")

    print(f"InterProScan: {len(confirmed)} confirmed, {len(rejected)} rejected", file=sys.stderr)
    """
}
