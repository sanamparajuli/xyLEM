// blast.nf
// ── BLAST: reciprocal best hit (RBH) validation ───────────────

process MAKEBLASTDB {
    tag "makeblastdb_${db_name}"
    label 'low_cpu'
    publishDir "${params.outdir}/03_blast/${db_name}", mode: 'copy'

    input:
    path fasta
    val  db_name

    output:
    tuple val(db_name), path("${db_name}.*"), emit: db
    path "${db_name}.log",                    emit: log

    script:
    """
    makeblastdb \\
        -in          ${fasta} \\
        -dbtype      prot \\
        -out         ${db_name} \\
        -parse_seqids \\
        -logfile     ${db_name}.log
    """
}

process BLAST_FORWARD {
    tag "blast_forward"
    label 'medium_cpu'
    publishDir "${params.outdir}/03_blast", mode: 'copy'

    input:
    path query_fasta
    tuple val(db_name), path(db_files)

    output:
    path "blast_forward.tsv", emit: results

    script:
    def db_prefix = db_files.find { it.name.endsWith('.phr') }?.name?.replaceAll(/\.phr$/, '') ?: db_name
    """
    blastp \\
        -query           ${query_fasta} \\
        -db              ${db_prefix} \\
        -out             blast_forward.tsv \\
        -outfmt          "6 qseqid sseqid pident length qlen slen qcovs evalue bitscore" \\
        -evalue          ${params.blast_evalue} \\
        -num_threads     ${task.cpus} \\
        -max_target_seqs 5
    """
}

process BLAST_REVERSE {
    tag "blast_reverse"
    label 'medium_cpu'
    publishDir "${params.outdir}/03_blast", mode: 'copy'

    input:
    path reference_fasta
    tuple val(db_name), path(db_files)

    output:
    path "blast_reverse.tsv", emit: results

    script:
    def db_prefix = db_files.find { it.name.endsWith('.phr') }?.name?.replaceAll(/\.phr$/, '') ?: db_name
    """
    blastp \\
        -query           ${reference_fasta} \\
        -db              ${db_prefix} \\
        -out             blast_reverse.tsv \\
        -outfmt          "6 qseqid sseqid pident length qlen slen qcovs evalue bitscore" \\
        -evalue          ${params.blast_evalue} \\
        -num_threads     ${task.cpus} \\
        -max_target_seqs 5
    """
}

process RBH_ORTHOLOGY {
    tag "rbh_orthology"
    label 'low_cpu'
    publishDir "${params.outdir}/03_blast", mode: 'copy'

    input:
    path forward_blast
    path reverse_blast
    path confirmed_ids

    output:
    path "all_confirmed_ids.txt", emit: all_ids
    path "orthology_table.tsv",   emit: orthology
    path "rbh_pairs.tsv",         emit: rbh_pairs
    path "rbh_summary.txt",       emit: summary

    script:
    """
    python3 << 'PYEOF'
forward_blast   = '${forward_blast}'
reverse_blast   = '${reverse_blast}'
confirmed_ids   = '${confirmed_ids}'
evalue_thresh   = float("${params.blast_evalue}")
identity_thresh = float("${params.blast_identity}")
coverage_thresh = float("${params.blast_coverage}")

import sys
from collections import defaultdict

def parse_blast(filepath):
    hits = defaultdict(list)
    with open(filepath) as fh:
        for line in fh:
            cols = line.strip().split("\\t")
            if len(cols) < 9:
                continue
            qid      = cols[0]
            sid      = cols[1]
            pident   = float(cols[2])
            qcovs    = float(cols[6])
            evalue   = float(cols[7])
            bitscore = float(cols[8])
            if evalue   > evalue_thresh:   continue
            if pident   < identity_thresh: continue
            if qcovs    < coverage_thresh: continue
            hits[qid].append({
                'sid': sid, 'evalue': evalue,
                'pident': pident, 'bitscore': bitscore, 'qcovs': qcovs
            })
    best = {}
    for qid, hit_list in hits.items():
        best[qid] = sorted(hit_list, key=lambda x: -x['bitscore'])[0]
    return best

fwd = parse_blast(forward_blast)
rev = parse_blast(reverse_blast)

with open(confirmed_ids) as fh:
    all_ids = [line.strip() for line in fh if line.strip()]

rbh_pairs      = []
orthology_rows = []

for seq_id in all_ids:
    if seq_id in fwd:
        ref_hit  = fwd[seq_id]['sid']
        fwd_eval = fwd[seq_id]['evalue']
        fwd_pid  = fwd[seq_id]['pident']
        fwd_cov  = fwd[seq_id]['qcovs']

        is_rbh = (ref_hit in rev and rev[ref_hit]['sid'] == seq_id)

        if is_rbh:
            orthology_type = "RBH_ORTHOLOG"
            rev_eval = rev[ref_hit]['evalue']
            rev_pid  = rev[ref_hit]['pident']
            rbh_pairs.append((seq_id, ref_hit, fwd_eval, fwd_pid,
                              fwd_cov, rev_eval, rev_pid))
        else:
            orthology_type = "PUTATIVE_HOMOLOG"
            rev_eval = rev.get(ref_hit, {}).get('evalue', 'NA')
            rev_pid  = rev.get(ref_hit, {}).get('pident', 'NA')

        orthology_rows.append((
            seq_id, orthology_type, ref_hit,
            f"{fwd_eval:.2e}", f"{fwd_pid:.1f}", f"{fwd_cov:.1f}",
            str(rev_eval) if isinstance(rev_eval, str) else f"{rev_eval:.2e}",
            str(rev_pid)  if isinstance(rev_pid,  str) else f"{rev_pid:.1f}"
        ))
    else:
        orthology_rows.append((
            seq_id, "NO_BLAST_HIT", "NA", "NA", "NA", "NA", "NA", "NA"
        ))

with open("all_confirmed_ids.txt", "w") as out:
    for sid in all_ids:
        out.write(sid + "\\n")

with open("orthology_table.tsv", "w") as out:
    out.write("seq_id\\torthology_type\\tbest_ref_hit\\t"
              "fwd_evalue\\tfwd_pident\\tfwd_qcovs\\t"
              "rev_evalue\\trev_pident\\n")
    for row in orthology_rows:
        out.write("\\t".join(row) + "\\n")

with open("rbh_pairs.tsv", "w") as out:
    out.write("candidate\\treference_hit\\tfwd_evalue\\tfwd_pident\\t"
              "fwd_qcovs\\trev_evalue\\trev_pident\\n")
    for row in rbh_pairs:
        out.write("\\t".join(str(x) for x in row) + "\\n")

type_counts = defaultdict(int)
for row in orthology_rows:
    type_counts[row[1]] += 1

with open("rbh_summary.txt", "w") as out:
    out.write("RBH Orthology Assignment Summary\\n")
    out.write("=================================\\n")
    out.write(f"Total sequences (all pass downstream) : {len(all_ids)}\\n")
    out.write(f"  RBH orthologs                       : {type_counts['RBH_ORTHOLOG']}\\n")
    out.write(f"  Putative homologs (fwd hit only)    : {type_counts['PUTATIVE_HOMOLOG']}\\n")
    out.write(f"  No significant BLAST hit            : {type_counts['NO_BLAST_HIT']}\\n")
    out.write("\\nNote: ALL sequences proceed to alignment/phylogenetics.\\n")
    out.write("Orthology assignments are annotations only.\\n")

print(f"Orthology: {type_counts['RBH_ORTHOLOG']} RBH, "
      f"{type_counts['PUTATIVE_HOMOLOG']} putative, "
      f"{type_counts['NO_BLAST_HIT']} no-hit -- all {len(all_ids)} proceed",
      file=sys.stderr)
PYEOF
    """
}
