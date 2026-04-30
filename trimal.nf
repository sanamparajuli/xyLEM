// modules/trimal.nf
// TrimAl: alignment trimming

process TRIMAL_TRIM {
    tag "trimal_trim"
    label 'low_cpu'
    publishDir "${params.outdir}/05_trimming", mode: 'copy'

    input:
    path alignment

    output:
    path "trimmed_alignment.fasta",  emit: trimmed_aln
    path "trimmed_alignment.html",   emit: html_report
    path "trimal_stats.tsv",         emit: stats

    script:
    def trimal_flags = params.trimal_method == "automated1" ? "-automated1"
                     : params.trimal_method == "gappyout"   ? "-gappyout"
                     : params.trimal_method == "strict"     ? "-strict"
                     : params.trimal_method == "strictplus" ? "-strictplus"
                     : "-automated1"
    """
    trimal \
        -in      ${alignment} \
        -out     trimmed_alignment.fasta \
        -htmlout trimmed_alignment.html \
        ${trimal_flags} \
        -fasta

    python3 - <<'EOF'
def parse_fasta(path):
    seqs = {}
    current = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                current = line[1:].split()[0]
                seqs[current] = ""
            elif current:
                seqs[current] += line
    return seqs

pre  = parse_fasta("${alignment}")
post = parse_fasta("trimmed_alignment.fasta")

pre_len  = len(next(iter(pre.values()),  ""))
post_len = len(next(iter(post.values()), ""))
retained = round(post_len / pre_len * 100, 1) if pre_len > 0 else 0.0

with open("trimal_stats.tsv", "w") as out:
    out.write("metric\tvalue\n")
    out.write(f"sequences_input\t{len(pre)}\n")
    out.write(f"sequences_output\t{len(post)}\n")
    out.write(f"columns_before_trim\t{pre_len}\n")
    out.write(f"columns_after_trim\t{post_len}\n")
    out.write(f"columns_retained_pct\t{retained}\n")
    out.write(f"trim_method\t${params.trimal_method}\n")

print(f"TrimAl: {pre_len} -> {post_len} columns ({retained}% retained)")
EOF
    """
}
