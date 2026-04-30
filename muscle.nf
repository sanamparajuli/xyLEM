// modules/muscle.nf
// MUSCLE v5: multiple sequence alignment (with optional outgroup)

process MUSCLE_ALIGN {
    tag "muscle_align"
    label 'medium_cpu'
    publishDir "${params.outdir}/04_alignment", mode: 'copy'

    input:
    tuple path(candidate_fasta), path(reference_fasta)
    path  outgroup_fasta

    output:
    path "combined_aligned.fasta",  emit: alignment
    path "sequence_manifest.tsv",   emit: manifest
    path "muscle.log",              emit: log

    script:
    """
    #!/usr/bin/env python3
    import os, subprocess, sys
    from collections import OrderedDict

    def read_fasta(filepath):
        seqs = OrderedDict()
        current_id  = None
        current_seq = []
        with open(filepath) as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith(">"):
                    if current_id:
                        seqs[current_id] = "".join(current_seq)
                    current_id  = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
        if current_id:
            seqs[current_id] = "".join(current_seq)
        return seqs

    candidates = read_fasta("${candidate_fasta}")
    references = read_fasta("${reference_fasta}")

    outgroup      = {}
    outgroup_file = "${outgroup_fasta}"
    if os.path.exists(outgroup_file) and os.path.getsize(outgroup_file) > 0:
        outgroup = read_fasta(outgroup_file)
        print(f"Outgroup sequences loaded: {len(outgroup)}", file=sys.stderr)
    else:
        print("No outgroup provided - proceeding without outgroup.", file=sys.stderr)

    merged = OrderedDict()
    source = {}
    for sid, seq in candidates.items():
        if sid not in merged:
            merged[sid] = seq
            source[sid] = "candidate"
    for sid, seq in references.items():
        if sid not in merged:
            merged[sid] = seq
            source[sid] = "reference"
    for sid, seq in outgroup.items():
        if sid not in merged:
            merged[sid] = seq
            source[sid] = "outgroup"

    print(f"Total for alignment: {len(merged)} "
          f"({len(candidates)} candidates, {len(references)} references, "
          f"{len(outgroup)} outgroup)", file=sys.stderr)

    with open("sequence_manifest.tsv", "w") as out:
        out.write("seq_id\tsource\tis_outgroup\n")
        for sid in merged:
            is_og = "true" if source[sid] == "outgroup" else "false"
            out.write(f"{sid}\t{source[sid]}\t{is_og}\n")

    with open("combined_dedup.fasta", "w") as out:
        for sid, seq in merged.items():
            out.write(f">{sid}\n{seq}\n")

    seq_count = len(merged)
    algo_flag = "-super5" if seq_count > 1000 else "-align"

    with open("muscle.log", "w") as log:
        log.write(f"Input sequences  : {seq_count}\n")
        log.write(f"  Candidates     : {len(candidates)}\n")
        log.write(f"  References     : {len(references)}\n")
        log.write(f"  Outgroup       : {len(outgroup)}\n")
        log.write(f"Algorithm        : {algo_flag}\n\n")

    cmd = ["muscle", algo_flag, "combined_dedup.fasta",
           "-output", "combined_aligned.fasta",
           "-threads", "${task.cpus}"]
    with open("muscle.log", "a") as log:
        result = subprocess.run(cmd, stderr=log, stdout=log)
    if result.returncode != 0:
        print("ERROR: MUSCLE failed - check muscle.log", file=sys.stderr)
        sys.exit(result.returncode)

    print("Alignment complete: combined_aligned.fasta", file=sys.stderr)
    """
}
