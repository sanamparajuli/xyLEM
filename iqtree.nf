// modules/iqtree.nf
// IQ-TREE3: maximum-likelihood phylogenetic tree
//
// Outgroup note:
//   --outgroup_taxon is entirely optional. If not set, IQ-TREE3 produces an
//   unrooted tree — which is standard practice. Users can root the tree
//   at their outgroup in any downstream tool (FigTree, iTOL, ggtree, ETE3 etc.)
//   using the sequence_manifest.tsv produced by MUSCLE_ALIGN to identify
//   which sequences are outgroup members.

process IQTREE3 {
    tag "iqtree3"
    label 'high_cpu'
    publishDir "${params.outdir}/06_phylogeny", mode: 'copy'

    input:
    path trimmed_alignment

    output:
    path "phylogeny.treefile",          emit: treefile
    path "phylogeny.iqtree",            emit: iqtree_log
    path "phylogeny.log",               emit: run_log
    path "phylogeny.contree",           emit: contree,  optional: true
    path "phylogeny.mldist",            emit: mldist,   optional: true
    path "phylogeny.model.gz",          emit: model,    optional: true
    path "phylogeny_ufboot.splits.nex", emit: splits,   optional: true

    script:
    def model_flag     = params.iqtree_model == "TEST" ? "-m TEST" : "-m ${params.iqtree_model}"
    def bootstrap_flag = params.iqtree_bootstrap > 0   ? "-B ${params.iqtree_bootstrap} --bnni" : ""
    // outgroup_taxon is optional — if not set the tree is left unrooted.
    // Users can root in FigTree / iTOL / ggtree using the outgroup sequences
    // identified in 04_alignment/sequence_manifest.tsv
    def outgroup_flag  = (params.outgroup_taxon != null && params.outgroup_taxon != "")
                         ? "-o ${params.outgroup_taxon}"
                         : ""
    """
    seq_count=\$(grep -c "^>" ${trimmed_alignment})
    echo "Building ML tree with \${seq_count} sequences..." | tee phylogeny.log

    if [ "\${seq_count}" -lt 4 ]; then
        echo "ERROR: need >=4 sequences for IQ-TREE; found \${seq_count}" | tee -a phylogeny.log
        exit 1
    fi

    if [ -n "${outgroup_flag}" ]; then
        echo "Outgroup: ${params.outgroup_taxon}" | tee -a phylogeny.log
    else
        echo "No outgroup_taxon set — tree will be unrooted. Root manually in your preferred viewer." | tee -a phylogeny.log
    fi

    iqtree3 \\
        -s      ${trimmed_alignment} \\
        --prefix phylogeny \\
        ${model_flag} \\
        ${bootstrap_flag} \\
        ${outgroup_flag} \\
        -T      ${task.cpus} \\
        --runs  5 \\
        -alrt   1000 \\
        --redo \\
        2>&1 | tee -a phylogeny.log

    python3 - <<'EOF'
import sys, os
treefile = "phylogeny.treefile"
if not os.path.exists(treefile):
    print("ERROR: treefile not produced", file=sys.stderr)
    sys.exit(1)
with open(treefile) as fh:
    tree_str = fh.read().strip()
leaf_count = tree_str.count(",") + 1
print(f"Tree produced with ~{leaf_count} taxa", file=sys.stderr)
EOF
    """
}
