target = "ITS1"
target_aliases = 'ITS1 OR "internal transcribed spacer 1" OR "ITS-1"'
primer_f = "GATATCCGTTGCCGAGAGTC"
primer_r = "CCGAAGGCGTCAAGGAACAC"
max_len = "90"
min_len = "50"

# Iterations: 0, 1, 2
ITER = range(0, 3)

############################################
# FINAL TARGET
############################################

rule all:
    input:
        # final culled outputs for all iterations
        expand(f"sequences/{target}_seq_segments_derep_cull_{{i}}.qza", i=ITER),

        # extracted segments for iterations 1 and 2
        expand(f"sequences/{target}_extracted_seq_segments_{{i}}.qza", i=[1,2]),

        # visualization for all iterations
        expand(f"visualizations/{{target}}_seq_segments_derep_cull_{{i}}.qzv", target=target, i=ITER)

############################################
# ONE-OFF PREPROCESSING STEPS
############################################

rule get_ncbi:
    conda: "/users/wwilber/.conda/envs/qiime2-amplicon-2026.1"
    output:
        sequences = f"sequences/{target}_ref_seq.qza",
        taxonomy = f"taxonomy/{target}_ref_tax.qza",
    threads: 1
    shell:
        r"""
        mkdir -p taxonomy sequences
        qiime rescript get-ncbi-data \
            --p-query 'txid35493[ORGN] AND ({target_aliases}) NOT environmental sample[Filter] NOT environmental samples[Filter] NOT environmental[Title] NOT uncultured[Title] NOT unclassified[Title] NOT unidentified[Title] NOT unverified[Title]' \
            --p-ranks kingdom phylum class order family genus species \
            --p-rank-propagation \
            --p-n-jobs 1 \
            --o-sequences {output.sequences} \
            --o-taxonomy {output.taxonomy} \
            --verbose
        """

rule rescript_dereplicate:
    conda: "/users/wwilber/.conda/envs/qiime2-amplicon-2026.1"
    input:
        sequences = f"sequences/{target}_ref_seq.qza",
        taxonomy = f"taxonomy/{target}_ref_tax.qza",
    output:
        dereplicated_sequences = f"sequences/{target}_ref_seq_derep.qza",
        dereplicated_taxonomy = f"taxonomy/{target}_ref_tax_derep.qza",
    threads: 8
    shell:
        r"""
        qiime rescript dereplicate \
            --i-sequences {input.sequences} \
            --i-taxa {input.taxonomy} \
            --p-mode uniq \
            --p-threads {threads} \
            --o-dereplicated-sequences {output.dereplicated_sequences} \
            --o-dereplicated-taxa {output.dereplicated_taxonomy}
        """

rule extract_reads:
    conda: "/users/wwilber/.conda/envs/qiime2-amplicon-2026.1"
    input:
        sequences = f"sequences/{target}_ref_seq_derep.qza"
    output:
        extracted_reads = f"sequences/{target}_seq_segments.qza"
    threads: 8
    shell:
        r"""
        qiime feature-classifier extract-reads \
            --i-sequences {input.sequences} \
            --p-f-primer {primer_f} \
            --p-r-primer {primer_r} \
            --p-min-length {min_len} \
            --p-max-length {max_len} \
            --p-n-jobs {threads} \
            --o-reads {output.extracted_reads}
        """

############################################
# ITERATIVE REFINEMENT STEPS (i = 0,1,2)
############################################

# extract-seq-segments only runs for i = 1 and i = 2
rule extract_segments:
    conda: "/users/wwilber/.conda/envs/qiime2-amplicon-2026.1"
    wildcard_constraints:
        i = "1|2"
    input:
        ref = f"sequences/{target}_ref_seq_derep.qza",
        refseg = lambda wc: (
            f"sequences/{target}_seq_segments_derep_cull_0.qza"
            if int(wc.i) == 1
            else f"sequences/{target}_seq_segments_derep_cull_{int(wc.i)-1}.qza"
        )
    output:
        seg = f"sequences/{target}_extracted_seq_segments_{{i}}.qza",
        unmatched = f"sequences/{target}_unmatched_sequences_{{i}}.qza"
    threads: 8
    shell:
        r"""
        qiime rescript extract-seq-segments \
            --i-input-sequences {input.ref} \
            --i-reference-segment-sequences {input.refseg} \
            --p-perc-identity 0.7 \
            --p-min-seq-len 10 \
            --p-threads {threads} \
            --o-extracted-sequence-segments {output.seg} \
            --o-unmatched-sequences {output.unmatched}
        """

# dereplicate runs for i = 0,1,2
rule dereplicate_iter:
    conda: "/users/wwilber/.conda/envs/qiime2-amplicon-2026.1"
    wildcard_constraints:
        i = "0|1|2"
    input:
        seq = lambda wc: (
            f"sequences/{target}_seq_segments.qza"
            if int(wc.i) == 0
            else f"sequences/{target}_extracted_seq_segments_{wc.i}.qza"
        ),
        tax = f"taxonomy/{target}_ref_tax_derep.qza"
    output:
        seq_d = f"sequences/{target}_seq_segments_derep_{{i}}.qza",
        tax_d = f"taxonomy/{target}_tax_segments_derep_{{i}}.qza"
    threads: 8
    shell:
        r"""
        qiime rescript dereplicate \
            --i-sequences {input.seq} \
            --i-taxa {input.tax} \
            --p-mode uniq \
            --p-threads {threads} \
            --o-dereplicated-sequences {output.seq_d} \
            --o-dereplicated-taxa {output.tax_d}
        """

# cull runs for i = 0,1,2
rule cull_iter:
    conda: "/users/wwilber/.conda/envs/qiime2-amplicon-2026.1"
    wildcard_constraints:
        i = "0|1|2"
    input:
        seq_d = f"sequences/{target}_seq_segments_derep_{{i}}.qza"
    output:
        clean = f"sequences/{target}_seq_segments_derep_cull_{{i}}.qza"
    threads: 8
    shell:
        r"""
        qiime rescript cull-seqs \
            --i-sequences {input.seq_d} \
            --p-n-jobs {threads} \
            --p-num-degenerates 1 \
            --p-homopolymer-length 8 \
            --o-clean-sequences {output.clean}
        """

# tabulate runs for i = 0,1,2
rule tabulate_iter:
    conda: "/users/wwilber/.conda/envs/qiime2-amplicon-2026.1"
    wildcard_constraints:
        i = "0|1|2"
    input:
        clean = f"sequences/{target}_seq_segments_derep_cull_{{i}}.qza"
    output:
        viz = f"visualizations/{target}_seq_segments_derep_cull_{{i}}.qzv"
    shell:
        r"""
        mkdir -p visualizations
        qiime feature-table tabulate-seqs \
            --i-data {input.clean} \
            --o-visualization {output.viz}
        """

