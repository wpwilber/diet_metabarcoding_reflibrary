target = "ITS1"
target_aliases = 'ITS1 OR "internal transcribed spacer 1" OR "ITS-1"'
primer_f = "GATATCCGTTGCCGAGAGTC"
primer_r = "CCGAAGGCGTCAAGGAACAC"
max_len = "90"
min_len = "50"

rule all:
    input: 
        f"taxonomy/{target}_ref_tax.qza",
        f"sequences/{target}_ref_seq.qza",
        f"sequences/{target}_ref_seq_derep.qza",
        f"taxonomy/{target}_ref_tax_derep.qza",
        f"sequences/{target}_seq_segments.qza",

rule get_ncbi:
    conda: "/users/wwilber/.conda/envs/qiime2-amplicon-2026.1"
    output:
        sequences = f"sequences/{target}_ref_seq.qza",
        taxonomy = f"taxonomy/{target}_ref_tax.qza",
    threads : 1
    shell:
        r"""
        mkdir -p taxonomy
        mkdir -p sequences
        qiime rescript get-ncbi-data \
        --p-query 'txid35493[ORGN] AND ({target_aliases}) NOT environmental sample[Filter] NOT environmental samples[Filter] NOT environmental[Title] NOT uncultured[Title] NOT unclassified[Title] NOT unidentified[Title] NOT unverified[Title]' \
            --p-ranks  kingdom phylum class order family genus species \
            --p-rank-propagation \
            --p-n-jobs 1 \
            --o-sequences {output.sequences} \
            --o-taxonomy  {output.taxonomy} \
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
    threads : 8
    shell:
        r"""
        qiime rescript dereplicate \
            --i-sequences {input.sequences} \
            --i-taxa {input.taxonomy} \
            --p-mode 'uniq' \
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
    threads : 8
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
