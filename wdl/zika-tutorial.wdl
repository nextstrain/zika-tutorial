version 1.0

workflow ZikaTutorial {
    input {
        File input_fasta
        File input_metadata
        File exclude
        File reference
        File colors
        File lat_longs
        File auspice_config
    }

    call IndexSequences {
        input:
            input_fasta = input_fasta
    }
    call Filter {
        input:
            input_fasta = input_fasta,
            sequence_index = IndexSequences.sequence_index,
            input_metadata = input_metadata,
            exclude = exclude,
    }
    call Align {
        input:
            filtered_sequences = Filter.filtered_sequences,
            reference = reference,
    }
    call Tree {
        input:
            alignment = Align.alignment,
    }
    call Refine {
        input:
            tree = Tree.tree,
            alignment = Align.alignment,
            input_metadata = input_metadata,
    }
    call Ancestral {
        input:
            time_tree = Refine.time_tree,
            alignment = Align.alignment,
    }
    call Translate {
        input:
            time_tree = Refine.time_tree,
            nt_muts = Ancestral.nt_muts,
            reference = reference,
    }
    call Traits {
        input:
            time_tree = Refine.time_tree,
            input_metadata = input_metadata,
    }
    call Export {
        input:
            time_tree = Refine.time_tree,
            input_metadata = input_metadata,
            branch_lengths = Refine.branch_lengths,
            traits = Traits.traits,
            nt_muts = Ancestral.nt_muts,
            aa_muts = Translate.aa_muts,
            colors = colors,
            lat_longs = lat_longs,
            auspice_config = auspice_config,
    }

    output {
        File auspice_json = Export.auspice_json
    }
}

task IndexSequences {
    input {
        File input_fasta
    }
    command {
        augur index \
            --sequences ~{input_fasta} \
            --output sequence_index.tsv
    }
    output {
        File sequence_index = "sequence_index.tsv"
    }
    runtime {
        docker: "nextstrain/base:latest"
    }
}

task Filter {
    input {
        File input_fasta
        File sequence_index
        File input_metadata
        File exclude

        String? group_by
        Int? sequences_per_group
        String? min_date
    }
    command {
        augur filter \
            --sequences ~{input_fasta} \
            --sequence-index ~{sequence_index} \
            --metadata ~{input_metadata} \
            --exclude ~{exclude} \
            ~{"--group-by " + group_by} \
            ~{"--sequences-per-group " + sequences_per_group} \
            ~{"--min-date " + min_date} \
            --output filtered.fasta
    }
    output {
        File filtered_sequences = "filtered.fasta"
    }
    runtime {
        docker: "nextstrain/base:latest"
    }
}

task Align {
    input {
        File filtered_sequences
        File reference
    }
    command {
        augur align \
            --sequences ~{filtered_sequences} \
            --reference-sequence ~{reference} \
            --output aligned.fasta \
            --fill-gaps
    }
    output {
        File alignment = "aligned.fasta"
    }
    runtime {
        docker: "nextstrain/base:latest"
    }
}

task Tree {
    input {
        File alignment
    }
    command {
        augur tree \
            --alignment ~{alignment} \
            --output tree_raw.nwk
    }
    output {
        File tree = "tree_raw.nwk"
    }
    runtime {
        docker: "nextstrain/base:latest"
    }
}

task Refine {
    input {
        File tree
        File alignment
        File input_metadata

        String coalescent = "opt"
        String date_inference = "marginal"
        Int clock_filter_iqd = 4
    }
    command {
        augur refine \
            --tree ~{tree} \
            --alignment ~{alignment} \
            --metadata ~{input_metadata} \
            --timetree \
            --coalescent ~{coalescent} \
            --date-confidence \
            --date-inference ~{date_inference} \
            --clock-filter-iqd ~{clock_filter_iqd} \
            --output-tree tree.nwk \
            --output-node-data branch_lengths.json
    }
    output {
        File time_tree = "tree.nwk"
        File branch_lengths = "branch_lengths.json"
    }
    runtime {
        docker: "nextstrain/base:latest"
    }
}

task Ancestral {
    input {
        File time_tree
        File alignment

        String inference = "joint"
    }
    command {
        augur ancestral \
            --tree ~{time_tree} \
            --alignment ~{alignment} \
            --inference ~{inference} \
            --output-node-data nt_muts.json
    }
    output {
        File nt_muts = "nt_muts.json"
    }
    runtime {
        docker: "nextstrain/base:latest"
    }
}

task Translate {
    input {
        File time_tree
        File nt_muts
        File reference
    }
    command {
        augur translate \
            --tree ~{time_tree} \
            --ancestral-sequences ~{nt_muts} \
            --reference-sequence ~{reference} \
            --output-node-data aa_muts.json
    }
    output {
        File aa_muts = "aa_muts.json"
    }
    runtime {
        docker: "nextstrain/base:latest"
    }
}

task Traits {
    input {
        File time_tree
        File input_metadata

        String columns = "region country"
    }
    command {
        augur traits \
            --tree ~{time_tree} \
            --metadata ~{input_metadata} \
            --columns ~{columns} \
            --confidence \
            --output-node-data traits.json
    }
    output {
        File traits = "traits.json"
    }
    runtime {
        docker: "nextstrain/base:latest"
    }
}

task Export {
    input {
        File time_tree
        File input_metadata
        File branch_lengths
        File traits
        File nt_muts
        File aa_muts
        File colors
        File lat_longs
        File auspice_config
    }
    command {
        augur export v2 \
            --tree ~{time_tree} \
            --metadata ~{input_metadata} \
            --node-data ~{branch_lengths} ~{traits} ~{nt_muts} ~{aa_muts} \
            --colors ~{colors} \
            --lat-longs ~{lat_longs} \
            --auspice-config ~{auspice_config} \
            --output zika.json
    }
    output {
        File auspice_json = "zika.json"
    }
    runtime {
        docker: "nextstrain/base:latest"
    }
}
