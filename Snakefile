configfile: "./config.yaml"


# HELPER FUNCTIONS 
nullish_values = [None, 'null', '', 'None', False]

flag_alternatives = {
    "refine": {
        "covariance": "no_covariance",
    },
    "ancestral": {
        "infer_ambiguous": "keep_ambiguous"
    }
}

def format_param(rule, param): 
    value = config[rule][param]

    if value not in nullish_values:
        return "--%s %s"%(param.replace("_", "-"), value)
    else:
        return None

def format_flag(rule, param): 
    value = config[rule][param]
    param = param.split('_FLAG')[0]

    if value not in nullish_values: # any truth-y value -> include flag
        return "--%s"%param.replace("_", "-")

    elif rule in flag_alternatives.keys() and param in flag_alternatives[rule].keys(): # a couple of args have 'inverse' flags 
        alternative = flag_alternatives[rule][param]
        return "--%s"%alternative.replace("_", "-")
    
    else: # false-y value -> don't include flag
        return None

def format_config_params(rule): 
    params = []

    if rule not in config.keys(): ## check the rule exists in the config file
        raise Exception("Rule %s parameters not found in config file" % rule)

    for param in config[rule].keys():
        if param.endswith('_FLAG'): ## format flags 
            formattedFlag = format_flag(rule, param)
            if formattedFlag is not None:
                params.append(formattedFlag)

        else: ## format other parameters
            formattedParam = format_param(rule, param)
            if formattedParam is not None:
                params.append(formattedParam)

    return " ".join(params)

# MAIN WORKFLOW
rule all:
    input:
        auspice_json = "results/{config.virus}.json",

input_fasta = config['data']['sequences'],
input_metadata = config['data']['metadata'],
excluded_strains = config['data']['excluded_strains'],
included_strains = config['data']['included_strains'],
reference = config['data']['reference'],

rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering.
        """
    input:
        sequences = input_fasta
    output:
        sequence_index = "results/sequence_index.tsv"
    params:
        verbose = '--verbose' if config['verbose'] else '',
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index}
            {params.verbose}
        """

rule filter:
    message:
        """
        Filtering strains
        """
    input:
        excluded_strains = excluded_strains,
        included_strains = included_strains,        
        metadata = input_metadata,
        sequences = input_fasta,
        sequence_index = "results/sequence_index.tsv",
    output:
        log = "logs/filter.log",
        sequences = "results/filtered.fasta"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --include {input.included_strains} \
            --exclude {input.excluded_strains} \
            --output-log {output.log}
            --output {output.sequences} \
            """ + format_shell('filter')

rule align:
    message:
        """
        Aligning sequences
        """
    input:
        reference = reference,
        sequences = rules.filter.output.sequences,
    output:
        alignment = "results/aligned.fasta"
    params:
        debug = '--debug' if config['verbose'] else '',
    shell:
        """
        augur align \
            {params.debug} \
            --reference-sequence {input.reference} \
            --sequences {input.sequences} \
            --output {output.alignment} \
        """ + format_shell('align')

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree_raw.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
        """ + format_shell('tree')

rule refine:
    message:
        """
        Refining tree
          - refining branch lengths
          - [if selected] estimating treetime 
        """
    input:
        alignment = rules.align.output,
        metadata = input_metadata,
        tree = rules.tree.output.tree,
    output:
        branch_lengths = "results/branch_lengths.json",
        tree = "results/tree.nwk",
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.branch_lengths} \
            """ + format_shell('refine')

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "results/nt_muts.json"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
        """ + format_shell('ancestral')

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = reference
    output:
        node_data = "results/aa_muts.json"
    shell:
        """
        augur translate \
            --ancestral-sequences {input.node_data} \
            --output-node-data {output.node_data} \
            --reference-sequence {input.reference} \
            --tree {input.tree} \
        """

rule traits:
    message: "Inferring ancestral traits"
    input:
        tree = rules.refine.output.tree,
        metadata = input_metadata
    output:
        node_data = "results/traits.json",
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output-node-data {output.node_data} \
        """ + format_shell('traits')

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = input_metadata,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
    output:
        auspice_json = rules.all.input.auspice_json,
    shell:
        """
        augur export v2 \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} \
            --output {output.auspice_json}
            --tree {input.tree} \
        """ + format_shell('export')

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
