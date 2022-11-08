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
    if rule not in config.keys() or config[rule] is None: ## check the rule exists in the config file
        return ""

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
        auspice_json = f"results/{config['virus']}.json",

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
    params:
        allParams = format_config_params('filter')
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --include {input.included_strains} \
            --exclude {input.excluded_strains} \
            --output-log {output.log} \
            --output {output.sequences} \
            {params.allParams}
            """

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
        allParams = format_config_params('align')
    shell:
        """
        augur align \
            --reference-sequence {input.reference} \
            --sequences {input.sequences} \
            --output {output.alignment} \
            {params.debug} \
            {params.allParams}
            """

rule mask:
    message:
        """
        Mask bases in alignment by converting `-` to `N`
        """
    input:
        alignment = rules.align.output.alignment
    output:
        alignment = "results/masked_alignment.fasta"
    params:
        allParams = format_config_params('mask')
    shell:
        """
        python3 helper_scripts/mask-alignment.py \
            --alignment {input.alignment} \
            --output {output.alignment} \
            {params.allParams}
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.mask.output.alignment
    output:
        tree = "results/tree_raw.nwk"
    params:
        allParams = format_config_params('tree')
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            {params.allParams}
            """
rule refine:
    message:
        """
        Refining tree
          - refining branch lengths
          - [if selected] estimating treetime 
        """
    input:
        alignment = rules.mask.output,
        metadata = input_metadata,
        tree = rules.tree.output.tree,
    output:
        branch_lengths = "results/branch_lengths.json",
        tree = "results/tree.nwk",
    params:
        allParams = format_config_params('refine')
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.branch_lengths} \
            {params.allParams}
        """ 

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.mask.output
    output:
        nt_muts = "results/nt_muts.json"
    params:
        allParams = format_config_params('ancestral')
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.nt_muts} \
            {params.allParams}
        """ 

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        nt_muts = rules.ancestral.output.nt_muts,
        reference = reference
    output:
        aa_muts = "results/aa_muts.json"
    # params:
    #     allParams = format_config_params('translate')
    shell:
        """
        augur translate \
            --ancestral-sequences {input.nt_muts} \
            --output-node-data {output.aa_muts} \
            --reference-sequence {input.reference} \
            --tree {input.tree} \
        """
        # {params.allParams}

rule traits:
    message: "Inferring ancestral traits"
    input:
        tree = rules.refine.output.tree,
        metadata = input_metadata
    output:
        traits = "results/traits.json",
    params:
        allParams = format_config_params('traits')
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output-node-data {output.traits} \
            {params.allParams}
        """

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = input_metadata,
        branch_lengths = rules.refine.output.branch_lengths,
        traits = rules.traits.output.traits,
        nt_muts = rules.ancestral.output.nt_muts,
        aa_muts = rules.translate.output.aa_muts,
    output:
        auspice_json = rules.all.input.auspice_json,
    params:
        allParams = format_config_params('export')
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} \
            --output {output.auspice_json} \
            {params.allParams}
        """ 

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
