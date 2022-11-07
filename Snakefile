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
        return "--%s"%param.replace("_", "-")
    else:
        return None

def format_flag(rule, param): 
    value = config[rule][param]
    param = param.split('_FLAG')[0]
    if value not in nullish_values:
        return "--%s"%param.replace("_", "-")
    elif rule in flag_alternatives.keys() and param in flag_alternatives[rule].keys():
        alternative = flag_alternatives[rule][param]
        return "--%s"%alternative.replace("_", "-")
    else:
        return None

def format_shell(rule): 
    params = []
    if rule not in config.keys():
        raise Exception("Rule %s parameters not found in config file" % rule)
    for param in config[rule].keys():
        if param.endswith('_FLAG'):
            formattedFlag = params.append(format_flag(rule, param))
            if formattedFlag:
                params.append(formattedFlag)
        else:
            formattedParam = params.append(format_param(rule, param))
            if formattedParam:
                params.append(formattedParam)
    returnValue = " ".join(params)
    print(returnValue)
    return returnValue

# MAIN WORKFLOW

rule all:
    input:
        auspice_json = "results/{config.virus}.json",

input_fasta = config['data']['sequences'],
input_metadata = config['data']['metadata'],
excluded_strains = config['data']['excluded_strains'],
included_strains = config['data']['included_strains'],
reference = config['data']['reference'],
colors = config['export']['colors'],
lat_longs = config['export']['lat_longs'],
auspice_config = config['export']['auspice_config']

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
        exclude_where = config['filter']['exclude_where'],
        exclude_all = format_flag("filter", "exclude_all"),
        exclude_ambiguous_dates_by = config['filter']['exclude_ambiguous_dates_by'],
        include_where = config['filter']['include_where'],
        max_date = config['filter']['max_date'],
        metadata_chunk_size = config['filter']['metadata_chunk_size'],
        metadata_id_columns = config['filter']['metadata_id_columns'],
        min_date = config['filter']['min_date'],
        min_length = config['filter']['min_length'],
        non_nucleotide = format_flag("filter", "non_nucleotide"),
    shell:
        """
        echo {params} &&  \       
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --include {input.included_strains} \
            --include-where {params.include_where} \
            --exclude {input.excluded_strains} \
            {params.exclude_all} \
            --exclude-ambiguous-dates-by {params.exclude_ambiguous_dates_by} \
            --exclude-where {params.exclude_where} \
            --max-date {params.max_date} \
            --metadata-chunk-size {params.metadata_chunk_size} \
            --metadata-id-columns {params.metadata_id_columns} \
            --min-date {params.min_date} \
            --min-length {params.min_length} \
            {params.non_nucleotide} \
            --output-log {output.log}
            --output {output.sequences} \
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
        existing_alignment = config['align']['existing_alignment'],
        debug = '--debug' if config['verbose'] else '',
        fill_gaps = format_flag('align', 'fill_gaps'),
        method = config['align']['method'],
        nthreads = config['nthreads'],
        reference_name = config['data']['reference_name'],
        remove_reference = format_flag('align', 'remove_reference'),
    shell:
        """
        augur align \
            {params.debug} \
            --existing-alignment {params.existing_alignment} \
            {params.fill_gaps} \
            --method {params.method} \
            --nthreads {params.nthreads}
            --output {output.alignment} \
            --reference-sequence {input.reference} \
            --reference-name {params.reference_name} \
            {params.remove_reference} \
            --sequences {input.sequences} \
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree_raw.nwk"
    params:
        exclude_sites = config['tree']['exclude_sites'],
        method = config['tree']['method'],
        nthreads = config['nthreads'],
        override_default_args = config['tree']['override_default_args'],
        substitution_model = config['tree']['substitution_model'],
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --exclude-sites {params.exclude_sites} \
            --method {params.method} \
            --nthreads {params.nthreads} \
            --output {output.tree}
            --override-default-args {params.override_default_args} \
            --substitution-model {params.substitution_model} \
        """

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
        node_data = "results/branch_lengths.json",
        tree = "results/tree.nwk",
    params:
        branch_length_inference = config['refine']['branch_length_inference'],
        clock_filter_iqd = config['refine']['clock_filter_iqd'],
        clock_rate = config['refine']['clock_rate'],
        clock_std_dev = config['refine']['clock_std_dev'],
        coalescent = config['refine']['coalescent'],
        covariance = format_flag('refine', 'no_covariance', '--covariance'),
        date_confidence = format_flag('refine', 'date_confidence'),
        date_format = config['refine']['date_format'],
        date_inference = config['refine']['date_inference'],
        divergence_units = config['refine']['divergence_units'],
        gen_per_year = config['refine']['gen_per_year'],
        keep_polytomies = format_flag('refine', 'keep_polytomies'),
        keep_root = format_flag('refine', 'keep_root'),
        precision = config['refine']['precision'],
        root = config['refine']['root'],
        seed = config['refine']['seed'],
        timetree = format_flag('refine', 'timetree'),
        use_fft = format_flag('refine', 'use_fft'),
        year_bounds = config['refine']['year_bounds'],
        verbosity = 6 if config['verbose'] else 0
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --branch_length_inference {params.branch_length_inference} \
            --clock-filter-iqd {params.clock_filter_iqd} \
            --clock-rate {params.clock_rate} \
            --clock-std-dev {params.clock_std_dev} \
            --coalescent {params.coalescent} \
            {params.covariance} \
            {params.date_confidence} \
            --date-format {params.date_format} \
            --date-inference {params.date_inference} \
            --divergence-units {params.divergence_units} \
            --gen-per-year {params.gen_per_year} \
            {params.keep_polytomies} \
            {params.keep_root} \
            --precision {params.precision} \
            --root {params.root} \
            --seed {params.seed} \
            {params.timetree} \
            {params.use_fft} \
            --year-bounds {params.year_bounds} \
            --verbosity {params.verbosity} \
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "results/nt_muts.json"
    params:
        inference = config['ancestral']['inference'],
        infer_ambiguous = format_flag('ancestral', 'infer_ambiguous', '--keep-ambiguous'),
        keep_overhangs = format_flag('ancestral', 'keep_overhangs'),

    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
            {params.infer_ambiguous} \
            {params.keep_overhangs} \
        """

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
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = input_metadata
    output:
        node_data = "results/traits.json",
    params:
        columns = config['traits']['columns'],
        confidence = format_flag('traits', 'confidence'),
        sampling_bias_correction = config['traits']['sampling_bias_correction'],
        weights = config['traits']['weights'],

    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output-node-data {output.node_data} \
            --columns {params.columns} \
            {params.confidence} \
            --sampling-bias-correction {params.sampling_bias_correction} \
            --weights {params.weights} \
        """

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
    params: 
        auspice_config = config['export']['auspice_config'],
        build_url = config['export']['build_url'],
        colors = config['export']['colors'],
        color_by_metadata = config['export']['color_by_metadata'],
        description = config['export']['description'],
        geo_resolutions = config['export']['geo_resolutions'],
        include_root_sequence = config['export']['include_root_sequence'],
        lat_longs = config['export']['lat_longs'],
        maintainers = config['export']['maintainers'],
        minify_json = format_flag('export', 'minify_json'),
        panels = config['export']['panels'],
        skip_validation = format_flag('export', 'skip_validation'),
        title = config['export']['title'],
    shell:
        """
        augur export v2 \
            --auspice-config {params.auspice_config} \
            --build-url {params.build_url} \
            --colors {params.colors} \
            --color-by-metadata {params.color_by_metadata} \
            --description {params.description} \
            --geo-resolutions {params.geo_resolutions} \
            --include-root-sequence \
            --lat-longs {params.lat_longs} \
            --maintainers {params.maintainers} \
            --metadata {input.metadata} \
            {params.minify_json} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} \
            --output {output.auspice_json}
            --panels {params.panels} \
            {params.skip_validation} \
            --title {params.title} \
            --tree {input.tree} \
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
