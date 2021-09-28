# Nextstrain build for Zika virus tutorial

This repository provides the data and scripts associated with the [Zika virus tutorial](https://nextstrain.org/docs/getting-started/zika-tutorial).

See the [original Zika build repository](https://github.com/nextstrain/zika) for more about the public build.

## Run the workflow with Cromwell

Install Java 11 and Auspice.

```
mamba create -n java -c conda-forge -c bioconda openjdk=11 auspice=2.29.1
conda activate java
```

[Download Cromwell and WOMtool binaries](https://github.com/broadinstitute/cromwell/releases).

```
mkdir -p bin
cd bin
curl -OL https://github.com/broadinstitute/cromwell/releases/download/69/cromwell-69.jar
curl -OL https://github.com/broadinstitute/cromwell/releases/download/69/womtool-69.jar
cd ..
```

Start Docker. Run the workflow.

```
java -jar bin/cromwell-69.jar run wdl/zika-tutorial.wdl -i wdl/inputs.json
```

View the resulting Auspice JSON in the execution directory of the Cromwell workflow (replacing `<workflow-id>` with the id for your run).

```
auspice view --datasetDir cromwell-executions/ZikaTutorial/<workflow-id>/call-Export/execution/
```
