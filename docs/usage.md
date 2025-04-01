# lcmb: Usage


## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Samplesheet input

The input sample sheet should be in a tab delimited format (extension must be .tsv). The required columns depend on which subworkflows you want to run (column names must be accurate but no need to be in this order, redundant columns will be ignored). Examples of a samplesheet that meets the minimal column requirements for relevant scenarios are as follows

| Scenario | Example |
| -------- | ------- |
| Running conpair only | [here](../assets/samplesheet_conpair.tsv) |
| Running filter SNV only or Running filter SNV + phylogenetics for SNVs | [here](../assets/samplesheet_filter-snv_phylogenetics.tsv) |
| Running filter INDELs only | [here](../assets/samplesheet_filter-indel.tsv) |
| Running filter INDELs + phylogenetics without running filter SNVs | [here](../assets/samplesheet_filter-indel_phylogenetics.tsv) |
| Running filter SNVs + filter INDELS or filter SNVs + filter INDELS + phylogenetics or full pipeline | [here](../assets/samplesheet_full_pipeline.tsv) |
| Running phylogenetics pipeline (no filtering) for SNVs then use topology from SNVs to run phylogenetics pipeline for INDELs | [here](../assets/samplesheet_phylogenetics_snv-then-indel.tsv) |
| Running phylogenetics pipeline (no filtering) with building tree topology (usually used for SNVs) | [here](../assets/samplesheet_phylogenetics_no_topology.tsv) |
| Running phylogenetics pipeline (no filtering) without building tree topology (topology should be provided, usually used for INDELs) | [here](../assets/samplesheet_phylogenetics_provided_topology.tsv) |


Detailed descriptions of the requirements are as follows

##### Conpair pipeline

Conpair pipeline (`params.run_conpair==true`) requires the following columns `sample_id`, `match_normal_id`, `pdid`, `bam`, `bai`, `bam_match `, and `bai_match`. An example of a samplesheet for a conpair workflow only can be found [here](../assets/samplesheet_conpair.tsv)

##### Filter SNVs

If Filter SNVs is run (`params.run_filter_snv==true`), samplesheet must include the additional columns `snv_vcf`, `snv_vcf_tbi`, `bas`, and `met`. It should look something like this:

| sample_id       | match_normal_id | pdid    | bam_match                    | bai_match                        | snv_vcf                                                | snv_vcf_tbi                                                | bam                                 | bai                                     | bas                                     | met                                        |
| --------------- | --------------- | ------- | ---------------------------- | -------------------------------- | ------------------------------------------------------ | ---------------------------------------------------------- | ----------------------------------- | --------------------------------------- | --------------------------------------- | ------------------------------------------ |
| PD47151n_lo0002 | PD47151b        | PD47151 | path/to/PD47151/PD47151b.bam | path/to/PD47151/PD47151b.bam.bai | path/to/PD47151/PD47151n_lo0002.caveman_c.annot.vcf.gz | path/to/PD47151/PD47151n_lo0002.caveman_c.annot.vcf.gz.tbi | path/to/PD47151/PD47151n_lo0002.bam | path/to/PD47151/PD47151n_lo0002.bam.bai | path/to/PD47151/PD47151n_lo0002.bam.bas | path/to/PD47151/PD47151n_lo0002.bam.met.gz |
| PD47151n_lo0004 | PD47151b        | PD47151 | path/to/PD47151/PD47151b.bam | path/to/PD47151/PD47151b.bam.bai | path/to/PD47151/PD47151n_lo0004.caveman_c.annot.vcf.gz | path/to/PD47151/PD47151n_lo0004.caveman_c.annot.vcf.gz.tbi | path/to/PD47151/PD47151n_lo0004.bam | path/to/PD47151/PD47151n_lo0004.bam.bai | path/to/PD47151/PD47151n_lo0004.bam.bas | path/to/PD47151/PD47151n_lo0004.bam.met.gz |
| PD52103n_lo0002 | PD52103b        | PD52103 | path/to/PD52103/PD52103b.bam | path/to/PD52103/PD52103b.bam.bai | path/to/PD52103/PD52103n_lo0002.caveman_c.annot.vcf.gz | path/to/PD52103/PD52103n_lo0002.caveman_c.annot.vcf.gz.tbi | path/to/PD52103/PD52103n_lo0002.bam | path/to/PD52103/PD52103n_lo0002.bam.bai | path/to/PD52103/PD52103n_lo0002.bam.bas | path/to/PD52103/PD52103n_lo0002.bam.met.gz |


##### Filter Indels

If Filter Indels is run (`params.run_filter_indel==true`), the samplesheet must include the additional columns `indel_vcf`, `indel_vcf_tbi`, `bas` and `met`. It should look something like this:

| sample_id       | match_normal_id | pdid    | bam_match                    | bai_match                        | indel_vcf                                                | indel_vcf_tbi                                                | bam                                 | bai                                     | bas                                     | met                                        |
| --------------- | --------------- | ------- | ---------------------------- | -------------------------------- | ------------------------------------------------------ | ---------------------------------------------------------- | ----------------------------------- | --------------------------------------- | --------------------------------------- | ------------------------------------------ |
| PD47151n_lo0002 | PD47151b        | PD47151 | path/to/PD47151/PD47151b.bam | path/to/PD47151/PD47151b.bam.bai | path/to/PD47151/PD47151n_lo0002.pindel.annot.vcf.gz | path/to/PD47151/PD47151n_lo0002.pindel.annot.vcf.gz.tbi | path/to/PD47151/PD47151n_lo0002.bam | path/to/PD47151/PD47151n_lo0002.bam.bai | path/to/PD47151/PD47151n_lo0002.bam.bas | path/to/PD47151/PD47151n_lo0002.bam.met.gz |
| PD47151n_lo0004 | PD47151b        | PD47151 | path/to/PD47151/PD47151b.bam | path/to/PD47151/PD47151b.bam.bai | path/to/PD47151/PD47151n_lo0004.pindel.annot.vcf.gz | path/to/PD47151/PD47151n_lo0004.pindel.annot.vcf.gz.tbi | path/to/PD47151/PD47151n_lo0004.bam | path/to/PD47151/PD47151n_lo0004.bam.bai | path/to/PD47151/PD47151n_lo0004.bam.bas | path/to/PD47151/PD47151n_lo0004.bam.met.gz |
| PD52103n_lo0002 | PD52103b        | PD52103 | path/to/PD52103/PD52103b.bam | path/to/PD52103/PD52103b.bam.bai | path/to/PD52103/PD52103n_lo0002.pindel.annot.vcf.gz | path/to/PD52103/PD52103n_lo0002.pindel.annot.vcf.gz.tbi | path/to/PD52103/PD52103n_lo0002.bam | path/to/PD52103/PD52103n_lo0002.bam.bai | path/to/PD52103/PD52103n_lo0002.bam.bas | path/to/PD52103/PD52103n_lo0002.bam.met.gz |



##### PHYLOGENETICS

If phylogenetics is run (`params.run_phylogenetics==true`), there are three scenarios:

1. If Filter SNV has been previously run, phylogenetics is run for SNVs, there is no additional column requirements than the requirements in Filter SNV
2. If Filter INDEL has been previously run, phylogenetics is run for INDELs.
- If Filter SNV has also been previously run, same as 1
- If Filter SNV has not been run, a column for tree topology is needed

| sample_id       | match_normal_id | pdid    | bam_match                    | bai_match                        | indel_vcf                                                | indel_vcf_tbi                                                | bam                                 | bai                                     | bas                                     | met                                        | topology |
| --------------- | --------------- | ------- | ---------------------------- | -------------------------------- | ------------------------------------------------------ | ---------------------------------------------------------- | ----------------------------------- | --------------------------------------- | --------------------------------------- | ------------------------------------------ | ------------------------ |
| PD47151n_lo0002 | PD47151b        | PD47151 | path/to/PD47151/PD47151b.bam | path/to/PD47151/PD47151b.bam.bai | path/to/PD47151/PD47151n_lo0002.pindel.annot.vcf.gz | path/to/PD47151/PD47151n_lo0002.pindel.annot.vcf.gz.tbi | path/to/PD47151/PD47151n_lo0002.bam | path/to/PD47151/PD47151n_lo0002.bam.bai | path/to/PD47151/PD47151n_lo0002.bam.bas | path/to/PD47151/PD47151n_lo0002.bam.met.gz | path/to/PD47151n.treefile |
| PD47151n_lo0004 | PD47151b        | PD47151 | path/to/PD47151/PD47151b.bam | path/to/PD47151/PD47151b.bam.bai | path/to/PD47151/PD47151n_lo0004.pindel.annot.vcf.gz | path/to/PD47151/PD47151n_lo0004.pindel.annot.vcf.gz.tbi | path/to/PD47151/PD47151n_lo0004.bam | path/to/PD47151/PD47151n_lo0004.bam.bai | path/to/PD47151/PD47151n_lo0004.bam.bas | path/to/PD47151/PD47151n_lo0004.bam.met.gz | path/to/PD47151n.treefile |
| PD52103n_lo0002 | PD52103b        | PD52103 | path/to/PD52103/PD52103b.bam | path/to/PD52103/PD52103b.bam.bai | path/to/PD52103/PD52103n_lo0002.pindel.annot.vcf.gz | path/to/PD52103/PD52103n_lo0002.pindel.annot.vcf.gz.tbi | path/to/PD52103/PD52103n_lo0002.bam | path/to/PD52103/PD52103n_lo0002.bam.bai | path/to/PD52103/PD52103n_lo0002.bam.bas | path/to/PD52103/PD52103n_lo0002.bam.met.gz | path/to/PD52103n.treefile |


3. If neither Filter SNV nor Filter Indel has been run,
- If phylogenetics need to be run for both SNVs, then use the topology obtained from SNV to assign INDEL mutations to (`params.snv_then_indel==true`), the required columns are `pdid`, `nr_path_snv`, `nv_path_snv`, `genotype_bin_path_snv`, `nr_path_indel`, `nv_path_indel` and `genotype_bin_path_indel`

| pdid | nr_path_snv | nv_path_snv | genotype_bin_path_snv | nr_path_indel | nv_path_indel | genotype_bin_path_indel |
| ---- | ------- | ------- | ------------ | ------- | ------- | ------------ |
| PD47151 | path/to/snv/PD47151/NR_somatic_noartefacts.txt | path/to/snv/PD47151/NV_somatic_noartefacts.txt | path/to/snv/PD47151/genotype_bin.txt | path/to/indel/PD47151/NR_somatic_noartefacts.txt | path/to/indel/PD47151/NV_somatic_noartefacts.txt | path/to/indel/PD47151/genotype_bin.txt |
| PD52103 | path/to/snv/PD52103/NR_somatic_noartefacts.txt | path/to/snv/PD52103/NV_somatic_noartefacts.txt | path/to/snv/PD52103/genotype_bin.txt | path/to/indel/PD52103/NR_somatic_noartefacts.txt | path/to/indel/PD52103/NV_somatic_noartefacts.txt | path/to/indel/PD52103/genotype_bin.txt |


- If a tree topology needs to be calculated before assigning mutations to the branches (`params.with_topology==false`, this is assumed for SNV), the required columns are `pdid`, `nr_path`, `nv_path` and `genotype_bin_path`. For example

| pdid | nr_path | nv_path | genotype_bin_path |
| ---- | ------- | ------- | ------------ |
| PD47151 | path/to/PD47151/NR_somatic_noartefacts.txt | path/to/PD47151/NV_somatic_noartefacts.txt | path/tp/PD47151/genotype_bin.txt |
| PD52103 | path/to/PD52103/NR_somatic_noartefacts.txt | path/to/PD52103/NV_somatic_noartefacts.txt | path/tp/PD52103/genotype_bin.txt |

- If a tree topology is specified to not be calculated (`params.with_topology==false`, this is assumed for Indels), the required columns are `pdid`, `nr_path`, `nv_path`, `genotype_bin` and `topology`. For example

| pdid | nr_path | nv_path | genotype_bin_path | topology |
| ---- | ------- | ------- | ------------ | -------- |
| PD47151 | path/to/PD47151/NR_somatic_noartefacts.txt | path/to/PD47151/NV_somatic_noartefacts.txt | path/tp/PD47151/genotype_bin.txt | path/to/PD47151.treefile |
| PD52103 | path/to/PD52103/NR_somatic_noartefacts.txt | path/to/PD52103/NV_somatic_noartefacts.txt | path/tp/PD52103/genotype_bin.txt | path/to/PD47151.treefile |



```bash
--input '[path to samplesheet file]'
```

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| REQUIRED COLUMNS FOR ALL SUBWORKFLOWS ---------------------------------------------------- |
| `pdid` |  Donor ID for your sample |                                                            |
| REQUIRED COLUMNS FOR CONPAIR, FILTER_SNV and/or FILTER_INDEL (`--run_conpair true` and/or `--run_filter_snv true` and/or `--run_filter_indel true`)---------------------------------------------------- |
| `sample_id`  | sample ID, must be unique |
| `match_normal_id` |  ID for your match normal sample |                                                            |
| `bam` | bam file for `sample_id`, must exist |                                                       |
| `bai` | tabix index file for `bam`, must exist |
| `bam_match` | bam file for `match_normal_id`, must exist |
| `bai_match` | tabix index file for `bam_match`, must exist |
| REQUIRED COLUMNS FOR FILTER_SNV and/or FILTER_INDEL PIPELINE (`--run_filter_snv true` and/or `--run_filter_indel true`) ---------------------------------------------------- |
| `bas` | bam status file for `bam`, must exist |
| `met` | met (samtools markedup) file for `bam`, must exist |
| REQUIRED COLUMNS FOR FILTER_SNV (`--run_filter_snv true`) ---------------------------------------------------- |
| `snv_vcf` | VCF file for the SNVs of `sample_id`, must exist |
| `snv_vcf_tbi` | tabix index file for VCF file for `snv_vcf_tbi`, must exist |
| REQUIRED COLUMNS FOR FILTER_INDEL (`--run_filter_indel true`) ---------------------------------------------------- |
| `indel_vcf` | VCF file for the indels of `sample_id`, must exist |
| `indel_vcf_tbi` | tabix index file for VCF file for `indel_vcf_tbi`, must exist |
| REQUIRED COLUMNS FOR PHYLOGENETICS FOR BOTH SNVs AND INDELs WITHOUT FILTERING SNVs OR INDELs (`--run_phylogenetics true --run_filter_snv false --run_filter_indel false --snv_then_indel true`) ---------------------------------------------------- |
| `nr_path_snv` | NR file (depth at the variant locus) for SNVs for `pdid`, must exist |
| `nv_path_snv` | NV file (reads supporting the variant) for SNVs for `pdid`, must exist |
| `genotype_bin_path_snv` | binary genotype file for SNVs for `pdid`, must exist |
| `nr_path_indel` | NR file (depth at the variant locus) for INDELs for `pdid`, must exist |
| `nv_path_indel` | NV file (reads supporting the variant) for INDELs for `pdid`, must exist |
| `genotype_bin_path_indel` | binary genotype file for INDELs for `pdid`, must exist |
| REQUIRED COLUMNS FOR PHYLOGENETICS or PHYLOGENETICS-GIVEN-TREE-TOPOLOGY WITHOUT RUNNNING FILTERING SNVs OR INDELs (`--run_phylogenetics true --run_filter_snv false --run_filter_indel false --snv_then_indel [false/null]`) ---------------------------------------------------- |
| `nr_path` | NR file (depth at the variant locus) for `pdid`, must exist |
| `nv_path` | NV file (reads supporting the variant) for `pdid`, must exist |
| `genotype_bin_path` | binary genotype file for `pdid`, must exist |
| REQUIRED COLUMNS FOR PHYLOGENETICS-GIVEN-TREE-TOPOLOGY WITHOUT RUNNNING FILTERING SNVs (`--run_phylogenetics true --run_filter_snv false --run_filter_indel true`, or `--run_phylogenetics true --run_filter_snv false --run_filter_indel false --with_topology true`) ---------------------------------------------------- |
| `topology` | tree topology file for `pdid`, must exist |

## Other inputs

| Input Parameters    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| REQUIRED PARAMS (must not be null, ie if no default, need to be defined by users in `nextflow run`) |
| `outdir`  | type: directory-path, The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure |
| `with_match_normal` |  type: boolean, default: true, whether the samples are provided with real (non-synthetic) match normal  |                                                            |
| `run_conpair` |  type: boolean, default: true, whether to run conpair to check for contamination  |                                                            |
| `run_filter_snv` |  type: boolean, default: true, whether to run filtering SNVs  |                                                            |
| `run_filter_indel` |  type: boolean, default: true, whether to run filtering INDELs  |                                                            |
| `run_phylogenetics` |  type: boolean, default: true, whether to run the phylogenetics pipeline (ie build phylogenetic tree and/or assign mutations to the branches/nodes)  |                                                            |
| REQUIRED PARAMS IF RUNNING CONPAIR |                                                          |
| `concordance_threshold` |  type: number, default: 90, threshold above which 2 samples are considered from the same donor  |                                                            |
| `contamination_threshold_samples` |  type: number, default: 0.3, contamination threshold above which a sample of interest is considered contaminated   |                                                            |
| `contamination_threshold_match` |  type: number, default: 5, contamination threshold above which a match normal sample is considered contaminated (this won't be used if `with_match_normal == false`)  |                                                            |
| REQUIRED PARAMS IF RUNNING FILTER_SNV |                                                          |
| `snv_vcfilter_config` |  type: file-path, default: "$projectDir/data/default/snv_default.filter", Path to vcfilter config file, default to data/default/snv_default.filter  |                                                            |
| `snv_rho_threshold` |  type: number, default: 0.1, rho threshold in the beta-binomial test below which variant is filtered out  |                                                            |
| `sigprofiler_genome` |  type: string, genome to input to sigprofiler matrix generator   |                                                            |
| REQUIRED PARAMS IF RUNNING FILTER_INDEL |                                                          |
| `indel_vcfilter_config` |  type: file-path, default: "$projectDir/data/default/indel_default.filter", Path to vcfilter config file, default to data/default/indel_default.filter  |                                                            |
| `indel_rho_threshold` |  type: number, default: 0.2, rho threshold in the beta-binomial test below which variant is filtered out  |                                                            |
| `sigprofiler_genome` |  type: string, genome to input to sigprofiler matrix generator   |                                                            |
| REQUIRED PARAMS IF RUNNING PHYLOGENETICS WITHOUT RUNNING FILTER_SNV OR FILTER_INDEL |                                                          |
| `snv_then_indel` |  type: boolean, whether running phylogenetics for SNVs (phylogenetics) followed by phylogenetics for INDELS (phylogenetics_provided_topology)   |                                                            |
| `provided_topology` |  type: boolean, when only one phylogenetics pipeline is run, whether it is provided a topology   |                                                            |
| `phylogenetics_outdir_basename` |  type: string, default: phylogenetics_snp_out if `with_topology == false` and phylogenetics_indel_out if `with_topology == true`, directory name to publish results when only one phylogenetics pipeline is run   |                                                            |
| OPTIONAL PARAMS |
| `image_cachedir`  | type: directory-path, /path/to/lcmb/image if a container runtime profile or is enabled or conf/container.config is loaded, eg with --profile singularity |
| `use_custom_genome` |  type: boolean, default: false, whether to use custom genome as specified in conf/custom_ref_genomes.config  |
| REFERENCE GENOME REQUIRED PARAMS |                                                          |
| `fasta` | type: file-path, Path to FASTA genome file |
| `fai` | type: file-path, Path to FASTA genome index (fa.fai) file. |
| `fasta_dict` | type: file-path, Path to FASTA genome dictionary (fa.dict) file, used in conpair pileup step |
| `marker_bed` | type: file-path, Path to marker bed file, used in the conpair step for contamination checking. |
| `marker_txt` | type: file-path, Path to marker txt file, used in the conpair step for contamination checking. |
| `hidepth` | type: file-path, Path to high depth region, used in cgpVaf under FILTER_SNV and FILTER_INDEL subworkflow |
| `hidepth_tbi` | type: file-path, Path to high depth region index (bed.tbi), used in cgpVaf under FILTER_SNV and FILTER_INDEL subworkflow |
| REFERENCE GENOME OPTIONAL PARAMS |                                                          |
| `genome` | type: string, If using a reference genome configured in the pipeline using `--use_custom_genome true` (ie `conf/custom_ref_genomes.config` is loaded) or iGenomes, use this parameter to give the ID/label for the reference. This is then used to build the full paths for all required reference genome files e.g. `--genome GRCh38`.  |
| `igenomes_ignore` | type: boolean, default: true, Do not load `igenomes.config` when running the pipeline. You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config` |
| `igenomes_base` | type: directory-path, default: "s3://ngi-igenomes/igenomes/", The base path to the igenomes reference files, ignored when `--igenomes-ignore` is enabled |
| `custom_genome_base` | type: directory-path, default: null, The base path to the custom genome reference files (REQUIRED IF `--use_custom_genome true` AND if used to build the full path in the reference genome profile in [conf/custom_ref_genomes.config](conf/custom_ref_genomes.config)) |

You can then run the pipeline like this

```bash

nextflow run /path/to/lcmb/main.nf \
  -profile <docker/singularity/.../institute> \
   --input /path/to/samplesheet.tsv \
   --with_match_normal true \
   --run_conpair true \
   --run_filter_snv true \
   --run_filter_indel true \
   --run_phylogenetics true \
   --hairpin_genome hg38 \
   --sigprofiler_genome GRCh38 \
   --use_custom_genome true \
   --genome genome_label_in_custom_genome_config \ #eg your_genome_label
   --outdir /path/to/outdir
```

## Running the pipeline

You can specify a `-profile`  option as follows

```bash

nextflow run /path/to/lcmb/main.nf \
   -profile <docker/singularity/.../institute> \
   --input /path/to/samplesheet.tsv \
   ...
   --outdir /path/to/outdir \
   -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run /path/to/lcmb/main.nf -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.tsv'
outdir: './results/'
genome: 'GRCh38'
<...>
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
cd /path/to/lcmb/dir
git pull
```

<!-- ### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/lcmb releases page](https://github.com/nf-core/lcmb/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter. -->

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

