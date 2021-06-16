**A reproducible pipeline for polishing DNA reads, with a particular focus on sequence data for historical specimens**.

[TOC]

## Introduction

**nf-polish** is a bioinformatics pipeline for NGS sequencing data and is particularly suited to process reads from historical (museum) specimens.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple computing infrastructures in a very portable manner. It currently creates a custom Conda environment from scratch and therefore does not require any pre-compiled software packages (besides Nextflow itself). Future development may include containerized versions as well to further enhance reproducibility. If possible, e.g. when running on a HPC cluster, the pipeline will process reads in parallel and all batch submission jobs are handled internally through the Nextflow workflow. The pipeline converts raw sequencing data (FASTQ inputs) into fully polished reads (FASTQ output), and records changes in read number and length along the way.

## Quick Start

1. Install [`nextflow`](https://www.nextflow.io/) (version >= 19.04) 
2. Install [`Conda`](https://conda.io/miniconda.html) (version >= 4.10) 
3. Download the pipeline, create a nextflow config profile that matches your cluster set-up ( [`profile`]( https://www.nextflow.io/docs/latest/config.html#config-profiles) ) and start running your own analysis! If you want to use the existing `rackham.config` remember to specify your SNIC project ID (format: `snic20XX-XX-XXX`) as well as the path to `nf-polish/environment.yml`

    ```bash
    nextflow run nf-polish/main.nf -profile mfn --reads /path/to/'*_R{1,2}_L001.fastq.gz'
    ```
4. Once your run has completed successfully, clean up the intermediate files.

    ```bash
    nextflow clean -f -k
    ```

**N.B.** Please make sure that:
* The scripts in `./bin/xxx.py` can be executed by any user. You can check the user permissions via:

    ```bash
    ls -l ./bin/
    ```
And if need be changed via:

    ```bash
    cd ./bin/
    chmod 755 *
    ```

* There is **sufficient storage capacity** on the location where you will execute the pipeline from (and where the nextflow work directory will be stored). On the MfN cluster, running on `/home/` will easily lead to problems and all pipelines will need to be executed from a `/data4/` project folder.

* Please use the `nextflow.config` file or the `--adapters ` option to specify the location of a fasta file with sequencing adapters that may have been used. E.g. the adapter file as listed in the `Trimmomatic` adapter folder.

## Pipeline Summary

### Default Steps

By default the pipeline currently performs the following:

* Sequencing quality control (`FastQC`)
* Removal of PCR duplicates (`HTStream/hts_SuperDeduper`) [optional]
* Adapter Trimming (`Trimmomatic`)
* Read merging (`PEAR`)
* Quality Trimming (`Trimmomatic`)
* Removal of low-complexity reads (`./bin/remove_low_complex.py`)
* Calculating processing statistics (`seqkit`)
* Gather and visualise processing statistics (`./bin/parse_visualize_stats.py`)

### Input
The pipeline takes in raw, demultiplexed, sequence reads which should all be stored in a single folder. It can currently only polish paired-end reads and reads do not need to adhere to a specific naming format. Reads that belong to the same paired-end read set should have identical names, except for using `_R1_` and `_R2_` (or `_R1.`, `.R1_`, etc.) to identify the forward and reverse read of the same library.

Run specific options can  be specified in the `./nextflow.config` script. All other values for programs are set at default. Boolean values can be specified as `true` or `false` and path names need to be absolute. If preferred, the same options listed in the config file can also be directly modified by using a parameter flag when initiating the nextflow run. Example given:

    ```bash
    nextflow run nf-polish/main.nf -profile mfn --reads /path/to/'*_R{1,2}_L001.fastq.gz' --outdir /path/to/results --skip_dedup false
    ```

### Output
The pipeline currently stores the FASTQ read files during each step (which may be changed in future versions) and the processed output folders are sequentially numbered. Within the folder for each processing step, reads are stored by library ID (e.g. the prefix of the original read file) and a short summary overview of read statistics is included. Following read polishing, processing statistics (read number and length) are gathered across all processing steps, and visualised in a summary overview.
