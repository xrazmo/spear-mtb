## SPEAR-MTB: en*S*emble *P*r*E*diction of *A*ntibiotic *R*esistance in *M*ycobacterium *T*u*B*erculosis

This repository has employed the following pipelines in order to arrive at a consensus on drug-resistant _Mycobacterium tuberculosis_ predictions:

- **[CRyPTIC](https://github.com/iqbal-lab-org)**: The pipeline involves decontaminating reads, performing quality controls and variant calling using various tools and the the H37Rv reference genome. We have incorporated the pipeline into the nextflow DSL2 and employed precompiled reference genome databases (e.g., H37Rv, nontuberculous Mycobacterium, human, etc.) to facilitate analyzing input genomes. Moreover, the [gnomonicus](https://github.com/oxfordmmm/gnomonicus) repository is used to predict durg reistance based on the detected mutations and catalogues available. The spear-MTB utilizes _WHO_ and _CRyPTIC_ catalogues, which are presented in [GARC1](https://fowlerlab.org/2018/11/25/goarc-a-general-ontology-for-antimicrobial-resistance-catalogues/) format and maintained [here](https://github.com/oxfordmmm/tuberculosis_amr_catalogues).

- **[TB Profiler](https://github.com/jodyphelan/TBProfiler)**:
  The pipeline includes aligning reads to the H37Rv reference genome, utilizing a pairwise aligner and then identifying variations with the use of bcftools. The catalogue of mutations is in [hgvs nomenclature](http://varnomen.hgvs.org/bg-material/simple/) and is maintained [here](https://github.com/jodyphelan/tbdb).

  > &#128712;**Info:**
  > A catalogue containing suggested mutaions by public health sweden (Fohm) will be added soon.

> &#x26A0; **Warning:**
> The spear-MTB pipeline is only tested using Illumina paired-end reads and on Unix-based system having Singularity.

## **Installation**

Spear-MTB needs you to install the following programs:

- Nextflow
- Singularity
- conda

After the successful installations, proceed with running the following commands:

```
git clone https://github.com/xrazmo/spear-mtb.git
cd ./spear-mtb | chmod +x ./setup.sh | ./setup.sh
```

These commands perform the following steps:

- Downloading the repository from Github.
- Creating a conda environment named _spear-mtb_.
- Downloading and extracting the precompiled reference database for the CRyPTIC pipeline.
- Downlaoding the latest version of catalogues.

## **Usage**

### - Nextflow config file:

Setting up the config file...

```

```

### - Running the pipeline:

```
conda activate spear-mtb
nextflow run ../spear-mtb/main.nf -profile [slurm or interactive] -w [working-directory] --input_dir [Reads-directory] --out_dir [output-directory]
```

### - Merging outputs:

```
nextflow run ../spear-mtb/main.nf -profile [slurm or interactive] -w [working-directory] --out_dir [output-directory] -entry generate_report
```


## **Report**

HTML template...
