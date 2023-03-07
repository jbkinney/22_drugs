# Alternative splicing changes induced by risdiplam and branaplam in HeLa cells

This section of the repository contains the [nextflow](https://www.nextflow.io/) pipeline used to analyze the alternative splicing changes induced by risdiplam and branaplam across the whole transcriptome using RNA-seq data to study their 5'splice site sequence dependent effect. 

It employs a general pipeline to count reads mapping to alternative splicing events based on [rMATS](https://rnaseq-mats.sourceforge.net) and performs Quality Control on the sequenced samples using [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [rseqc](https://pythonhosted.org/RSeQC/). Next, we perform more adhoc quantitative analysis to infer exon and splice site sequence level drug effects using [AS_quant](https://bitbucket.org/cmartiga/as_quant/src). 

### Dependencies

These are the required tools to re-run our analysis:

- nextflow 21.10.0
- STAR 2.8.8a
- rMATS 4.1
- rseqc 4.0.0
- fastqc 0.11.9
- gffread 0.12.1
- multiqc 1.0.dev0
- [AS_quant](https://bitbucket.org/cmartiga/as_quant/src) with commit id 5bdc0284713bbe4a4c9341f4dc2de6fcf347e5f8


Create a new python environment and install the package manually for [AS_quant](https://bitbucket.org/cmartiga/as_quant/src)

```bash
conda create -n as_quant python=3.7
git clone git@bitbucket.org:cmartiga/as_quant.git
cd as_quant
git checkout 5bdc0284713bbe4a4c9341f4dc2de6fcf347e5f8 # if necessary
python setup.py install
```

Create other environments for the other tools, either manually or using the provided file. Note that it may take some time to solve all the dependencies

```bash
conda create -f environment.yml
```

### Download raw data

Download the data prior to running the analysis

```bash
cd fastq
bash download.sh
```

### Run the pipeline

The pipeline downloads genome references directly and performs all analysis automatically. However, the paths to the environments with the different tools will need to be updated beforehand. After that, it can be run simply as follows:

```bash
conda activate nextflow
nextflow pipeline.nf -c nextflow.config
```

The pipeline can be run locally as well as in a computer cluster thanks to the workflow manager easily, but we strongly recommend using a computer cluster due to high computational requirements. In our system, the whole pipeline took 1-2 days depending on resource availability. We used the configuration to run on the cluster specified by the `nextflow.config` file on a SGE scheduler. You can check [nextflow](https://www.nextflow.io/) documentation for more details such as configuration in other job schedulers.


Producing as main output files

- exons.psi.csv: file containing the estimated PSI and their posterior credible intervals in each of the experimental conditions for each exon with sufficient read coverage
- allelic_manifolds.csv: file containin the estimated effect of each drug for every 5'splice site sequence that was observed at least 10 times in the data

