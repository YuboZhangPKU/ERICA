# ScriptsForModelEvaluation
Below we describe each file in ScriptsForModelEvaluation.zip

## 0DataSimulation.sh
The commands are used for data simulation, data labeling, and format conversion. The scripts used in the workflow are listed below: 

### FourTaxonSimulation.py, AsymmetricFiveTaxonSimulation.py, and SymmetricFiveTaxonSimulation.py
The Python scripts are used for simulating multiple sequence alignment (MSA) datasets of four-taxon, asymmetric five-taxon, and symmetric five-taxon, respectively. A set of demographic scenarios of various phylogenies, divergence times, and gene flow events are included. Optional parameters include mutation rates, recombination rates, the number of haplotypes, and the length of sequences.

The script requires external software `ms` and [`Seq-Gen`](http://tree.bio.ed.ac.uk/software/seqgen/).

### FourTaxonSimulationf.py
The Python script is used for simulating datasets under scenarios of an instant admixture event, with different intensities, times, and directions of gene flow.

The script requires external software `ms` and [`Seq-Gen`](http://tree.bio.ed.ac.uk/software/seqgen/).

### FourTaxonSimulationAI.py, and SymmetricFiveTaxonSimulationAI.py
The Python scripts are used for simulating datasets of four-taxon, and symmetric five-taxon, respectively. Selection is introduced in the models.

The script requires external software [`msms`](https://www.mabs.at/publications/software-msms/) and [`Seq-Gen`](http://tree.bio.ed.ac.uk/software/seqgen/).

### Genealogy2Weights.py
The Python script is used for calculating the proportion of each topology using the topology weighting method.

The script requires external software [`Twisst`](https://github.com/simonhmartin/twisst/tree/master).

### ChangeSequenceOrder.py
The Python script is used for reordering the MSAs.

### Phylip2Geno.py
The Python script is used for converting the MSAs into geno format for subsequent analysis.

### ErrorGeneration.py
The Python script randomly introduces sequencing errors and missing data into the MSAs.


## 1ModelTrainingEvaluation.sh
The shell script contains commands used for ERICA model training, evaluation, and data analysis using other software, including window-tree-based topology weighting, *D* and *fd* statistics of ABBA-BABA test, DFOIL, IBDmix, SPrime, and genomatnn. 

The script requires software [`Twisst`](https://github.com/simonhmartin/twisst/tree/master), [`genomics_general`](https://github.com/simonhmartin/genomics_general), [`DFOIL`](https://github.com/jbpease/dfoil), [`vcftools`](https://vcftools.github.io/), [`IBDmix`](https://github.com/PrincetonUniversity/IBDmix), [`SPrime`](https://github.com/browning-lab/sprime), and [`genomatnn`](https://github.com/grahamgower/genomatnn).

The scripts used in the workflow are listed below:
### windows_MBB.py, and moving_block_bootstrap.R
The scripts are used to filter windows with *D* and *fd* values significantly different from zero. The standard errors are calculated using a moving block bootstrap approach, and the p values are calculated using two-tailed z-tests.

### sliding_window_dfoil.py
The Python script is used for running DFOIL analysis with sliding windows.

### geno2vcf.py
The Python script is used for converting the data from geno format to vcf format for subsequent analysis.

### RandomSample.py
The Python script is used for randomly sampling one individual from the donor population to satisfy the assumptions of the method IBDmix。

### AdaptiveTest.map
The genetic map is used for SPrime analysis. The map distance is calculated based on the recombination rate used in data simulation.  


## 2RealGenomicDataProcessing.sh
The shell script contains commands used for processing real genomic data of genus *Heliconius* and genus *Oryza*. Specifically, the commands include genotype calling using genome re-sequencing or pan-genome data, demographic modeling using G-PhoCS, phylogenetic analysis, and gene expression analysis.
The script requires software [`Trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic), [`Bowtie 2`](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [`picard-tools 1.96`](https://broadinstitute.github.io/picard/), [`GATK 3.7`](https://gatk.broadinstitute.org/hc/en-us), [`G-PhoCS`](http://compgen.cshl.edu/GPhoCS/), [`vcf2phylip`](https://github.com/edgardomortiz/vcf2phylip), [`RAxML`](https://cme.h-its.org/exelixis/web/software/raxml/index.html), [`tophat2`](https://github.com/DaehwanKimLab/tophat2), and [`htseq-count`](https://htseq.readthedocs.io/en/master/).

### run_gphocs.py, rice_G-PhoCS_dataset_7_indv_Omeri.txt, and gphocs_rice_7_pop_nomig_h.ctl
The data file `rice_G-PhoCS_dataset_7_indv_Omeri.txt` is the sequence of neutral and independent loci used for demographic modeling. The control file `gphocs_rice_7_pop_nomig_h.ctl` is the G-PhoCS model without any migration band. The Python script is used for running G-PhoCS analysis to test all possible migration bands between current species independently.

### gphocs_rice_7_pop_mig_3.ctl
The control file `gphocs_rice_7_pop_mig_3.ctl` is the full model with all significant migration bands.
