# PredIG: an interpretable predictor of T-cell epitope immunogenicity

Roc Farriol-Duran<sup>1,2</sup>*, Miguel Vázquez<sup>1</sup>, Eduard Porta-Pardo<sup>1,2</sup>, Víctor Guallar<sup>1,3</sup>ª
1. Barcelona Supercomputing Center (BSC), 08034 Barcelona, Spain.
2. Josep Carreras Leukaemia Research Institute (IJC), Badalona, Spain
3. Institució Catalana de Recerca i Estudis Avançats (ICREA), Barcelona, Spain

\* First author <br>
ª Corresponding author <br>
For scientific or usage enquires refer to: roc.farriol@bsc.es and/or victor.guallar@bsc.es

### Abstract
Cytotoxic T cells are key effectors in the immune response against pathogens and cancer. Hence, their activation, driven by the recognition of immunogenic epitopes, consitutes a coveted goal for immunotherapies. However, the epitope landscape, both in cancer and infection, is too large to test due to the immense number of candidates versus the high cost and low throughput of experimental techniques. Enabling larger throughtputs, immunoinformatic models prioritize the candidates with greater potential but their success rate has remained incremental and their explainability limited. Here we present PredIG, a predictor of T-cell epitope immunogenicity that integrates antigenic and physicochemical properties of 17448 pHLA-I using XGBoost, a decision-tree-based algorithm that boosts explainability. PredIG outperforms state-of-the-art methods in two pathogen and non-canonical cancer antigen held-out sets. In cancer neoantigens, PredIG increases the success rate of binding affinity predictions and identifies alternative immunogenic epitopes. Our XAI scheme pinpoints the importance of antigenic and physicochemical epitope properties and their differences in each antigen type. Overall, PredIG can increase the immunogenicity success rates in vaccine design for cancer and infection and displays an unprecedented interpretability to build community trust. Plus, its containerized environments and a user-friendly webserver grant PredIG's accessibility at https://horus.bsc.es/predig

### Graphical Abstract
![Alt text](images/predig_graph_abstract.jpg)

### Usage Scheme
##### PredIG usage modes in a user-friendly webserver implementation and in containerized environments for high-throughput reproducibility in HPC environments.

<p align="center">
    <img src="images/fig6_usage_scheme.jpg" alt="Usage Scheme" width="60%">
</p>

**Exploration Modes (Inputs)** <br>
A) **"CSV-Uniprot” mode**: input a .CSV file with pairs of peptide and HLA-I allele and the Uniprot ID of the corresponding parental protein. <br>
B) **"CSV-Recombinant" mode**: input a .CSV file with pairs of peptide and HLA-I allele and the amino acid sequence of the protein of origin. This mode is designed to support (recombinant) proteins without Uniprot ID but can also work with any protein sequence. <br>
C) **"FASTA" mode**: input a FASTA file with the target protein sequence and a .CSV file with a list of HLA-I alleles of interest ("HLA_allele" column). By default, PredIG will generate all possible epitopes of 8 to 14 AA of length and will calculate against the input HLA-I alleles. <br>

**PredIG Model Selection: antigen-type specific** <br>
D) The user can choose between three PredIG predictive models: PredIG-NeoA optimized for cancer neoantigens, PredIG-Non-Can for non-canonical cancer antigens and PredIG-Path for pathogen antigens. E) PredIG's output is a CSV with one pHLA-I per row containing PredIG score and all the predictors in PredIG feature space. <br>

Find PredIG's models in the data folder of this repo.  <br>

### Datasets
Find PredIG's train, test and held-out datasets in the data folder of this repo. <br>
These include: <br>
- predig_train_modf.csv > PredIG train set <br>
- predig_test_modf.csv > PredIG test set for model validation. <br>
- predig_i1_modf.csv > PredIG held-out for cancer neoantigen generalization assessment. i1 refers to Independent 1. <br>
- predig_i2_modf.csv > PredIG held-out for non-canonical cancer antigen generalization assessment. i2 refers to Independent 2. <br>
- predig_i3_modf.csv > PredIG held-out for pathogen generalization assessment. Contains epitopes from SARS-CoV-2. i3 refers to Independent 3. <br>



### Tutorial


