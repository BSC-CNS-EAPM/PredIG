# PredIG: an interpretable predictor of T-cell epitope immunogenicity

Roc Farriol-Duran1,2*, Miguel Vázquez1, Eduard Porta-Pardo1,2, Víctor Guallar1,3ª
1. Barcelona Supercomputing Center (BSC), 08034 Barcelona, Spain.
2. Josep Carreras Leukaemia Research Institute (IJC), Badalona, Spain
3. Institució Catalana de Recerca i Estudis Avançats (ICREA), Barcelona, Spain

\* First author <br>
ª Corresponding author <br>
For scientific or usage enquires refer to: roc.farriol@bsc.es and/or victor.guallar@bsc.es

## Abstract
Cytotoxic T cells are key effectors in the immune response against pathogens and cancer. Hence, their activation, driven by the recognition of immunogenic epitopes, consitutes a coveted goal for immunotherapies. However, the epitope landscape, both in cancer and infection, is too large to test due to the immense number of candidates versus the high cost and low throughput of experimental techniques. Enabling larger throughtputs, immunoinformatic models prioritize the candidates with greater potential but their success rate has remained incremental and their explainability limited. Here we present PredIG, a predictor of T-cell epitope immunogenicity that integrates antigenic and physicochemical properties of 17448 pHLA-I using XGBoost, a decision-tree-based algorithm that boosts explainability. PredIG outperforms state-of-the-art methods in two pathogen and non-canonical cancer antigen held-out sets. In cancer neoantigens, PredIG increases the success rate of binding affinity predictions and identifies alternative immunogenic epitopes. Our XAI scheme pinpoints the importance of antigenic and physicochemical epitope properties and their differences in each antigen type. Overall, PredIG can increase the immunogenicity success rates in vaccine design for cancer and infection and displays an unprecedented interpretability to build community trust. Plus, its containerized environments and a user-friendly webserver grant PredIG's accessibility at https://horus.bsc.es/

## Graphical Abstract
Embed image


## Tutorial


