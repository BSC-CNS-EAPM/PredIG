### Datasets
Antigen types, experimental Sources & Immunogenicity Label Harmonization
To obtain a sizeable dataset for machine-learning modelling, we retrieved 15K pairs of epitope and HLA-I allele (pHLAs) with validated immunogenicity, from public databases47,49,63 and from recent literature4,36,48,64 As for disease of origin, our dataset contains epitopes derived from pathogens and from cancer including neoantigens, tumor-associated antigens and non-canonical cancer antigens. Experimental validations comprise T-cell reactivity assays such as IFN-y Elispot, 41BB staining or TNFa-release and T-cell binding assays such as MHC-tetramers and multimers. The experimental results are considered ground truth and used as immunogenicity label. Pairs of epitope and HLA allele (pHLA) validated as as "positive", "positive-low", "positive-intermediate" and "positive-high" are harmonized as "positive" for binarization of the label. pHLAs labeled as "negative" maintain the same label. Among pHLA cases with more than a single validation experiment, those with at least one responsive subject were considered positive and the remaining negative instances were discarded to avoid discording/discrepant immunogenicity labels for the same pHLA. This criterion is based on the fact that each T-cell assay replicate/experiment uses a unique T-cell sample with polyclonal TCRs and diverse T-cell states. Thus, one finding of immunogenicity validates the ground truth of a pHLA instance whereas a negative finding does not invalidate it because it can imply the absence of the responsive TCR or the presence of non-effector T-cell states in the sample.

###Data Curation
All datasets were curated in a series of general quality filters and a subset of dataset-specific filters. The general filters include HLA nomenclature standardization to discard alleles/cases with missing or insufficient information (explained below); epitope length between 8 and 15 amino acids and at least one patient tested for all cases and one patient responding for positive cases. Duplicates between datasets were removed based on pHLA identity (ie. concatenating epitope sequence and standardized HLA allele) maintaining the positive instances in cases of discording annotation.
We adapted HLA allele ontologies66 to focus on experiments performed on mono-allelic cell lines or on epitopes predicted as binders for a single HLA-I allele of the antigen presenting cells used in the assay. IEDB49 data was curated as follows: host organism was restricted to “Homo sapiens” to only include epitopes tested against human T-cells; MHCType was restricted to "MHC-I"; the dataset was split by Epitope.Relationship into "neo-epitope" for cancer datasets and "pathogen" for pathogenic datasets and MHC allele restriction was refined using HLA Allele Ontologies66; Allele.Evidence.Code was used to discard "not determined" instances. Data from LANL84 and HCV databases85 was retrieved using Repitope62). and curated as in IEDB. TANTIGEN v1.063 and v2.047 data only contains positive cases, the number of subjects tested and the experimental source are annotated in the original publication. In PRIME v1.036, we discarded data coming from the studies "Calis", "Dengue" and "Random" for HLA resolution limitations; the number of subjects tested, responded and experimental source depend on the specific study36. The non-canonical tumoral antigen dataset was retrieved from Gros et al64 and curated by the "nonC-TL" antigen category to include alternative reading frames (OffFrame), intronic, intergenic, non-coding regions (non-CDS) and 5' and 3' (UTR5, UTR3). This dataset was validated by IFN- y elispots and 41BB stainings. SARS-CoV-2 data retrieved from Schumacher et al.4 was curated using general filters exclusively and validated by pHLA multimer assays. TESLA data for cancer neoantigens was curated using general filters and validated using pHLA multimers48. Tran et al data including cancer neoantigens from gastrointestinal tumors was curated using general filters and validated by IFN- y elispot and 41BB staining65. 

###HLA nomenclature standardization
We only retained those pHLA points with 4-digits of HLA resolution, necessary to perform reliable binding predictions (Fig 2B) and their nomenaclature was standardized to the HLA allele nomenclature established by the WHO Nomenclature Committee for Factors of the HLA System86, which is as follows: HLA-(Gene)*(Allele group) : (specific HLA protein or allotype). For instance, HLA-A*02:01. Data points with insufficient allelic annotation were discarded (ie. HLA-A2). Data points with deeper resolution were restricted to non-synonymous modifications specified by 4-digits (ie. HLA-A*01:01:03 > HLA-A*01:01). All annotations that did not follow the standardized symbols were harmonized into 4-digits (ie. HLA-A0101 > HLA-A*01:01).

### Data Splits
The training set includes data curated from the IEDB49, LANL84, HCV85,87, TANTIGEN v1.063 and PRIME v1.036. This set contains epitopes of mixed origin including cancer associated antigens, neoantigens and pathogens. The independent sets (held-outs) include data from recent publications4,64, TESLA consortium48 and TANTIGEN v2.047. These include a cancer neoantigen dataset (Independent1), a non-canonical cancer antigen dataset (Independent 2) and a pathogen dataset containing SARS-CoV-2 T-cell epitopes (Independent 3). To avoid data leakage, none of the independent sets contains any pHLA found in the training set.