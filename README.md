# MMOSurv

MMOSurv: Meta-learning for few-shot survival analysis with multi-omics data

MMOSurv enables to learn an effective multi-omics survival prediction model from a very few training samples of a specific cancer type, with the meta-knowledge across tasks from relevant cancer types. This is an implementation of MMOSurv in Python 3.9.5 under Linux with CPU Intel(R) Core(TM) i9-12900K and GPU NVIDIA GeForce RTX 3090. It follows a modern deep learning design and is implemented by PyTorch platform.

# Environment Installation

$ conda env create -f environment.yml

$ source activate survival_analysis

# Data Download

The TCGA dataset are publicly available and can be obtained from: https://www.cancer.gov/ccg/research/genome-sequencing/tcga.  We downloaded genomic data and survival data using the R package TCGAbiolinks.

Enter the R environment

$ setwd("Data_preprocessing/R")

$ source("miRNA_download.R")

$ source("tcga_gene_expression_download.R")

$ source("tcga_clinical_download.R")

# Data Preprocessing

#### cd MMOSurv/Data_preprocessing

The python program (extract_RNASeq_expression.py) summarized the gene expression of patients with the same type of cancer into a csv file. 

#### python extract_RNASeq_expression.py

#### python extract_miRNA_expression.py

The python program (RNASeq_variable_filter.py) imputed missing data, filtered out noise-sensitive features (those whose values ​​remained almost constant across samples), and performed a logarithmic transformation.

#### python RNASeq_variable_filter.py

#### python microRNA_variable_filter.py

The python program (clinical_preprocess.py)  aggregated the survival data of cancer patients (such as censoring indicator δ and observed time O).

#### python clinical_preprocess.py

The python program (matched_patients.py) only kept the patients with matched gene expression, microRNA expression data and survival data.

#### python matched_patients.py

# Run the main routine

MMOSurv first learns a suitable initialization of parameters for the multi-omics survival model from multi-omics data of relevant cancers, and then adapts the parameters quickly and efficiently for the target cancer task with a few training samples.

#### git clone https://github.com/LiminLi-xjtu/MMOSurv.git

#### cd MMOSurv

#### python mvml_survival_analysis.py
