# MMOSurv

MMOSurv: Meta-learning for few-shot survival analysis with multi-omics data

This is an implementation of MMOSurv in Python 3.9.5 under Linux with CPU Intel(R) Core(TM) i9-12900K and GPU NVIDIA GeForce RTX 3090. It follows a modern deep learning design and is implemented by PyTorch platform.

# Environment Installation

$ conda create -n survival-analysis

$ source activate survival-analysis

$ pip 


# Data and Preprocessing

The TCGA dataset are publicly available and can be obtained from: https://www.cancer.gov/ccg/research /genome-sequencing/tcga. We put raw data (gene/microRNA expression data and clinical data) of some cancer patients in the Data_preprocessing/Origin_Data directory .

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

git clone https://github.com/LiminLi-xjtu/MMOSurv.git

cd MMOSurv

python mvml_survival_analysis.py
