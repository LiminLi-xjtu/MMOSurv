# MMOSurv

MMOSurv: Meta-learning for few-shot survival analysis with multi-omics data

This is an implementation of MMOSurv in Python 3.9.5 under Linux with CPU Intel(R) Core(TM) i9-12900K and GPU NVIDIA GeForce RTX 3090. It follows a modern deep learning design and is implemented by PyTorch platform.

# Installation

$ pip install pandas=1.3.5, numpy=1.20.2, scikit-learn=1.0.2, pytorch=1.9.0, lifelines==0.27.0


# Data and Preprocessing

The TCGA dataset are publicly available and can be obtained from: https://www.cancer.gov/ccg/research /genome-sequencing/tcga. We put raw data (gene/microRNA expression data and clinical data) of cancer patients in the Data_preprocessing/Origin_Data directory .

#### cd MMOSurv/Data_preprocessing

we summarize the gene expression of patients with the same type of cancer into a csv file in the file extract_RNASeq_expression.py

#### python extract_RNASeq_expression.py

#### python extract_miRNA_expression.py

In the file RNASeq_variable_filter.py, we imputed missing data, filtered out noise-sensitive features (those whose values ​​remained almost constant across samples), and performed a logarithmic transformation.

#### python RNASeq_variable_filter.py

#### python microRNA_variable_filter.py

#### python clinical_preprocess.py

# Run the main routine

git clone https://github.com/LiminLi-xjtu/MMOSurv.git

cd MMOSurv

python mvml_survival_analysis.py
