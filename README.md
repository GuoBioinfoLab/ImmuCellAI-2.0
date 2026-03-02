# ImmuCellAI 2.0
ImmuCellAI 2.0 is an advanced bulk transcriptome deconvolution tool for **immunogenomic research**, optimized to accurately infer the relative proportions of immune cell populations from bulk RNA-seq data. It supports the quantification of **9 major immune cell types and 53 minor functional subtypes** (including rare subsets like exhausted B/T cells, tissue-resident memory T cells, and tumor-associated macrophages).

# ImmuCellAI 2.0 Website
The ImmuCellAI 2.0 Website is a valuable online resource for immunology research, supporting the estimation of immune cell abundance across 9 major immune cell types and their corresponding 53 minor subtypes.
Access it here: https://guolab.wchscu.cn/immucellai/

# ImmuCellAI 2.0 Python Package
## Install using pip
pip install immucellai2

# Usage
ImmuCellAI 2.0 expects a TPM matrix as input and can be implemented with only two lines of code in Python.
```
import immucellai2
reference_data = immucellai2.load_tumor_reference_data()
```
```
result = immucellai2.run_ImmuCellAI2(
    reference_file=reference_data,
    sample_file=<file_path>,    #  User-defined
    output_file=<file_path>,    #  User-defined
    thread_num=8
)
```
# Citation
