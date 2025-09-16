# ImmuCellAI 2.0
ImmuCellAI 2.0 is an advanced bulk transcriptome deconvolution tool that accurately infers the proportions of a comprehensive range of immune cell types, including rare subsets, from bulk RNA-seq data.

# ImmuCellAI 2.0 Website
The ImmuCellAI 2.0 Website is a valuable online resource for immunology research, supporting the estimation of immune cell abundance across 9 major immune cell types and their corresponding 53 minor subtypes.
Access it here: https://guolab.wchscu.cn/immucellai/

# Install ImmuCellAI 2.0
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
