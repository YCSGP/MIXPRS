# MIXPRS

MIXPRS is a data fission-based multi-population PRS integration framework designed to effectively combine PRS derived from multiple populations and methods. The MIXPRS pipeline requires **only GWAS summary statistics and LD reference panels**, and involves three main steps (Figure 1):

* **Step 1: MIX-GWAS subsampling**
  Generate independent subsampled training and tuning GWAS datasets.

* **Step 2: MIX-PRS combining weights**
  Estimate combining weights using two component methods, **JointPRS-auto** and **SDPRX**, each requiring only GWAS summary statistics and LD reference panels.

* **Step 3: Obtain the final MIXPRS**
  Calculate the integrated MIXPRS using weights estimated in Step 2.

Detailed implementations of these component methods are available in their respective repositories:

* [JointPRS-auto](https://github.com/LeqiXu/JointPRS)
* [SDPRX](https://github.com/eldronzhou/SDPRX)

<p align="center">
  <img src="https://github.com/user-attachments/files/20256268/Figure1.pdf" alt="MIXPRS Workflow"/>
</p>

## Getting Started
In this section, we will offer step-by-step guidance on MIXPRS implementation.

### 1. MIXPRS Installation
For the first time, you need to use the following code to install MIXPRS:
```
git clone https://github.com/LeqiXu/MIXPRS.git
cd MIXPRS
conda env create -f environment.yml
conda activate MIXPRS
```

After this, you only need to use
```
conda activate MIXPRS
```

### 2. LD Reference Panel Download
We use reference panels from [PRS-CSx](https://github.com/getian107/PRScsx#getting-started) and you can follow their instructions to download them. It is strongly recommended to create two subfolders within your reference directory:

- **1KG**: This subfolder should contain LD reference panels constructed using the 1000 Genomes Project phase 3 samples and the corresponding SNP information file.
- **UKBB**: This subfolder should contain LD reference panels constructed using the UK Biobank data and the corresponding SNP information file.

Then place the downloaded LD reference panels and the SNP information file into their corresponding subfolders.
