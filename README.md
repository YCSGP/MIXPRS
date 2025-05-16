# MIXPRS

MIXPRS is a data fission-based multi-population PRS integration framework designed to effectively combine PRS derived from multiple populations and methods. The MIXPRS pipeline requires **only GWAS summary statistics and LD reference panels**, and involves three main steps (Figure 1):

* **Step1: MIX-GWAS subsampling**
  Generate independent subsampled training and tuning GWAS datasets.

* **Step2: MIX-PRS combining weights**
  Estimate combining weights using two component methods, **JointPRS-auto** and **SDPRX**, each requiring only GWAS summary statistics and LD reference panels.

* **Step3: Obtain MIXPRS**
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

### 3. Pruned SNP List Preparation
In the `MIXPRS` repository, we provide precomputed pruned SNP lists for five populations (EUR, EAS, AFR, SAS, and AMR), derived from genotype data of the Phase-3 1000 Genome Project. These SNP lists were generated using PLINK2 with parameters: window size = 250, step size = 5, and correlation threshold $r^2 = 0.5$.

To replicate our results or generate customized pruned SNP lists, use the following PLINK2 commands:
```bash
module load PLINK/2

i=1
pval=1
r2=0.5
wc=250

for pop in EUR EAS AFR SAS AMR; do

plink2 --bfile ./1000g_phase3_data/geno_data/${pop} \
       --indep-pairwise ${wc} 5 ${r2} \
       --out ./snplist/${pop}_prune_pval${pval}_r2${r2}_wc${wc}_${i}
done
```
You can adjust these parameters or apply the commands to other genotype datasets according to your analytical needs.

### 4. Summary Statistics Preparation
We require the following format for summary statistics input (including the header line):
```
SNP         A1  A2  BETA        SE          Z                P      N
rs3934834   T   C   0.00492590  0.00761858  0.646564057869   0.518  100061
rs3766192   C   T   0.00144597  0.00580197  0.249220523374   0.803  100953
rs9442372   A   G   0.00261076  0.00586576  0.445084694907   0.656  98059
...
```

Here
- `SNP`: SNP rsID.
-  `A1`: Effect allele.
-  `A2`: Alternative allele.
-  `BETA`: Effect size of allele A1, which is only used to determine the direction of an association.
-  `SE`: Standard error of the effect size.
-  `Z`: Z-score (BETA divided by SE).
-  `P`: P-value of the effect, which is used to calculate the standardized effect size.
-  `N`: GWAS sample size.

### 5. MIXPRS Implementation
Below is a step-by-step guide to implementing MIXPRS:

#### Step 0: MIXPRS data preparation

Before running MIXPRS, ensure you have prepared the following datasets (see previous sections for detailed information):
* **LD Reference Panel** (2. LD Reference Panel Download)
* **Pruned SNP List** (3. Pruned SNP List Preparation)
* **GWAS Summary Statistics** (4. Summary Statistics Preparation)

In the following example commands, we assume:
* You are using the **1KG LD reference panel** downloaded as recommended.
* The precomputed pruned SNP list provided within the cloned MIXPRS repository (`snplist` folder) is used.
* GWAS summary statistics are formatted as `${trait}_${pop}_MIXPRS_sumstat.txt`, following the required format specified previously. We provide an example summary statistics dataset for EAS HDL (500 SNPs on chromosome 1) within the cloned MIXPRS repository (`example_data` folder), obtained from [GLGC](https://csg.sph.umich.edu/willer/public/glgc-lipids2021/).

#### Step1: MIX-GWAS subsampling
This step partitions a single original GWAS dataset into independent training and tuning GWAS datasets using data fission principles.

The default parameters for subsampling are:
* Requires a pruned SNP list (`--prune_snplist`) provided in the MIXPRS repository (`snplist` folder).
* `indep_approx=TRUE`: uses an identity covariance matrix for the pruned SNPs instead of LD reference panels (the `--ref_dir` flag is still required).
* `train_tune_ratio=3`: divides the original GWAS dataset into training (3/4 of samples) and tuning (1/4 of samples) subsets.
* `repeat=4`: generates four independent training-tuning GWAS summary statistics pairs for robust evaluation.

Below is an example command illustrating this step. Replace `${MIXPRS_path}`, `${ref_data_path}`, `${summary_stats}`, `${pop}`, `${trait}`, and `${output_path}` with your actual paths and filenames:
```bash
conda activate MIXPRS

python ${MIXPRS_path}/MIX_subsample2.py \
  --ref_dir=${ref_data_path}/1KG \
  --sst_file=${summary_stats}/${trait}_${pop}_MIXPRS_sumstat.txt \
  --pop=${pop} \
  --prune_snplist=${MIXPRS_path}/snplist/${pop}_prune_pval1_r20.5_wc250_1.snplist \
  --indep_approx=TRUE \
  --train_tune_ratio=3 \
  --repeat=4 \
  --out_dir=${output_path}/subsample/clean \
  --out_name=${trait}_prune_snplist_1
```
Ensure all placeholders match your actual directory structure and filenames.


#### Step2: MIX-PRS combining weights

#### Step3: Obtain MIXPRS
