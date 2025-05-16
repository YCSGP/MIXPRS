# MIXPRS
MIXPRS is a data fission-based multi-population PRS integration framework that effectively combines PRS derived from multiple populations and methods. The MIXPRS pipeline involves three main steps (Figure 1) and requires only GWAS summary statistics:

(1) MIX-GWAS subsampling to generate independent training and tuning GWAS datasets.
(2) Estimation of MIX-PRS weights via two component methods, **JointPRS-auto** and **SDPRX**, each requiring only GWAS summary statistics and LD reference panels.
(3) Calculation of the final integrated MIXPRS.

Detailed implementations of these two component methods are available at their respective repositories:

* [JointPRS-auto](https://github.com/LeqiXu/JointPRS)
* [SDPRX](https://github.com/eldronzhou/SDPRX)

![MIXPRS Workflow](https://github.com/user-attachments/files/20256268/Figure1.pdf)
