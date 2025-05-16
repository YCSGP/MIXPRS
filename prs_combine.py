#!/usr/bin/env python

"""
PRS weight.

"""

import numpy as np

def combine_weight(weight_file):
    all_weights = []
    for wfile in weight_file:
        with open(wfile, 'r') as f:
            line = f.readline().strip()       # read the single line
            vals = line.split()              # split on whitespace
            float_vals = [float(v) for v in vals]
            all_weights.append(float_vals)
            
    # Convert to numpy array for easy column-wise averaging
    all_weights_arr = np.array(all_weights)  # shape: (num_files, k)
    avg_weights = np.mean(all_weights_arr, axis=0)

    return avg_weights


import numpy as np

def convert_to_per_allele(snp_dict, frq_dict, prs_beta_matrix, combine_weight, n_prs_beta, pop, out_dir, out_name): 
    snp_all = []
    a1_all = []
    frq_all = []

    combine_prs = []
    separate_prs = []

    for chrom in range(1, 23):
        snp_all.extend(snp_dict[chrom]['SNP'])
        a1_all.extend(snp_dict[chrom]['A1'])
        frq_all.extend(frq_dict[chrom])

        frq_chr = np.array(frq_dict[chrom])
        het_chr = np.sqrt(2.0 * frq_chr * (1.0 - frq_chr))

        combine_prs_chr = prs_beta_matrix[chrom] @ combine_weight
        combine_prs_chr = combine_prs_chr.flatten()
        combine_prs_chr = combine_prs_chr / het_chr
        combine_prs.extend(combine_prs_chr)

        sep_prs_chr = prs_beta_matrix[chrom] * combine_weight.reshape(1,-1)
        sep_prs_chr = sep_prs_chr / het_chr[:, None]
        separate_prs.append(sep_prs_chr)

    # 1) Write out single combined file
    combine_prs_file = out_dir + '/' + '%s_%s_MIXPRS.txt' % (out_name, pop)
    combine_header_cols = ["SNP", "A1", "MIXPRS"]

    with open(combine_prs_file, 'w') as ff:
        ff.write("\t".join(combine_header_cols) + "\n")
        for snp, a1, cprs in zip(snp_all, a1_all, combine_prs):
            ff.write('%s\t%s\t%.6e\n' % (snp, a1, cprs))

    # 2) Stack separate PRS arrays and write them out
    separate_prs = np.vstack(separate_prs)
    separate_prs_file = out_dir + '/' + '%s_%s_MIXPRS_separate.txt' % (out_name, pop)
    separate_header_cols = ["SNP", "A1"] + [f"MIXPRS_separate_{i+1}" for i in range(n_prs_beta)]

    with open(separate_prs_file, 'w') as ff:
        ff.write("\t".join(separate_header_cols) + "\n")
        for idx, (snp, a1) in enumerate(zip(snp_all, a1_all)):
            sprs_row = separate_prs[idx]
            sprs = "\t".join(f"{val:.6e}" for val in sprs_row)
            ff.write('%s\t%s\t%s\n' % (snp, a1, sprs))
