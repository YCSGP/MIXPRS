#!/usr/bin/env python

"""
GWAS subsampling algorithm for MIX.

"""

import numpy as np
from scipy import linalg
from scipy.stats import norm
from numpy import random

def eigen_decomp(ldblk):
    eigenvalues, eigenvectors = linalg.eigh(ldblk)  # Use eigh for symmetric matrices
    
    # Compute sqrt with thresholding
    threshold = 1e-10
    Lambda_transformed = np.sqrt(np.maximum(eigenvalues, threshold))
    
    ldblk_transformed = eigenvectors @ np.diag(Lambda_transformed) @ eigenvectors.T
    
    return ldblk_transformed


def subsample2(snp_dict, beta_std_dict, se_dict, n_dict, frq_dict, idx_dict, ld_blk, blk_size, indep_approx, train_tune_ratio, repeat, pop, chrom, out_dir, out_name, seed):
    print('... MIX GWAS subsampling ...')

    # seed
    if seed != None:
        random.seed(seed)

    # derived stats
    p = len(snp_dict['SNP'])
    n_blk = len(ld_blk)
    train_n_dict = n_dict * train_tune_ratio / (train_tune_ratio + 1)
    tune_n_dict = n_dict / (train_tune_ratio + 1)

    # initialization
    train_beta_std_dict = np.zeros(p)
    tune_beta_std_dict = np.zeros(p)
    train_z_dict = np.zeros(p)
    tune_z_dict = np.zeros(p)
    train_p_dict = np.zeros(p)
    tune_p_dict = np.zeros(p)
    train_se_dict = np.zeros(p)
    tune_se_dict = np.zeros(p)
    train_beta_dict = np.zeros(p)
    tune_beta_dict = np.zeros(p)
    
    # snp information
    snp_pp = [snp_dict['SNP'][ii] for ii in idx_dict]
    chr_pp = [snp_dict['CHR'][ii] for ii in idx_dict]
    bp_pp = [snp_dict['BP'][ii] for ii in idx_dict]
    a1_pp = [snp_dict['A1'][ii] for ii in idx_dict]
    a2_pp = [snp_dict['A2'][ii] for ii in idx_dict]
    frq_pp = [frq_dict[ii] for ii in idx_dict]

    # GWAS subsampling
    for rpt in range(1,repeat+1):
        print('--- repeat-' + str(rpt) + ' ---')

        mm = 0

        for kk in range(n_blk):  # number of LD blocks in that population
            if blk_size[kk] == 0:
                continue
            else:                        
                idx_blk = range(mm,mm+blk_size[kk])

                if indep_approx == 'TRUE':
                    if blk_size[kk] > 1:
                        ld_blk_half = np.eye(blk_size[kk])
                    else:
                        ld_blk_half = 1
                else:
                    if blk_size[kk] > 1:
                        ld_blk_half = eigen_decomp(ld_blk[kk])
                    else:
                        ld_blk_half = np.sqrt(ld_blk[kk])

                z_matrix = random.randn(blk_size[kk])
                train_beta_std_dict[idx_blk] = beta_std_dict[idx_blk] + np.sqrt(tune_n_dict[idx_blk] / (n_dict[idx_blk] * train_n_dict[idx_blk])) * ld_blk_half @ z_matrix
                tune_beta_std_dict[idx_blk] = beta_std_dict[idx_blk] - np.sqrt(train_n_dict[idx_blk] / (n_dict[idx_blk] * tune_n_dict[idx_blk])) * ld_blk_half @ z_matrix

                train_z_dict[idx_blk] = train_beta_std_dict[idx_blk] * np.sqrt(train_n_dict[idx_blk])
                tune_z_dict[idx_blk] = tune_beta_std_dict[idx_blk] * np.sqrt(tune_n_dict[idx_blk])
                train_p_dict[idx_blk] = 2 * (1 - norm.cdf(np.abs(train_z_dict[idx_blk])))
                tune_p_dict[idx_blk] = 2 * (1 - norm.cdf(np.abs(tune_z_dict[idx_blk])))
                train_se_dict[idx_blk] = se_dict[idx_blk] * np.sqrt(n_dict[idx_blk] / train_n_dict[idx_blk])
                tune_se_dict[idx_blk] = se_dict[idx_blk] * np.sqrt(n_dict[idx_blk] / tune_n_dict[idx_blk])
                train_beta_dict[idx_blk] = train_z_dict[idx_blk] * train_se_dict[idx_blk]
                tune_beta_dict[idx_blk] = tune_z_dict[idx_blk] * tune_se_dict[idx_blk]

                mm += blk_size[kk]
        
        # write subsampled GWAS
        eff_file_train = out_dir + '/' + '%s_%s_train_GWAS_approx%s_ratio%.2f_repeat%d.txt' % (out_name, pop, indep_approx, train_tune_ratio, rpt)
        eff_file_tune = out_dir + '/' + '%s_%s_tune_GWAS_approx%s_ratio%.2f_repeat%d.txt' % (out_name, pop, indep_approx, train_tune_ratio, rpt)
        
        with open(eff_file_train, 'a') as ff:
            if chrom == 1:
                ff.write('SNP\tCHR\tBP\tA1\tA2\tA1_Frq\tBETA\tSE\tZ\tP\tN\n')

            for snp, chr, bp, a1, a2, a1_frq, beta, se, z, p, n in zip(snp_pp, chr_pp, bp_pp, a1_pp, a2_pp, frq_pp, train_beta_dict, train_se_dict, train_z_dict, train_p_dict, train_n_dict):
                ff.write('%s\t%d\t%d\t%s\t%s\t%.6e\t%.6e\t%.6e\t%.6e\t%.12e\t%d\n' % (snp, chr, bp, a1, a2, a1_frq, beta, se, z, p, n))

        with open(eff_file_tune, 'a') as ff:
            if chrom == 1:
                ff.write('SNP\tCHR\tBP\tA1\tA2\tA1_Frq\tBETA\tSE\tZ\tP\tN\n')

            for snp, chr, bp, a1, a2, a1_frq, beta, se, z, p, n in zip(snp_pp, chr_pp, bp_pp, a1_pp, a2_pp, frq_pp, tune_beta_dict, tune_se_dict, tune_z_dict, tune_p_dict, tune_n_dict):
                ff.write('%s\t%d\t%d\t%s\t%s\t%.6e\t%.6e\t%.6e\t%.6e\t%.12e\t%d\n' % (snp, chr, bp, a1, a2, a1_frq, beta, se, z, p, n))
