#!/usr/bin/env python

"""
Parse the reference panel, summary statistics, and prune set.

"""


import os
import scipy as sp
from scipy.stats import norm
from scipy import linalg
import numpy as np
import h5py

def parse_ref(ref_file, chrom, ref):
    print('... parse reference file: %s ...' % ref_file)

    if ref == '1kg' or ref == 'ukbb':
        ref_dict = {'CHR':[], 'SNP':[], 'BP':[], 'A1':[], 'A2':[], 
                    'FRQ_AFR':[], 'FRQ_AMR':[], 'FRQ_EAS':[], 'FRQ_EUR':[], 'FRQ_SAS':[],
                    'FLP_AFR':[], 'FLP_AMR':[], 'FLP_EAS':[], 'FLP_EUR':[], 'FLP_SAS':[]}
        with open(ref_file) as ff:
            header = next(ff)
            for line in ff:
                ll = (line.strip()).split()
                if int(ll[0]) == chrom:
                    ref_dict['CHR'].append(chrom)
                    ref_dict['SNP'].append(ll[1])
                    ref_dict['BP'].append(int(ll[2]))
                    ref_dict['A1'].append(ll[3])
                    ref_dict['A2'].append(ll[4])
                    ref_dict['FRQ_AFR'].append(float(ll[5]))
                    ref_dict['FRQ_AMR'].append(float(ll[6]))
                    ref_dict['FRQ_EAS'].append(float(ll[7]))
                    ref_dict['FRQ_EUR'].append(float(ll[8]))
                    ref_dict['FRQ_SAS'].append(float(ll[9]))
                    ref_dict['FLP_AFR'].append(int(ll[10]))
                    ref_dict['FLP_AMR'].append(int(ll[11]))
                    ref_dict['FLP_EAS'].append(int(ll[12]))
                    ref_dict['FLP_EUR'].append(int(ll[13]))
                    ref_dict['FLP_SAS'].append(int(ll[14]))

    print('... %d SNPs on chromosome %d read from %s ...' % (len(ref_dict['SNP']), chrom, ref_file))
    return ref_dict


def parse_prune_snplist(prune_snplist):
    print('... parse prune snplist file: %s ...' % prune_snplist)

    prune_dict = {'SNP':[]}
    with open(prune_snplist) as ff:
        for line in ff:
            ll = (line.strip()).split()
            prune_dict['SNP'].append(ll[0])

    print('... %d SNPs read from %s ...' % (len(prune_dict['SNP']), prune_snplist))
    return prune_dict


def parse_sumstats(ref_dict, prune_dict, sst_file, pop):
    print('... parse ' + pop.upper() + ' sumstats file: %s ...' % sst_file)

    ATGC = ['A', 'T', 'G', 'C']
    sst_dict = {'SNP':[], 'A1':[], 'A2':[]}
    with open(sst_file) as ff:
        header = next(ff)
        for line in ff:
            ll = (line.strip()).split()
            if ll[1] in ATGC and ll[2] in ATGC:
                sst_dict['SNP'].append(ll[0])
                sst_dict['A1'].append(ll[1])
                sst_dict['A2'].append(ll[2])

    print('... %d SNPs read from %s ...' % (len(sst_dict['SNP']), sst_file))


    idx = [ii for (ii,frq) in enumerate(ref_dict['FRQ_'+pop.upper()]) if frq>0]
    snp_ref = [ref_dict['SNP'][ii] for ii in idx]
    a1_ref = [ref_dict['A1'][ii] for ii in idx]
    a2_ref = [ref_dict['A2'][ii] for ii in idx]

    mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    ref_snp = set(zip(snp_ref, a1_ref, a2_ref)) | set(zip(snp_ref, a2_ref, a1_ref)) | \
              set(zip(snp_ref, [mapping[aa] for aa in a1_ref], [mapping[aa] for aa in a2_ref])) | \
              set(zip(snp_ref, [mapping[aa] for aa in a2_ref], [mapping[aa] for aa in a1_ref]))

    sst_snp = set(zip(sst_dict['SNP'], sst_dict['A1'], sst_dict['A2']))

    comm_snp = ref_snp & sst_snp

    if prune_dict != None:
        prune_snp = set(prune_dict['SNP'])
        comm_snp = {entry for entry in comm_snp if entry[0] in prune_snp}

    print('... %d common SNPs in the %s reference, %s sumstats, and provided prune set ...' % (len(comm_snp), pop.upper(), pop.upper()))


    sst_eff = {}
    with open(sst_file) as ff:
        header = (next(ff).strip()).split()
        header = [col.upper() for col in header]
        for line in ff:
            ll = (line.strip()).split()
            snp = ll[0]; a1 = ll[1]; a2 = ll[2]; n_subj = int(float(ll[7])); n_sqrt = sp.sqrt(n_subj)
            if a1 not in ATGC or a2 not in ATGC:
                continue
            if (snp, a1, a2) in comm_snp or (snp, mapping[a1], mapping[a2]) in comm_snp:
                if 'BETA' in header:
                    beta = float(ll[3]); se = float(ll[4])
                elif 'OR' in header:
                    beta = sp.log(float(ll[3])); se = float(ll[4])

                p = max(float(ll[6]), 1e-323)
                beta_std = sp.sign(beta)*abs(norm.ppf(p/2.0))/n_sqrt
                sst_eff.update({snp: (beta_std, a1, a2, se, n_subj)})
            elif (snp, a2, a1) in comm_snp or (snp, mapping[a2], mapping[a1]) in comm_snp:
                if 'BETA' in header:
                    beta = float(ll[3])
                elif 'OR' in header:
                    beta = sp.log(float(ll[3]))

                p = max(float(ll[6]), 1e-323)
                beta_std = -1*sp.sign(beta)*abs(norm.ppf(p/2.0))/n_sqrt
                sst_eff.update({snp: (beta_std, a1, a2, se, n_subj)})


    sst_dict = {'SNP':[], 'A1':[], 'A2':[], 'BETA_STD':[], 'SE':[], 'N':[], 'FRQ':[], 'FLP':[]}
    for (ii,snp) in enumerate(ref_dict['SNP']):
        if snp in sst_eff:
            beta_std, a1, a2, se, n_subj = sst_eff[snp]
            sst_dict['SNP'].append(snp)
            sst_dict['A1'].append(a1)
            sst_dict['A2'].append(a2)
            sst_dict['BETA_STD'].append(beta_std)
            sst_dict['SE'].append(se)
            sst_dict['N'].append(n_subj)

            a1_ref = ref_dict['A1'][ii]; a2_ref = ref_dict['A2'][ii]
            if (snp, a1_ref, a2_ref) in comm_snp or (snp, mapping[a1_ref], mapping[a2_ref]) in comm_snp:
                sst_dict['FRQ'].append(ref_dict['FRQ_'+pop.upper()][ii])
                sst_dict['FLP'].append(ref_dict['FLP_'+pop.upper()][ii])
            elif (snp, a2_ref, a1_ref) in comm_snp or (snp, mapping[a2_ref], mapping[a1_ref]) in comm_snp:
                sst_dict['FRQ'].append(1-ref_dict['FRQ_'+pop.upper()][ii])
                sst_dict['FLP'].append(-1*ref_dict['FLP_'+pop.upper()][ii])

    return sst_dict


def parse_prs_beta(sst_dict, prs_beta_file, pp):
    print('... parse prs beta file: %s ...' % prs_beta_file)

    sst_dict["PRS_BETA_STD" + pp] = [0] * len(next(iter(sst_dict.values())))

    with open(prs_beta_file) as ff:
        header = (next(ff).strip()).split()
        header = [col.upper() for col in header]

        sst_snp_index_map = {snp: idx for idx, snp in enumerate(sst_dict['SNP'])}

        for line in ff:
            ll = (line.strip()).split()
            snp = ll[0]; a1 = ll[1]; beta = float(ll[2])
            
            if snp in sst_snp_index_map:
                sst_idx = sst_snp_index_map[snp]

                sst_a1 = sst_dict['A1'][sst_idx]
                sst_a2 = sst_dict['A2'][sst_idx]
                sst_frq = sst_dict['FRQ'][sst_idx]
                sst_het = np.sqrt(2.0 * sst_frq * (1.0 - sst_frq))

                if a1 in sst_a1:
                    beta_std = beta * sst_het
                    sst_dict["PRS_BETA_STD" + pp][sst_idx] = beta_std
                elif a1 in sst_a2:
                    beta_std = -1 * beta * sst_het
                    sst_dict["PRS_BETA_STD" + pp][sst_idx] = beta_std

    return sst_dict


def parse_ldblk(ldblk_dir, sst_dict, pop, chrom, ref):
    print('... parse %s reference LD on chromosome %d ...' % (pop.upper(), chrom))

    if ref == '1kg':
        chr_name = ldblk_dir + '/ldblk_1kg_' + pop.lower() + '/ldblk_1kg_chr' + str(chrom) + '.hdf5'
    elif ref == 'ukbb':
        chr_name = ldblk_dir + '/ldblk_ukbb_' + pop.lower() + '/ldblk_ukbb_chr' + str(chrom) + '.hdf5'

    hdf_chr = h5py.File(chr_name, 'r')
    n_blk = len(hdf_chr)
    ld_blk = [sp.array(hdf_chr['blk_'+str(blk)]['ldblk']) for blk in range(1,n_blk+1)]

    snp_blk = []
    for blk in range(1,n_blk+1):
         snp_blk.append([bb.decode("UTF-8") for bb in list(hdf_chr['blk_'+str(blk)]['snplist'])])

    blk_size = []
    mm = 0
    for blk in range(n_blk):
        idx = [ii for (ii,snp) in enumerate(snp_blk[blk]) if snp in sst_dict['SNP']]
        blk_size.append(len(idx))
        if idx != []:
            idx_blk = range(mm,mm+len(idx))
            flip = [sst_dict['FLP'][jj] for jj in idx_blk]
            ld_blk[blk] = ld_blk[blk][sp.ix_(idx,idx)]*sp.outer(flip,flip)

            _, s, v = linalg.svd(ld_blk[blk])
            h = sp.dot(v.T, sp.dot(sp.diag(s), v))
            ld_blk[blk] = (ld_blk[blk]+h)/2

            mm += len(idx)
        else:
            ld_blk[blk] = sp.array([])

    return ld_blk, blk_size


def align_ldblk(ref_dict, sst_dict, chrom):
    print('... align reference LD on chromosome %d ...' % chrom)

    snp_dict = {'CHR':[], 'SNP':[], 'BP':[], 'A1':[], 'A2':[]}
    for (ii,snp) in enumerate(ref_dict['SNP']):
        if snp in sst_dict['SNP']:
            snp_dict['SNP'].append(snp)
            snp_dict['CHR'].append(ref_dict['CHR'][ii])
            snp_dict['BP'].append(ref_dict['BP'][ii])

            idx = sst_dict['SNP'].index(snp)
            snp_dict['A1'].append(sst_dict['A1'][idx])
            snp_dict['A2'].append(sst_dict['A2'][idx])

    n_snp = len(snp_dict['SNP'])
    print('... %d valid SNPs ...' % n_snp)

    beta_std_dict = np.array(sst_dict['BETA_STD'])
    se_dict = np.array(sst_dict['SE'])
    n_dict = np.array(sst_dict['N'])
    frq_dict = np.array(sst_dict['FRQ'])
    idx_dict = [ii for (ii, snp) in enumerate(snp_dict['SNP']) if snp in sst_dict['SNP']]

    return snp_dict, beta_std_dict, se_dict, n_dict, frq_dict, idx_dict


def align_prs_beta(sst_dict, n_prs_beta):
    prs_beta_matrix = []

    for pp in range(n_prs_beta):
        column_name = "PRS_BETA_STD" + str(pp)
        prs_beta_matrix.append(sst_dict[column_name])

    prs_beta_matrix = np.array(prs_beta_matrix).T

    return prs_beta_matrix


def standardize_prs_beta_for_each_method(ld_blk, blk_size, prs_beta_matrix, n_prs_beta, indep_approx):
    total_per_beta_ldblk_beta = np.zeros(n_prs_beta)

    for prs_num in range(n_prs_beta):

        sum_for_this_prs = 0.0

        for chrom in range(1, 23):
            n_blk = len(ld_blk[chrom])
            mm = 0 

            for kk in range(n_blk):
                if blk_size[chrom][kk] == 0:
                    continue

                idx_blk = np.arange(mm, mm + blk_size[chrom][kk])

                if indep_approx == 'TRUE':
                    block_value = prs_beta_matrix[chrom][idx_blk, prs_num].T @ prs_beta_matrix[chrom][idx_blk, prs_num]
                else:
                    block_ld = ld_blk[chrom][kk]
                    block_value = prs_beta_matrix[chrom][idx_blk, prs_num].T @ block_ld @ prs_beta_matrix[chrom][idx_blk, prs_num]

                sum_for_this_prs += block_value

                mm += blk_size[chrom][kk]

        total_per_beta_ldblk_beta[prs_num] = sum_for_this_prs

    denom = np.sqrt(total_per_beta_ldblk_beta)

    prs_beta_matrix_scaled = {}
    for chrom in range(1, 23):
        prs_beta_matrix_scaled[chrom] = prs_beta_matrix[chrom] / denom

    return prs_beta_matrix_scaled
