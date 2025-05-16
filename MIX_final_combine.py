#!/usr/bin/env python

import os
import sys
import getopt

import parse_genet
import prs_combine

import numpy as np

def parse_param():
    long_opts_list = ['ref_dir=', 'sst_file=', 'pop=', 'prs_beta_file=', 'weight_file=', 'indep_approx=', 'out_dir=', 'out_name=', 'help']

    param_dict = {'ref_dir': None, 'sst_file': None, 'pop': None, 
                  'prs_beta_file': None, 'weight_file': None, 'indep_approx': 'TRUE',
                  'out_dir': None, 'out_name': None}

    print('\n')

    if len(sys.argv) > 1:
        try:
            opts, args = getopt.getopt(sys.argv[1:], "h", long_opts_list)          
        except:
            print('* Option not recognized.')
            print('* Use --help for usage information.\n')
            sys.exit(2)

        for opt, arg in opts:
            if opt == "-h" or opt == "--help":
                print(__doc__)
                sys.exit(0)
            elif opt == "--ref_dir": param_dict['ref_dir'] = arg
            elif opt == "--sst_file": param_dict['sst_file'] = arg
            elif opt == "--pop": param_dict['pop'] = arg
            elif opt == "--prs_beta_file": param_dict['prs_beta_file'] = arg.split(',')
            elif opt == "--weight_file": param_dict['weight_file'] = arg.split(',')
            elif opt == "--indep_approx": param_dict['indep_approx'] = arg.upper()
            elif opt == "--out_dir": param_dict['out_dir'] = arg
            elif opt == "--out_name": param_dict['out_name'] = arg
    else:
        print(__doc__)
        sys.exit(0)


    if param_dict['ref_dir'] == None:
        print('* Please specify the directory to the reference panel using --ref_dir\n')
        sys.exit(2)
    elif param_dict['sst_file'] == None:
        print('* Please provide at least one summary statistics file using --sst_file\n')
        sys.exit(2)
    elif param_dict['pop'] == None:
        print('* Please specify the population of the GWAS sample using --pop\n')
        sys.exit(2)
    elif param_dict['prs_beta_file'] == None:
        print('* Please provide at least two prs beta file using --prs_beta_file\n')
        sys.exit(2)
    elif param_dict['weight_file'] == None:
        print('* Please provide the weight file using --weight_file\n')
        sys.exit(2)
    elif param_dict['out_dir'] == None:
        print('* Please specify the output directory using --out_dir\n')
        sys.exit(2)
    elif param_dict['out_name'] == None:
        print('* Please specify the prefix of the output file using --out_name\n')
        sys.exit(2)

    for key in param_dict:
        print('--%s=%s' % (key, param_dict[key]))

    print('\n')
    return param_dict


def main():
    param_dict = parse_param()
    n_prs_beta = len(param_dict['prs_beta_file'])

    combine_weight = prs_combine.combine_weight(param_dict['weight_file'])
    combine_weight = combine_weight.reshape(-1,1)
    
    if n_prs_beta != combine_weight.shape[0]:
        print('* Please ensure the number of PRS beta files matches the entries in the weight file. Use --prs_beta_file and --weight_file to provide them\n')
        sys.exit(2)
    else:
        print('*** %d prs beta files detected ***\n' % n_prs_beta)

    ld_blk = {}
    blk_size = {}
    snp_dict = {}
    frq_dict = {}
    prs_beta_matrix = {}

    for chrom in range(1,23):
        print('##### process chromosome %d #####' % int(chrom))

        if os.path.isfile(param_dict['ref_dir'] + '/snpinfo_mult_1kg_hm3'):
            ref = '1kg'
            ref_dict = parse_genet.parse_ref(param_dict['ref_dir'] + '/snpinfo_mult_1kg_hm3', int(chrom), ref)
        elif os.path.isfile(param_dict['ref_dir'] + '/snpinfo_mult_ukbb_hm3'):
            ref = 'ukbb'
            ref_dict = parse_genet.parse_ref(param_dict['ref_dir'] + '/snpinfo_mult_ukbb_hm3', int(chrom), ref)
        
        prune_dict = None
        sst_dict = parse_genet.parse_sumstats(ref_dict, prune_dict, param_dict['sst_file'], param_dict['pop'])

        for pp in range(n_prs_beta):
            if pp == 0:
                sst_prs_dict = parse_genet.parse_prs_beta(sst_dict, param_dict['prs_beta_file'][pp], str(pp))
            else:
                sst_prs_dict = parse_genet.parse_prs_beta(sst_prs_dict, param_dict['prs_beta_file'][pp], str(pp))

        ld_blk[chrom], blk_size[chrom] = parse_genet.parse_ldblk(param_dict['ref_dir'], sst_prs_dict, param_dict['pop'], int(chrom), ref)

        snp_dict[chrom], _, _, _, frq_dict[chrom], _ = parse_genet.align_ldblk(ref_dict, sst_prs_dict, int(chrom))
        prs_beta_matrix[chrom] = parse_genet.align_prs_beta(sst_prs_dict, n_prs_beta)

        print('##### finish chromosome %d #####' % int(chrom))
        print('\n')
    
    print('##### process across 22 chromosomes #####')

    prs_beta_matrix = parse_genet.standardize_prs_beta_for_each_method(ld_blk, blk_size, prs_beta_matrix, n_prs_beta, param_dict['indep_approx'])
    prs_combine.convert_to_per_allele(snp_dict, frq_dict, prs_beta_matrix, combine_weight, n_prs_beta, param_dict['pop'], param_dict['out_dir'], param_dict['out_name'])

    print('##### finish across 22 chromosomes #####')
    print('\n')

if __name__ == '__main__':
    main()
