#!/usr/bin/env python

import os
import sys
import getopt

import parse_genet
import gwas_subsample

def parse_param():
    long_opts_list = ['ref_dir=', 'sst_file=', 'pop=', 'prune_snplist=', 'indep_approx=', 'train_tune_ratio=', 'repeat=', 'out_dir=', 'out_name=', 'seed=', 'help']

    param_dict = {'ref_dir': None, 'sst_file': None, 'pop': None,
                  'prune_snplist': None, 'indep_approx': 'TRUE', 'train_tune_ratio': 3, 'repeat': 4, 
                  'out_dir': None, 'out_name': None, 'seed': None}

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
            elif opt == "--prune_snplist": param_dict['prune_snplist'] = arg
            elif opt == "--indep_approx": param_dict['indep_approx'] = arg.upper()
            elif opt == "--train_tune_ratio": param_dict['train_tune_ratio'] = float(arg)
            elif opt == "--repeat": param_dict['repeat'] = int(arg)
            elif opt == "--out_dir": param_dict['out_dir'] = arg
            elif opt == "--out_name": param_dict['out_name'] = arg
            elif opt == "--seed": param_dict['seed'] = int(arg)
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

    for chrom in range(1,23):
        print('##### process chromosome %d #####' % int(chrom))

        if os.path.isfile(param_dict['ref_dir'] + '/snpinfo_mult_1kg_hm3'):
            ref = '1kg'
            ref_dict = parse_genet.parse_ref(param_dict['ref_dir'] + '/snpinfo_mult_1kg_hm3', int(chrom), ref)
        elif os.path.isfile(param_dict['ref_dir'] + '/snpinfo_mult_ukbb_hm3'):
            ref = 'ukbb'
            ref_dict = parse_genet.parse_ref(param_dict['ref_dir'] + '/snpinfo_mult_ukbb_hm3', int(chrom), ref)
        
        if param_dict['prune_snplist'] != None:
            prune_dict = parse_genet.parse_prune_snplist(param_dict['prune_snplist'])
        else:
            prune_dict = None

        sst_dict = parse_genet.parse_sumstats(ref_dict, prune_dict, param_dict['sst_file'], param_dict['pop'])

        ld_blk, blk_size = parse_genet.parse_ldblk(param_dict['ref_dir'], sst_dict, param_dict['pop'], int(chrom), ref)

        snp_dict, beta_std_dict, se_dict, n_dict, frq_dict, idx_dict = parse_genet.align_ldblk(ref_dict, sst_dict, int(chrom))

        gwas_subsample.subsample2(snp_dict, beta_std_dict, se_dict, n_dict, frq_dict, idx_dict, ld_blk, blk_size,
            param_dict['indep_approx'], param_dict['train_tune_ratio'], param_dict['repeat'], 
            param_dict['pop'], int(chrom), param_dict['out_dir'], param_dict['out_name'], param_dict['seed'])

        print('##### finish chromosome %d #####' % int(chrom))
        print('\n')


if __name__ == '__main__':
    main()


