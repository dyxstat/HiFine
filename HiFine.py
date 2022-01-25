#########The structure of the main script is modified from bin3C########
from Binning_refiner import FindShare,FindShare_bin3C
from Cluster import ClusterBin
from utils import load_object,save_object,make_dir,app_path,gen_bins
from Merge import merge,assign
import logging
import sys
import argparse
import os
import warnings
import pandas as pd

##Ignore the warning information of package deprecation##
warnings.filterwarnings("ignore")

__version__ = '1.0, released at 01/2022'

if __name__ == '__main__':
    
    def mk_version():
        return 'HiFine v{}'.format(__version__)

    def out_name(base, suffix):
        return '{}{}'.format(base, suffix)

    def ifelse(arg, default):
        if arg is None:
            return default
        else:
            return arg

    runtime_defaults = {
        'min_frac': 0.8,
        'alpha': 0.6,
        'beta': 0.3,
        'min_binsize': 500000,
        'min_complete': 500000
    }

    global_parser = argparse.ArgumentParser(add_help=False)
    global_parser.add_argument('-V', '--version', default=False, action='store_true', help='Show the application version')
    global_parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    global_parser.add_argument('--cover', default=False, action='store_true', help='Cover existing files')
    global_parser.add_argument('--log', help='Log file path [OUTDIR/HiFine.log]')


    parser = argparse.ArgumentParser(description='HiFine: integrating Hi-c-based and shotgun-based methods to reFine binning of metagenomic contigs')

    subparsers = parser.add_subparsers(title='commands', dest='command', description='Valid commands',
                                       help='choose an analysis stage for further options')

    cmd_pl = subparsers.add_parser('refine', parents=[global_parser],
                                      description='Refine the original bins from shotgun-based and Hi-C based binning methods.')


    '''
    pipeline subparser input
    '''
    cmd_pl.add_argument('--min-binsize', type=int,
                               help='Minimum bin size of shared bins [default: 500000]')
    cmd_pl.add_argument('--min-complete', type=int,
                               help='Minimum bin size of relatively complete bins [default: 500000]')
    cmd_pl.add_argument('--min-frac', type=float,
                               help='fraction to determine the relatively complete bins [default: 0.8]')
    cmd_pl.add_argument('--alpha', type=float,
                               help='hyperparameter in step2 [default: 0.6]')
    cmd_pl.add_argument('--beta', type=float,
                               help='hyperparameter in step3 [default: 0.3]')
    cmd_pl.add_argument('--bin3C', default=False, action='store_true',
                               help='Whether the Hi-C-based binning method is bin3C')
    cmd_pl.add_argument('HiC', help='Folder of bins constrcuted by Hi-C data')
    cmd_pl.add_argument('Shotgun', help='Folder of bins constrcuted by Shotgun data')
    cmd_pl.add_argument('MAP', help='Contact Map instance [HiCzin_normalized_contact.gz]')
    cmd_pl.add_argument('FASTA', help='Assembled contigs file [.fasta]')
    cmd_pl.add_argument('OUTDIR', help='Output directory')
 

    args = parser.parse_args()

    if args.version:
        print(mk_version())
        sys.exit(0)

    try:
        make_dir(args.OUTDIR, args.cover)
    except IOError:
        print('Error: cannot find out directory or the directory already exists')
        sys.exit(1)

    logging.captureWarnings(True)
    logger = logging.getLogger('main')

    # root log listens to everything
    root = logging.getLogger('')
    root.setLevel(logging.DEBUG)

    # log message format
    formatter = logging.Formatter(fmt='%(levelname)-8s | %(asctime)s | %(name)7s | %(message)s')

    # Runtime console listens to INFO by default
    ch = logging.StreamHandler()
    if args.verbose:
        ch.setLevel(logging.DEBUG)
    else:
        ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    root.addHandler(ch)

    # File log listens to all levels from root
    if args.log is not None:
        log_path = args.log
    else:
        log_path = os.path.join(args.OUTDIR, 'hifine.log')
    fh = logging.FileHandler(log_path, mode='a')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    root.addHandler(fh)

    # Add some environmental details
    logger.debug(mk_version())
    logger.debug(sys.version.replace('\n', ' '))
    logger.debug('Command line: {}'.format(' '.join(sys.argv)))


    if args.command == 'refine':
        cl = load_object(args.MAP)
        map_matrix = cl.seq_map
        contig_name = cl.name
        contig_len = cl.len

        ref_len = {}
        for i in range(contig_name.shape[0]):
            ref_len[contig_name[i]] = contig_len[i]

        ref_ind = {}
        for index , contig in enumerate(contig_name):
            ref_ind[contig] = index
        
        logger.info('Loading Hi-C contact matrix normalized by Hi-C contacts')
        logger.info('Hi-C contact matrix contain {} contigs with total length'.format(map_matrix.shape[0] , sum(ref_len.values())))
        logger.info('Begin Step 1 construct fragmented bins...')
        # Create a contact map for analysis by HiCzin
        if args.bin3C:
            fs = FindShare_bin3C( args.HiC,
                            args.Shotgun,
                            ref_len,
                            min_binsize=ifelse(args.min_binsize, runtime_defaults['min_binsize']),
                            min_frac=ifelse(args.min_frac, runtime_defaults['min_frac']),
                            min_complete_size=ifelse(args.min_complete, runtime_defaults['min_complete']))
        else:
            fs = FindShare( args.HiC,
                            args.Shotgun,
                            ref_len,
                            min_binsize=ifelse(args.min_binsize, runtime_defaults['min_binsize']),
                            min_frac=ifelse(args.min_frac, runtime_defaults['min_frac']),
                            min_complete_size=ifelse(args.min_complete, runtime_defaults['min_complete']))


        ##########Step2 merge fragmented bins#############
        logger.info('Begin Step 2 merge fragmented bins...')
        gp_merge = merge(fs.shared_bin_plus , ref_ind , map_matrix , ifelse(args.alpha, runtime_defaults['alpha']))
    

        ##########Step3 assign short contigs##############
        logger.info('Begin Step 3 recruit unbinned contigs into merged bins...')
        gp_final = assign(gp_merge , map_matrix , ref_ind , ifelse(args.beta , runtime_defaults['beta']))
        
        dist_cluster = {}
        for ci in range(len(gp_final)):
            for contig in gp_final[ci]:
                dist_cluster[contig] = 'group'+str(ci)
        
        with open(os.path.join(args.OUTDIR ,'hifine_cluster.txt'),'w') as out:
            for key , value in dist_cluster.items():
                out.write(str(key)+ '\t' +str(value))
                out.write('\n')
                
        gen_bins(args.FASTA , os.path.join(args.OUTDIR ,'hifine_cluster.txt') , os.path.join(args.OUTDIR , 'HIFINE_BIN') )






