import os

def call_ECC(args):
    import sys
    from collections import defaultdict
    from .logger import get_logger
    from .utils import check_file, check_dir
    from .find_ccs import find_ccs_reads, load_ccs_reads
    from .find_eccdna import scan_ccs_reads
    from .cluster import cluster_all_ret, corboundary_print, merge_subset_tree

    if args.input is None or args.output is None:
        sys.exit('Please provide input and output file.')
    if args.reference is None:
        sys.exit('Please specific FASTA of reference genome')

    # Check parameters
    in_file = check_file(args.input)
    out_dir = check_dir(args.output)
    ref_fasta = check_file(args.reference)

    check_dir(out_dir + '/tmp')
    prefix = args.prefix
    threads = int(args.threads)
    mapq = int(args.mapq)


    logger = get_logger('ECCTRY', fname='{}/{}.log'.format(out_dir, prefix))
    logger.info('----------------- Input paramters ------------------')
    logger.info('Input reads: ' + os.path.basename(in_file))
    logger.info('Output directory: ' + os.path.basename(out_dir))
    logger.info('Multi threads: {}'.format(args.threads))
    logger.info('----------------- Calling CCSs -----------------')

    reads_count = defaultdict(int)
    



    # Calling CCSs
    if os.path.exists('{}/tmp/{}.ccs.fa'.format(out_dir, prefix)) and os.path.exists('{}/tmp/{}.raw.fa'.format(out_dir, prefix)):
        logger.info('Step 1 - Loading CCSs in previous run')
        ccs_seq = load_ccs_reads(out_dir, prefix)
        reads_count['consensus'] = len(ccs_seq)
    else:
        logger.info('Step 1 - Calling CCSs')
        total_reads, ro_reads, ccs_seq = find_ccs_reads(in_file, out_dir, prefix, threads)
        reads_count['total'] = total_reads
        reads_count['consensus'] = ro_reads

    if 'total' in reads_count:
        logger.info('Total Reads: {}'.format(reads_count['total']))
    logger.info('Cyclic Consensus Reads: {}'.format(reads_count['consensus']))

    logger.info('----------------- Calling CCSs finish!-----------------')




    logger.info('Step 2 - Map CCSs to the reference genome')

    reads_count, all_ret = scan_ccs_reads(ccs_seq, ref_fasta, out_dir, prefix, threads, mapq)

    logger.info('----------------- Mapping finish!-----------------')




    logger.info('Step 3 - Cluster eccDNAs based on genomic positions')

    final_tree = cluster_all_ret(all_ret, threads)

    mergesub_tree = merge_subset_tree(final_tree)
    corboundary_print(mergesub_tree, out_dir, prefix)

    logger.info('----------------- Clustering finish!-----------------')

    return



def main():
    import argparse

    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(help='commands')

    CCS_parser = subparsers.add_parser('call_ECC')

    CCS_parser.add_argument('-i', '--in', dest='input', metavar='READS', default=None,help='Input reads.fq.gz', )
    CCS_parser.add_argument('-o', '--out', dest='output', metavar='DIR', default=None,help='Output directory, default: ./', )    
    CCS_parser.add_argument('-p', '--prefix', dest='prefix', metavar='PREFIX', default="ECCTRY",help='Output sample prefix, (default: %(default)s)', )
    CCS_parser.add_argument('-t', '--threads', dest='threads', metavar='INT', default=os.cpu_count(),help='Number of threads, (default: use all cores)', )
    CCS_parser.add_argument('-r', '--ref', dest='reference', metavar='REF', default=None,help='Reference genome FASTA file', )
    # 先不使用注释文件或者现有的数据库
    # CCS_parser.add_argument('-a', '--anno', dest='gtf', metavar='GTF', default=None,help='Genome reference gtf, (optional)', )

    CCS_parser.add_argument('-q', '--mapq', dest='mapq', metavar='MAPQ', default=30, help='Minimum mapping quality allowed in the alignments', )

    CCS_parser.set_defaults(func=call_ECC)


    args = parser.parse_args()

    # Run function
    try:
        func = args.func
    except AttributeError:
        parser.error("too few arguments")
    func(args)



if __name__ == '__main__':
    main()

    
