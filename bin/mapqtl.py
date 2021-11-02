#!/usr/bin/env python3

## CIS eQLT mapping with tensorQTL
## see https://github.com/broadinstitute/tensorqtl

import sys
import os
import argparse
from datetime import datetime
import torch
import numpy
import pandas
import tensorqtl

sys.path.insert(1, os.path.dirname(__file__))
from post import calculate_qvalues

SEED = 123456

def set_argument_parser():
    parser = argparse.ArgumentParser(description="Map eQTL with linear models in tensorQTL.")

    parser.add_argument("--output-dir", "-o",
        default=os.curdir,
        help="output directory",
        dest="output_dir")
    parser.add_argument('--output-prefix', default = "test_output",
                        help='Prefix for output file names',
                        dest = "oufn_prfx")
    parser.add_argument("--gene-expression-bed", "-e", required = True,
                        help="Normalized gene expression matrix in gzipped BED format.")
    parser.add_argument("--genotype-plink-prefix", "-g", required = True,
                        help="PLINK file prefix for genotypes [<prefix>.bed, <prefix>.bim, <prefix>.fam].",
                        dest = "plink_gt_fnprfx")
    parser.add_argument("--genotype-pcs-plink", "-p", required = True,
                        help="Principal components of genotypes (samples x PCs) [TSV].")
    parser.add_argument("--gene-expression-pcs", "-x", required = False, default = None,
                        help="Principal components of gene expression vectors (samples x PCs) [TSV].")
    parser.add_argument('--window', '-w', default=1000000, type=numpy.int32,
                        help='Cis-window size, in bases. Default: 1000000.')
    parser.add_argument('--seed', default=None, type=int,
                        help='Seed for permutations.')
    parser.add_argument('--fdr', default=float(0.05), type=float, help='FDR for cis-QTLs')
    parser.add_argument('--qvalue_lambda', default=None, type=float, help='lambda parameter for pi0est in qvalue.')
    parser.add_argument(
        "--cis-window-position", default="mid", dest="cis_win_pos",
        help="Position of CIS window in gene [mid|tss]. mid: mid-point of annotated gene. tss: transcription start site."
    )
    return parser.parse_args()

class TestArgs:
    def __init__(self):
        datadir = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/handover/franke_data4/aggrsce_out"
        self.gene_expression_bed = os.path.join(datadir, "./normcount_allchrAll.bed.gz")
        self.plink_gt_fnprfx = os.path.join(datadir, "./franke_gt_plink/franke_gt_plink_bed")
        self.genotype_pcs_plink = os.path.join(datadir, "./franke_gt_plink/franke_gtpca.tsv")
        self.gene_expression_pcs = os.path.join(datadir, "./normcount_all_PC20.tsv")
        self.oufn_prfx = 'test_out' # prefix for output files
        self.output_dir = os.curdir
        self.seed = SEED
        self.window = 1000000
        self.fdr = 0.05
        self.qvalue_lambda = None
        self.cis_win_pos = 'mid'

if __name__ == '__main__':
    print(sys.path)
    #help(tensorqtl)
    print('PyTorch {}'.format(torch.__version__))
    print('Pandas {}'.format(pandas.__version__))

    #args = TestArgs()
    args = set_argument_parser()

    # the following requires gzipped BED file
    logger = tensorqtl.SimpleLogger(os.path.join(args.output_dir, f'{args.oufn_prfx}.tensorQTL.cis.log'))
    logger.write(f'[{datetime.now().strftime("%b %d %H:%M:%S")}] Running TensorQTL: cis-QTL mapping')
    if torch.cuda.is_available():
        logger.write(f'  * using GPU ({torch.cuda.get_device_name(torch.cuda.current_device())})')
    else:
        logger.write('  * WARNING: using CPU!')
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if args.seed is not None:
        logger.write(f'  * using seed {args.seed}')

    # load inputs
    logger.write(f'  * reading phenotypes ({args.gene_expression_bed})')
    pheno_df, pheno_pos_df = tensorqtl.read_phenotype_bed(args.gene_expression_bed)
    # read_phenotype_bed() expects gzipped BED file with start = end - 1 defining the centre of the CIS window

    logger.write(f'  * reading covariates ({args.genotype_pcs_plink})')
    cov_gt_df = pandas.read_csv(args.genotype_pcs_plink, sep='\t', index_col=1).drop('#FID', axis=1)
    cov_gt_df.rename(lambda x: f'gt_{x}', axis = 1, inplace = True)

    if args.gene_expression_pcs:
        logger.write(f'  * adding covariates ({args.gene_expression_pcs})')
        cov_pheno_df = pandas.read_csv(args.gene_expression_pcs, sep='\t')
        cov_pheno_df.rename(lambda x: f'gx_{x}', axis = 1, inplace = True)
        cov_df = pandas.concat([cov_gt_df.loc[pheno_df.columns], cov_pheno_df], axis = 1)
    else:
        # cov_df = cov_gt_df[pheno_df.columns]
        cov_df = cov_gt_df.loc[pheno_df.columns]

    logger.write(f'  * reading genotypes from PLINK ({args.plink_gt_fnprfx})')
    plink_reader = tensorqtl.genotypeio.PlinkReader(args.plink_gt_fnprfx)
    genotype_df = plink_reader.load_genotypes()
    genotype_df = genotype_df[pheno_df.columns]
    variant_df = plink_reader.bim.set_index('snp')[['chrom', 'pos']]

    plink_chromosome_names = set(pandas.Categorical(variant_df["chrom"].values).categories)
    phenotype_chromosome_names = set(pandas.Categorical(pheno_pos_df["chr"].values).categories)
    common_chromosome_names = plink_chromosome_names & phenotype_chromosome_names
    if len(common_chromosome_names) < 1:
        logger.write(
            "ERROR: There are no common chromosome names between the PLINK genotypes"
            " and the gene annotation.\n"
            )
    print("common chromosome names:", common_chromosome_names)
    # remove X,Y-chromosome
    # all chromosomes must be present in variant_df !
    #sel = (pheno_pos_df['chr'] != 'chrY') & (pheno_pos_df['chr'] != 'chrX')
    sel = pheno_pos_df["chr"].isin(common_chromosome_names)
    pheno_df = pheno_df[sel]
    pheno_pos_df = pheno_pos_df[sel]

    # cis_df = tensorqtl.cis.map_nominal(
    #     genotype_df,
    #     variant_df,
    #     pheno_df,
    #     pheno_pos_df,
    #     covariates_df=covariates_df
    #     )

    logger.write(f'  * mapping cis eQTL (window = {args.window}')
    cis_df = tensorqtl.cis.map_cis(
        genotype_df,
        variant_df,
        pheno_df,
        pheno_pos_df,
        covariates_df=cov_df,
        window=args.window,
        nperm=1000,
        maf_threshold=0,
        logger = logger,
        seed=args.seed
        )
    out_file = os.path.join(args.output_dir, args.oufn_prfx+'.cis_qtl.tsv.gz')
    cis_df.to_csv(out_file, sep='\t', float_format='%.6g')

    calculate_qvalues(cis_df, fdr=args.fdr, qvalue_lambda=args.qvalue_lambda, logger=logger)

    out_file_qval = os.path.join(args.output_dir, args.oufn_prfx+'.cis_qtl_qval.tsv.gz')
    cis_df.to_csv(out_file_qval, sep='\t', float_format='%.6g')

    sys.exit(0)
