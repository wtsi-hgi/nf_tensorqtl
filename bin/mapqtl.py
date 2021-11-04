#!/usr/bin/env python3

## CIS eQLT mapping with tensorQTL
## see https://github.com/broadinstitute/tensorqtl

DEBUG = True

import sys
import os
import argparse
from datetime import datetime
import torch
import numpy
import scipy.stats as stats
import pandas
import tensorqtl

sys.path.insert(1, os.path.dirname(__file__))
from post import calculate_qvalues

SEED = 123456

FMT_FLOAT_TSV = "%.6g"
COLNAM_GENE_SYMBOLS = "hgnc_symbols"

def get_cis_window_centres_from_gene_annotation(gene_annot_df, cis_centre = 'mid', logger = None):
    if logger is None:
        logger = tensorqtl.SimpleLogger()
    cis_pos_logstr = ''
    has_tss = "tss" in gene_annot_df.columns
    if has_tss:
        tss_logstr = "have"
        cis_pos_logstr = 'TSS'
        colnames = ["chromosome", "tss"]
        colmapper = {"chromosome":"chr"}
    else:
        tss_logstr = "don't have"
        colnames = ["chromosome", "start"]
        colmapper = {"chromosome":"chr", "start": "tss"}
    pos_df = gene_annot_df[colnames].rename(columns = colmapper)
    if cis_centre == 'mid':
        cis_pos_logstr = "gene midpoint"
        # set start = end - 1 to gene mid-point
        pos_df["tss"] = ((pos_df["tss"] + gene_annot_df["end"])/2).round().astype('int64')
    elif cis_centre != 'tss' or not has_tss: # ref == 'start'
        cis_pos_logstr = "gene start"
        s = gene_annot_df["strand"] < 0 # start/end annotation is on forward strand
        pos_df["tss"][s] = gene_annot_df["end"][s]
    logger.write(f"  * {tss_logstr} TSS annotation. CIS window centred at '{cis_pos_logstr}'.")
    return pos_df

def load_genexpression_data(fn_gx, fn_ga, cis_centre = 'mid', logger = None):
    if logger is None:
        logger = tensorqtl.SimpleLogger()
    logger.write(f'  * reading phenotypes ({fn_gx})')
    df_gx = pandas.read_table(fn_gx)
    logger.write(f'  * reading gene boundaries ({fn_ga})')
    df_ga = pandas.read_table(fn_ga, index_col = 1)
    index_common = df_gx.index.intersection(df_ga.index)
    logger.write('  * {:d} of {:d} genes are annotated.'.format(len(index_common), df_gx.shape[0]))
    df_gx = df_gx.loc[index_common]
    df_ga = df_ga.loc[index_common]
    df_pos = get_cis_window_centres_from_gene_annotation(
        df_ga,
        cis_centre = cis_centre,
        logger = logger)
    df_names = None
    if COLNAM_GENE_SYMBOLS in df_ga:
        df_names = df_ga[COLNAM_GENE_SYMBOLS]
    return df_gx, df_pos, df_names

def load_PLINK_genotypes(plink_gt_fnprfx, logger):
    logger.write(f'  * reading genotypes from PLINK ({plink_gt_fnprfx})')
    plink_reader = tensorqtl.genotypeio.PlinkReader(plink_gt_fnprfx)
    genotype_df = plink_reader.load_genotypes()
    genotype_df = genotype_df[pheno_df.columns]
    variant_df = plink_reader.bim.set_index('snp')[['chrom', 'pos']]
    return variant_df, genotype_df

def check_chromosome_names_against_PLINK(plink_variant_df, pheno_df, cispos_df, logger):
    plink_chromosome_names = set(pandas.Categorical(plink_variant_df["chrom"].values).categories)
    phenotype_chromosome_names = set(pandas.Categorical(cispos_df["chr"].values).categories)
    common_chromosome_names = plink_chromosome_names & phenotype_chromosome_names
    if len(common_chromosome_names) < 1:
        sys.exit(
            "ERROR: There are no common chromosome names between the PLINK genotypes"
            " and the gene annotation."
            )
    logstr = "  * common chromosome names:\n"
    i = 0
    for chrnam in common_chromosome_names:
        i += 1
        logstr += f"    [{i}] {chrnam}\n"
    logger.write(logstr)

    # remove X,Y-chromosome
    # all chromosomes must be present in variant_df !
    #sel = (pheno_pos_df['chr'] != 'chrY') & (pheno_pos_df['chr'] != 'chrX')
    sel = cispos_df["chr"].isin(common_chromosome_names)
    return pheno_df[sel], cispos_df[sel]

def load_covariates(plink_gt_pcs, fn_gx_pcs = None, logger = None):
    if logger is None:
        logger = tensorqtl.SimpleLogger()
    logger.write(f'  * reading covariates ({plink_gt_pcs})')
    cov_gt_df = pandas.read_csv(plink_gt_pcs, sep='\t', index_col=1).drop('#FID', axis=1)
    cov_gt_df.rename(lambda x: f'gt_{x}', axis = 1, inplace = True)

    if fn_gx_pcs:
        logger.write(f'  * adding covariates ({fn_gx_pcs})')
        cov_pheno_df = pandas.read_csv(fn_gx_pcs, sep='\t')
        cov_pheno_df.rename(lambda x: f'gx_{x}', axis = 1, inplace = True)
        cov_df = pandas.concat([cov_gt_df.loc[pheno_df.columns], cov_pheno_df], axis = 1)
    else:
        # cov_df = cov_gt_df[pheno_df.columns]
        cov_df = cov_gt_df.loc[pheno_df.columns]
    return cov_df

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
        "--cis-window-position", default="mid", dest="cis_window_pos",
        help="Position of CIS window in gene [mid|start]. mid: mid-point of annotated gene. start: gene start as annotated (~transcription start site)."
    )
    parser.add_argument('--gene-annotation', required = True,
                        help="Gene annotation file with Ensembl Gene IDs and HGNC symbols",
                        dest = "gene_annot_tsv")

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
        self.cis_window_pos = 'mid'
        self.gene_annot_tsv= "/lustre/scratch123/hgi/mdt1/projects/ukbb_scrna/pipelines/Pilot_UKB/handover/franke_data4/gene_annot_chr.tsv"

if __name__ == '__main__':
    if DEBUG:
        print(sys.path)
        #help(tensorqtl)
        print('PyTorch {}'.format(torch.__version__))
        print('Pandas {}'.format(pandas.__version__))
        print("directory: ", os.path.abspath(os.curdir))
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
    #pheno_df, pheno_pos_df = tensorqtl.read_phenotype_bed(args.gene_expression_bed)
    # read_phenotype_bed() expects gzipped BED file with start = end - 1 defining the centre of the CIS window
    pheno_df, cispos_df, gene_symbols = load_genexpression_data(
            args.gene_expression_bed, args.gene_annot_tsv,
            cis_centre = args.cis_window_pos, logger = logger)

    cov_df = load_covariates(args.genotype_pcs_plink, fn_gx_pcs = args.gene_expression_pcs, logger = logger)

    variant_df, genotype_df = load_PLINK_genotypes(args.plink_gt_fnprfx, logger)
    pheno_df, cispos_df = check_chromosome_names_against_PLINK(variant_df, pheno_df, cispos_df, logger)

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
        cispos_df,
        covariates_df=cov_df,
        window=args.window,
        nperm=1000,
        maf_threshold=0,
        logger = logger,
        seed=args.seed
        )
    out_file_bin = os.path.join(args.output_dir, args.oufn_prfx+'.bin.gz')
    cis_df.to_pickle(out_file_bin)
    out_file = os.path.join(args.output_dir, args.oufn_prfx+'.cis_eqtl.tsv.gz')
    cis_df.to_csv(out_file, sep='\t', float_format=FMT_FLOAT_TSV)

    cis_df.insert(0,"chr", cispos_df.loc[cis_df.index]["chr"])
    cis_df.insert(1,"cis_window_pos", cispos_df.loc[cis_df.index]["tss"])

    if gene_symbols is not None:
        cis_df.insert(0,"gene_symbol", gene_symbols.loc[cis_df.index])


    # remove NAns in predicted values
    # otherwise the call scipy.stats.pearsonr(cis_df['pval_perm'], cis_df['pval_beta'])
    # in tensorqtl post.calculate_qvalues() throws an error
    out_file_dropped = os.path.join(args.output_dir, args.oufn_prfx+'.cis_eqtl_dropped.tsv.gz')
    sv = numpy.isnan(cis_df['pval_beta'])
    logger.write("  * dropping {:d} variants withouth Beta-approximated p-values to {:s}\n.".format(
        sum(sv), out_file_dropped))
    cis_df_dropped = cis_df.loc[sv]
    cis_df_dropped.to_csv(out_file_dropped, sep='\t', float_format=FMT_FLOAT_TSV)

    cis_df.drop(index = cis_df_dropped.index, inplace = True)

    logger.write(f'  * Number of phenotypes tested: {cis_df.shape[0]}')
    r = stats.pearsonr(cis_df['pval_perm'], cis_df['pval_beta'])[0]
    logger.write(f'  * Correlation between Beta-approximated and empirical p-values: : {r:.4f}')

    tensorqtl.calculate_qvalues(cis_df, fdr=args.fdr, qvalue_lambda=args.qvalue_lambda, logger=logger)
    out_file_qval= os.path.join(args.output_dir, args.oufn_prfx+'.cis_eqtl_qval.tsv.gz')
    cis_df.to_csv(out_file_qval, sep='\t', float_format=FMT_FLOAT_TSV)

    sys.exit(0)
