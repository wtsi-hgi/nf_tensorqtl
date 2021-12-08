#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def read_reftab(fnam):
    df = pd.read_csv(fnam, header = [0,1])
    # fix multi-indexed columns
    last_nam = ''
    tups = []
    for c in df.columns:
      if c[0][:8] != 'Unnamed:' or last_nam is None:
        last_nam = c[0]
      tups.append((last_nam, c[1]))
    df.columns = pd.MultiIndex.from_tuples(tups, names=["param", "type"])

    rds = df[[("eQTL", "Gene"), ("eQTL", "SNP"), ("P-values", "PBMC"), ("Significant at FDR < 0.05", "PBMC")]]
    rds.columns = pd.Index(["ensembl_id", "rsid_r", "pval_r", "signif05_r"])
    # there are multiple variants per gene, keep only one as in the tensorQTL output
    rds = rds.sort_values("pval_r", axis = 0) \
        .groupby("ensembl_id").first() \
        .sort_values("pval_r", axis = 0)
    rds.insert(2, "rank_r", pd.Series(range(1, rds.shape[0]+1), index = rds.index))
    return rds

def read_qtltab(fnam):
    qdf = pd.read_table(fnam)
    qds = qdf.rename(
        {"phenotype_id": "ensembl_id",
            "variant_id": "rsid_q",
            "pval_true_df": "pval_q",
            "qval": "qval_q"
        },
        axis = 1
    )
    # qds = qdf[["phenotype_id","variant_id", "pval_true_df", "qval"]]
    qds = qds[["ensembl_id", "rsid_q", "pval_q", "qval_q"]].sort_values("pval_q", axis = 0)
    qds.insert(2, 'rank_q', pd.Series(range(1, qds.shape[0]+1), index = qds.index))
    return qds

def join_and_filter(qdf, rdf):
    mdf = pd.merge(qdf, rdf, on = 'ensembl_id', how = "outer")

    # select only those rows that represent significant entries in either data set
    mf = mdf[(mdf["qval_q"] <= 0.05) | (mdf['signif05_r'] == '*')]
    # set NaNs to rank 0
    zr = np.isnan(mf["pval_r"])
    mf.loc[zr, "rank_r"] = 0
    zq = np.isnan(mf["pval_q"])
    mf.loc[zq, "rank_q"] = 0

    return mf.convert_dtypes(convert_floating = False).set_index("ensembl_id")

def make_rank_plots(mf, fnam):
    # make a scatter plot of ranks
    fig, axes = plt.subplots(nrows=2, ncols=2, constrained_layout=True)

    x = mf["rank_r"].values
    y = mf["rank_q"].values
    xmax = x.max()

    axes[0,0].set_ylabel('tensorQTL')
    axes[0,0].set_title("rank scatter of gene associations")
    axes[0,0].plot(x, y, 'bo')

    axes[1,0].set_xlabel('franke')
    axes[1,0].set_ylabel('tensorQTL')
    axes[1,0].plot(x, y, 'bo')
    axes[1,0].axis([0, xmax, 0, xmax])
    # make 2 line plots of ranks

    # m = mf[["rank_r", "rank_q"]].values
    axes[0,1].plot(x, label='franke')
    axes[0,1].plot(y, label='tensorQTL')
    axes[0,1].set_title("ranks of gene associations")
    axes[0,1].legend()


    axes[1,1].plot(x, label='franke')
    axes[1,1].plot(y, label='tensorQTL')
    axes[1,1].set_xlabel('franke')
    axes[1,1].axis([0, xmax, 0, xmax])

    fig.savefig("{:s}.pdf".format(fnam))
    return fig

if __name__ == '__main__':

    nargs = len(sys.argv)

    if nargs != 4:
        sys.exit("usage: {:s} <ref QTL table [CSV]> <new QTL table [TSV]> <plot filename [w/o extension]>".format(argv[0]))

    fnreftab = sys.argv[1]
    fnqtltab = sys.argv[2]
    oufn = sys.argv[3]

    rdf = read_reftab(fnreftab)
    qdf = read_qtltab(fnqtltab)
    mf = join_and_filter(qdf, rdf)
    fig = make_rank_plots(mf, oufn)

    exit(0)
