#!/usr/bin/env python3

DEBUG = True
## aggregate scRNA-seq expression levels per donor over all cells of one cell-type

import sys
import os
import argparse
import gzip
import numpy
import anndata
import pandas

# input file format for gene annotations:
# EnsEMBL gene region annotation downloaded as TSV file using BioMart
# Format:
# Gene stable ID  Gene stable ID version  Gene name       Chromosome/scaffold name        Gene start (bp) Gene end (bp)   Strand
# ENSG00000284662 ENSG00000284662.1       OR4F16  1       685679  686673  -1
# ENSG00000186827 ENSG00000186827.11      TNFRSF4 1       1211340 1214153 -1
# ...

STRTRANSL_SPACE_TO_UNDERSCORE = str.maketrans({' ':'_'})
FMT_BEDFNAM_COUNT_MARIX = "{:s}_{:s}_counts.bed.gz"
FMT_FNAM_CELLNUM_TABLE = "{:s}_{:s}_ncells.tsv"
FNAM_GENE_IDS = "gene_ids_without_annotation.txt"
FNAM_FILE_TABLE = "celltype_files.tsv"


COUNT_MATRIX_COLUMN_DICT = {
    "chromosome": "#chr",
    "query_id": "gene_id"
    }

class DataHandover:
    def __init__(self):
        self.data_dir = ''
        self.ftab = None # pandas data frame [experiment_id,donor_id,filename_h5ad,filename_annotation_tsv]
        self.anndata_objects = ()
        self.cell_types = set()
        self.gene_names = None
        self.gene_annot_df = None
        self.failed_gene_idx = set()

    def load_file_table(self, data_dir, donor_only = True):
        self.data_dir = data_dir
        fnam = os.path.join(data_dir, "files.tsv")
        self.ftab = pandas.read_table(fnam).set_index(["experiment_id", "donor_id"], drop = False)
        if donor_only:
            self.ftab.drop(index = ["doublet", "unassigned"], level=1, inplace = True)

        sys.stderr.write("# list of {:d} files loaded.\n".
                            format(len(self.ftab.shape))
                            )

        return self.ftab

    def get_cell_types(self):
        celltypes = set()
        for ad in self.anndata_objects:
            cts = set(pandas.Categorical(ad.obs['azimuth.celltyp.l2']).categories)
            celltypes = celltypes.union(cts)
        return celltypes

    def check_gene_ids(self):
        sys.stderr.write("# checking gene IDs across input files ...\n")
        gene_names_list = None
        for ad in self.anndata_objects:
            gene_names = ad.var.index
            if gene_names_list is None:
                gene_names_list = gene_names
            elif not (gene_names == gene_names_list).all():
                sys.exit(
                "ERROR: inconsistent gene name vector in file {}".
                format(row['filename_h5ad'])
                )
        return gene_names_list

    def load_count_files(self):
        anndata_objects = []
        for i in self.ftab.index:
            fpath = os.path.join(self.data_dir, self.ftab['filename_h5ad'][i])
            sys.stderr.write("# loading file {} ...\n".format(fpath))
            anndata_objects.append(anndata.read_h5ad(fpath))
        self.anndata_objects = tuple(anndata_objects)
        self.gene_names = self.check_gene_ids()
        self.cell_types = self.get_cell_types()
        return self.anndata_objects

    def check_gene_annot(self, gene_annotd):
        failed_gene_ids = []
        ctr = 0
        for gn in self.gene_names:
            if gn not in gene_annotd.index:
                ctr += 1
                sys.stderr.write("WARNING: [{:d}] gene {} not found in annotation.\n".format(ctr, gn))
                failed_gene_ids.append(gn)
        sys.stderr.write("# {:d} out of {:d} genes have no annotation.\n".format(ctr, len(self.gene_names)))
        return set(failed_gene_ids)

    def add_gene_annotation(self, fnam):
        sys.stderr.write(
            "# loading EnsEMBL/BioMart gene annotations from file {}\n".
            format(args.gene_annot_file))

        self.gene_annot_df = pandas.read_table(fnam, index_col = 1)

        sys.stderr.write(
            "# annotations loaded for {:d} genes.\n".
            format(self.gene_annot_df.shape[0])
            )

        self.failed_gene_ids = self.check_gene_annot(self.gene_annot_df)
        return self.failed_gene_ids

    def _f_sum_over_celltype(self, ad, celltyp = 'all'):
        sv_qcpass = ad.obs['qc.filter.pass'] == 'True'
        if celltyp is None or celltyp == 'all':
            # generate pseudo-bulk count
            n_cells_celltype = len(sv_qcpass)
        else:
            # cell-type specific count
            sv_celltype = ad.obs['azimuth.celltyp.l2'] == celltyp
            n_cells_celltype = sum(sv_celltype)
            sv_qcpass = sv_qcpass & sv_celltype
        n_cells_qc = sum(sv_qcpass)

        mr = numpy.sum(ad[sv_qcpass,].X, axis = 0, dtype = 'int')
        df = pandas.DataFrame(mr, columns = ad.var.index)
        # type(mr): <class 'numpy.matrix'>
        return df, n_cells_celltype, n_cells_qc

    def write_unannotated_gene_ids(self, oufn = FNAM_GENE_IDS):
        with open(oufn, 'w') as oufh:
            for gn in self.failed_gene_ids:
                oufh.write("{:s}\n".format(gn))
        return len(self.failed_gene_ids)

    def sum_counts_over_celltype(self, celltype = 'all'):
        count_df_lst = []
        n_cells_celltype_lst = []
        n_cells_qc_pass_lst = []
        for i in range(len(self.ftab)):
            row = self.ftab.iloc[i]
            sum_df, n_cells_celltype, n_cells_qc = self._f_sum_over_celltype(self.anndata_objects[i], celltype)
            count_df_lst.append(sum_df)
            n_cells_celltype_lst.append(n_cells_celltype)
            n_cells_qc_pass_lst.append(n_cells_qc)
        counts_df = pandas.concat(count_df_lst, axis = 0, join = 'inner', ignore_index = True).set_index(self.ftab["donor_id"])

        donor_df = pandas.DataFrame({
            "n_cells": n_cells_celltype_lst,
            "n_cells_qc_pass": n_cells_qc_pass_lst
            },
            index = self.ftab.index
        )

        return counts_df, donor_df

    def aggregate_over_celltype(self, celltype = 'all', n_cells_min = int(0), outdir = os.curdir, oufprfx = 'aggrsc'):
        n_lines = 0
        count_matrix, donor_df = self.sum_counts_over_celltype(celltype = celltype)
        celltypestr = celltype.translate(STRTRANSL_SPACE_TO_UNDERSCORE)
        oufn_matrix = FMT_BEDFNAM_COUNT_MARIX.format(oufprfx, celltypestr)
        oufpath_matrix = os.path.join(outdir, oufn_matrix)
        oufn_donors = FMT_FNAM_CELLNUM_TABLE.format(oufprfx, celltypestr)
        oufpath_donors = os.path.join(outdir, oufn_donors)
        if count_matrix is None:
            sys.stderr.write("WARNING: no counts for celltype {:s}\n".format(celltypestr))
        else:
            anno_mt = pandas.concat(
                [dh.gene_annot_df[["chromosome", "start", "end", "query_id"]].rename(COUNT_MATRIX_COLUMN_DICT, axis = 1),
                count_matrix.T],
                axis = 1,
                join = "inner" # only the genes that are annotated
                )
            n_lines = len(anno_mt)
            anno_mt.to_csv(oufpath_matrix, sep = "\t")
            sys.stderr.write("# {:d} lines written to file {:s} ...\n".format(n_lines, oufpath_matrix))

        donor_df.to_csv(oufpath_donors, sep = '\t')
        sys.stderr.write("# {:d} lines written to file {:s} ...\n".format(len(donor_df), oufpath_donors))

        n_donors = sum(donor_df["n_cells_qc_pass"] >= 0)
        n_donors_with_min_n_cells = sum(donor_df["n_cells_qc_pass"] >= n_cells_min)

        return n_lines, n_donors, n_donors_with_min_n_cells, oufn_matrix, oufn_donors,

def set_argument_parser():
    parser = argparse.ArgumentParser(description="Aggregate gene/transcript counts over cells of one class/type.")

    parser.add_argument("--output-dir", "-o",
        default=os.curdir,
        help="output directory",
        dest="outdir")

    parser.add_argument("--gene-annotations", required = True,
                        help="EnsEMBL/BioMart file of gene annotations [TSV]",
                        dest="gene_annot_file")
    parser.add_argument("--handover-data-dir", "-d", default = os.curdir,
                        help="Output directory of the Azimuth celltype-assignment pipeline.",
                        dest="datadir_handover")
    parser.add_argument("--output-file-list", default = None,
                        help="Write list of (filtered) celltypes and corresponding output files [CSV].",
                        dest="oufn_filtered_file_lst")
    parser.add_argument("--output-file-prefix", default = 'aggrsce',
                        help="Prefix for output files.",
                        dest="oufn_prfx")
    parser.add_argument(
        "--minium-cell-number", "-c", type = int, default = int(0), dest="n_cells_min",
        help="Minimum number of cells for donor to be included in analysis."
        )
    parser.add_argument(
        "--minimum-donor-number", "-n", type = int, default = int(0), dest="n_donors_min",
        help="Minimum number of donors with a number of cells at or above the threshold set by --minium-cell-number."
        )
    return parser.parse_args()

class TestArgs:
    def __init__(self):
        self.gene_annot_file = 'gene_annot.tsv'
        self.datadir_handover = 'franke_data4'
        self.oufn_gene_ids = 'gene_ids_unannotated.txt'

if __name__ == '__main__':
    # args = TestArgs()
    args = set_argument_parser()
    if args.outdir != os.curdir and not os.access(args.outdir, os.F_OK):
        os.mkdir(args.outdir)

    if DEBUG:
        sys.stderr.write("current directory: {}\n".format(os.path.abspath(os.curdir)))
        for fn in os.listdir():
            sys.stderr.write("- ls: {}\n".format(fn))
        for fp in os.listdir(args.datadir_handover):
            sys.stderr.write("-- ls: {}\n".format(fp))

    oufn_unannotated_genes = os.path.join(args.outdir,'{:s}_{:s}'.format(args.oufn_prfx, FNAM_GENE_IDS))
    oufn_file_list = os.path.join(args.outdir,'{:s}_{:s}'.format(args.oufn_prfx, FNAM_FILE_TABLE))

    dh = DataHandover()
    dh.load_file_table(args.datadir_handover, donor_only = True)
    dh.load_count_files()
    dh.add_gene_annotation(args.gene_annot_file)
    dh.write_unannotated_gene_ids(oufn_unannotated_genes)

    # print(celltypes)
    ctr = 0
    celltypes = ['all'] + list(dh.cell_types)
    bed_file_lst = []
    tsv_file_lst = []
    n_donors_lst = []
    n_donors_min_cells_lst = []

    for celltyp in celltypes:
        ctr += 1
        sys.stderr.write("processing celltype '{:s}' ...\n".format(celltyp))
        n_lines, n_donors, n_donors_with_min_n_cells, oufn_matrix, oufn_donor_table = \
            dh.aggregate_over_celltype(
                celltyp,
                n_cells_min = args.n_cells_min,
                outdir = args.outdir,
                oufprfx = args.oufn_prfx
                )
        sys.stderr.write(f"# [{ctr}]\t{celltyp}\t{n_lines}\t{n_donors}"
            f"\t{n_donors_with_min_n_cells}\t{oufn_matrix}\t{oufn_donor_table}\n")
        bed_file_lst.append(oufn_matrix)
        tsv_file_lst.append(oufn_donor_table)
        n_donors_lst.append(n_donors)
        n_donors_min_cells_lst.append(n_donors_with_min_n_cells)
    celltype_df = pandas.DataFrame({
        'celltype': celltypes,
        'n_donors': n_donors_lst,
        'n_donors_min_cells': n_donors_min_cells_lst,
        'count_file': bed_file_lst,
        'donor_file': tsv_file_lst})

    celltype_df.to_csv(oufn_file_list, sep = '\t')
    sys.stderr.write(f"# List of output file written to file {oufn_file_list}.\n")

    if args.oufn_filtered_file_lst:
        sel = celltype_df["n_donors_min_cells"] > args.n_donors_min
        celltype_df[sel].to_csv(args.oufn_filtered_file_lst, columns = ["celltype", "count_file", "donor_file"], index = False)
        sys.stderr.write("# Filtered list of {:d} celltype files written to file {:s}.\n"
            .format(sum(sel), args.oufn_filtered_file_lst))
        sys.stderr.write("# This list contains files for celltypes with at least {:d} donors"
            " each having a minimum of {:d} QC'ed cells".format(args.n_donors_min, args.n_cells_min))

    exit(0)
