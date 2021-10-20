#!/usr/bin/env python

DEBUG = True
## aggregate scRNA-seq expression levels per donor over all cells of one cell-type

import sys
import os
import csv
import argparse
import gzip
import numpy
import anndata
import pandas

## input files
DATA_DIR_HANDOVER = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/handover/franke_data4"
# franke data by default

GENE_ANNOT_FILE = os.path.join(DATA_DIR_HANDOVER, "gene_annot.tsv")
# EnsEMBL gene region annotation downloaded as TSV file using BioMart
# Format:
# Gene stable ID  Gene stable ID version  Gene name       Chromosome/scaffold name        Gene start (bp) Gene end (bp)   Strand
# ENSG00000284662 ENSG00000284662.1       OR4F16  1       685679  686673  -1
# ENSG00000186827 ENSG00000186827.11      TNFRSF4 1       1211340 1214153 -1
# ...

CHRNAM_PREFIX_BED = 'chr' # prefix for chromosome names in BED output file
STRTRANSL_SPACE_TO_UNDERSCORE = str.maketrans({' ':'_'})
FMT_BEDFNAM_COUNT_MARIX = "count_{:s}.bed.gz"
FMT_FNAM_CELLNUM_TABLE = "ncells_{:s}.tsv"
FNAM_GENE_IDS = "gene_ids_without_annotation.txt"
HEADER_BED_FILE = "chr\tstart\tend\tgene_id"
HEADER_CELLNUM_TABLE = "experiment_id\tdonor_id\tsample_name\tcell_type\tn_cells\tn_cells_pass_qc\tdonor_included\tfile_name"
FMT_CELLNUM_TABLE_ROW = "{0[0]:s}\t{0[1]:s}\t{0[2]:s}\t{0[3]:s}\t {0[4]:d}\t{0[5]:d}\t{0[6]}\t{0[7]:s}"

def set_argument_parser():
    parser = argparse.ArgumentParser(description="Aggregate gene/transcript counts over cells of one class/type.")

    parser.add_argument("--output-dir", "-o",
        default=os.curdir,
        help="output directory",
        dest="outdir")

    parser.add_argument("--gene-annotations", default = GENE_ANNOT_FILE,
                        help="EnsEMBL/BioMart file of gene annotations [TSV]",
                        dest="gene_annot_file")
    parser.add_argument("--handover-data-dir", "-d", default = DATA_DIR_HANDOVER,
                        help="Output directory of the Azimuth celltype-assignment pipeline.",
                        dest="datadir_handover")
    parser.add_argument("--gene-ids-output-file", default = None,
                        help="Output file for gene ids.",
                        dest="oufn_gene_ids")
    return parser.parse_args()



def load_gene_annotation_table(fnam):
    gad = {}
    if len(fnam) > 3 and fnam[-3:] == ".gz":
        import gzip
        infh = gzip.open(fnam, mode='rt')
    else:
        infh = open(fnam, 'r')
    for row in csv.DictReader(infh, delimiter='\t'):
        chrnam = row['chromosome']
        start = int(row['start'])
        end = int(row['end'])
        id = row['query_id']
        strand = int(row['strand'])
        dt = (chrnam, start, end, strand, row['ensembl_id'])
        try:
            gad[id].append(dt)
        except KeyError:
            gad[id] = [dt]
    infh.close()
    return gad


class DataHandover:
    def __init__(self):
        self.data_dir = ''
        self.ftab = () # contains tuple of row dictionaries, row <-> donor
        self.anndata_objects = ()
        self.celltypes = set()
        self.gene_names = ()
        self.failed_gene_idx = set()

    def load_file_table(self, data_dir, donor_only = True):
        self.data_dir = data_dir
        ftab = []
        fnam = os.path.join(self.data_dir, "files.tsv")

        sys.stderr.write(
            "# loading list of data files from {}\n".
            format(fnam))

        infh = open(fnam, 'r')
        for row in csv.DictReader(infh, delimiter='\t'):
            donor_id = row['donor_id']
            if not donor_only or (donor_id != 'unassigned' and donor_id != 'doublet'):
                ftab.append(row)
        infh.close()

        sys.stderr.write("# list of {:d} files loaded.\n".
            format(len(ftab)))
        self.ftab = tuple(ftab)

        return self.ftab

    def load_count_files(self):
        anndata_objects = []
        for row in self.ftab:
            fpath = os.path.join(self.data_dir, row['filename_h5ad'])
            sys.stderr.write("# loading file {} ...\n".format(fpath))
            anndata_objects.append(anndata.read_h5ad(fpath))
        self.anndata_objects = tuple(anndata_objects)
        return self.anndata_objects

    def get_cell_types(self):
        celltypes = set()
        for ad in self.anndata_objects:
            cts = set(pandas.Categorical(ad.obs['azimuth.celltyp.l2']).categories)
            celltypes = celltypes.union(cts)
        self.celltypes = celltypes
        return self.celltypes

    def check_gene_ids(self):
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
        self.gene_names = tuple(gene_names_list)
        return self.gene_names

    def check_gene_annot(self, gene_annotd):
        failed_gene_ids = []
        ctr = 0
        for gn in self.gene_names:
            if gn not in gene_annotd:
                ctr += 1
                sys.stderr.write("WARNING: [{:d}] gene {} not found in annotation.\n".format(ctr, gn))
                failed_gene_ids.append(gn)
        sys.stderr.write("# {:d} out of {:d} genes have no annotation.\n".format(ctr, len(self.gene_names)))
        self.failed_gene_ids = set(failed_gene_ids)
        return len(self.failed_gene_ids)

    def _f_sum_over_celltype(self, ad, celltyp = 'all'):
        selv_qcpass = ad.obs['qc.filter.pass'] == 'True'
        if celltyp is None or celltyp == 'all':
            # generate pseudo-bulk count
            n_cells_celltype = sum(selv_qcpass)
        else:
            # cell-type specific count
            selv_celltype = ad.obs['azimuth.celltyp.l2'] == celltyp
            n_cells_celltype = len(selv_celltype)
            selv_qcpass = selv_qcpass & selv_celltype
        n_cells_qc = sum(selv_qcpass)

        mr = numpy.sum(ad[selv_qcpass,].X, axis = 0, dtype = 'int')
        return mr, n_cells_celltype, n_cells_qc

    def write_failed_gene_ids(self, oufn = FNAM_GENE_IDS):
        with open(oufh, 'w') as oufh:
            for gn in self.failed_gene_ids:
                oufh.write("{:s}\n".format(gn))
        return len(self.failed_gene_ids)

    def write_cell_num_table(self, oufn, n_cells_per_donor_and_celltype):
        oufh = open(oufn, 'w')
        oufh.write(HEADER_CELLNUM_TABLE + '\n')
        linctr = 0
        for tup in n_cells_per_donor_and_celltype:
            oufh.write(FMT_CELLNUM_TABLE_ROW.format(tup) + '\n')
            linctr += 1
        oufh.close()
        return linctr

    def write_count_matrix(self, oufn, count_mtx, sample_names, gene_annot_dict = None):
        oufh = gzip.open(oufn, mode='wt')
        headerstr = HEADER_BED_FILE
        for snam in sample_names:
            headerstr += "\t{:s}".format(snam)
        oufh.write(headerstr + '\n')
        n_genes = len(self.gene_names)
        n_samples = len(sample_names)
        linctr = 0
        if count_mtx is None:
            sys.stderr.write("WARNING: empty count matrix in file {:s}.\n".format(oufn))
        else:
            for g in range(n_genes):
                gene_id = self.gene_names[g]
                if gene_annot_dict is None:
                    linstr = "{:s}".format(gene_id)
                else:
                    try:
                        ga = gene_annot_dict[gene_id]
                    except KeyError:
                        continue
                    else:
                        mcol = count_mtx[:,g]
                        if numpy.sum(mcol) < 1:
                            continue # skip gene if there are no counts
                        counts = list(mcol.flat)
                        linctr += 1
                        chrnam, start, end = ga[0][:3]
                        midpos = int((start + end)/2) # 0-based for BED file
                        linstr = "{:s}{}\t{:d}\t{:d}\t{:s}".format(CHRNAM_PREFIX_BED, chrnam, midpos, midpos+1, gene_id)
                for s in range(n_samples):
                    linstr += "\t{:d}".format(counts[s])
                oufh.write(linstr + '\n')
                linctr += 1
        oufh.close()
        return linctr

    def sum_counts_over_celltype(self, celltype = 'all',  n_cells_threshold = 5,
                                        gene_annotd = None, outdir = os.curdir):
        n_cells_per_donor_and_celltype = []
        celltypestr = celltype.translate(STRTRANSL_SPACE_TO_UNDERSCORE)
        oufn = os.path.join(outdir, FMT_BEDFNAM_COUNT_MARIX.format(celltypestr))
        count_matrix = None
        n_donors = len(self.ftab)
        sample_names = []
        for i in range(n_donors):
            r = self.ftab[i]
            experiment_id = r['experiment_id']
            donor_id = r['donor_id']
            # sample_name = "{:s}.{:s}".format(experiment_id, donor_id)
            sample_name = donor_id
            mr, n_cells_celltype, n_cells_qc = self._f_sum_over_celltype(self.anndata_objects[i], celltype)
            is_donor_included = n_cells_qc >= n_cells_threshold
            n_cells_per_donor_and_celltype.append(
                (experiment_id, donor_id, sample_name,
                    celltype, n_cells_celltype, n_cells_qc, is_donor_included, oufn)
                )
            if is_donor_included:
                sample_names.append(sample_name)
                if count_matrix is None:
                    count_matrix = mr
                else:
                    count_matrix = numpy.concatenate((count_matrix, mr), axis = 0)
        if count_matrix is None:
            sys.stderr.write("WARNING: no counts for celltype {:s}\n".format(celltypestr))
        linctr = self.write_count_matrix(oufn, count_matrix, sample_names, gene_annotd)
        sys.stderr.write("# {:d} lines writen to file {:s} ...\n".format(linctr, oufn))

        n_cells_per_donor_and_celltype = tuple(n_cells_per_donor_and_celltype)
        oufn = os.path.join(outdir, FMT_FNAM_CELLNUM_TABLE.format(celltypestr))
        linctr = self.write_cell_num_table(oufn, n_cells_per_donor_and_celltype)
        sys.stderr.write("# {:d} lines writen to file {:s} ...\n".format(linctr, oufn))

        return n_cells_per_donor_and_celltype

def set_argument_parser():
    parser = argparse.ArgumentParser(description="Aggregate gene/transcript counts over cells of one class/type.")

    parser.add_argument("--output-dir", "-o",
        default=os.curdir,
        help="output directory",
        dest="outdir")

    parser.add_argument("--gene-annotations", default = GENE_ANNOT_FILE,
                        help="EnsEMBL/BioMart file of gene annotations [TSV]",
                        dest="gene_annot_file")
    parser.add_argument("--handover-data-dir", "-d", default = DATA_DIR_HANDOVER,
                        help="Output directory of the Azimuth celltype-assignment pipeline.",
                        dest="datadir_handover")
    parser.add_argument("--gene-ids-output-file", default = None,
                        help="Output file for gene ids for which there was no annotation.",
                        dest="oufn_gene_ids")
    parser.add_argument("--min-cell-num", type = int, default = int(5),
                        help="Minimum number of cells per donor.",
                        dest="n_cells_min")
    return parser.parse_args()


if __name__ == '__main__':
    args = set_argument_parser()
    if args.outdir != os.curdir and not os.access(args.outdir, os.F_OK):
        os.mkdir(args.outdir)

    sys.stderr.write(
        "# loading EnsEMBL/BioMart gene annotations from file {}\n".
        format(args.gene_annot_file))
    gannotd = load_gene_annotation_table(args.gene_annot_file)
    sys.stderr.write(
        "# annotations loaded for {:d} genes.\n".
        format(len(gannotd)))

    dh = DataHandover()
    dh.load_file_table(args.datadir_handover, donor_only = True)
    dh.load_count_files()
    dh.check_gene_ids()
    dh.check_gene_annot(gene_annotd = gannotd)
    if (args.oufn_gene_ids):
        dh.write_failed_gene_ids(path.join(args.outdir, args.oufn_gene_ids))
    celltypes = dh.get_cell_types()
    # print(celltypes)
    dh.sum_counts_over_celltype(
        "all",
        n_cells_threshold = args.n_cells_min,
        gene_annotd = gannotd,
        outdir = args.outdir
        )

    for ct in celltypes:
        sys.stderr.write("processing celltype '{:s}' ...\n".format(ct))
        dh.sum_counts_over_celltype(
            ct,
            n_cells_threshold = args.n_cells_min,
            gene_annotd = gannotd,
            outdir = args.outdir
            )

    exit(0)
