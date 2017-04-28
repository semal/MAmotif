from optparse import OptionParser
from RNA_seq_analysis import *


def output_rpkm(cell_folder, refgenes_file):
    """
    using multiple process to get rpkm, process number is 4.
    saving the result in a txt file named by cell_folder name
    """
    if cell_folder.endswith(os.sep):
        cell_folder = cell_folder[:-1]

    from joblib import Parallel, delayed
    fipath_list = \
        [os.sep.join([cell_folder, fi]) for fi in os.listdir(cell_folder) if fi.endswith('.bed')]
    Parallel(n_jobs=4)(delayed(run1)(fi, refgenes_file) for fi in fipath_list)
    new_fipath_list = [fipath[:-4] + '_rpkm.txt' for fipath in fipath_list]

    result = open(cell_folder + '_rpkms.txt', 'w')
    header = \
        '\t'.join(['Transcipt-ID', 'Gene-Symbol'] +
                  [fi for fi in os.listdir(cell_folder) if fi.endswith('.bed')]) + '\n'
    result.write(header)
    rpkms = [[] for _ in range(len(fipath_list))]
    refgenes = []
    for i, fipath in enumerate(new_fipath_list):
        fi = open(fipath)
        fi.readline()  # skip the header line
        for line in fi:
            rpkms[i].append(line.split()[-1])  # rpkm value
            if i == 0:  # get transcript id and gene symbol form a file
                refgenes.append('\t'.join(line.split('\t')[0:2]))
        fi.close()
    # print rpkms
    for i, refgene in enumerate(refgenes):
        line = '\t'.join([refgene] + [rpkm[i] for rpkm in rpkms]) + '\n'
        result.write(line)
    result.close()

    # delete the temp files
    # for fipath in new_fipath_list:
    #     os.remove(fipath)


if __name__ == '__main__':
    opt_parser = OptionParser()
    opt_parser.add_option('--cell', help='the folder including RNA-seq bed files of cell')
    opt_parser.add_option('--refgenes', help='the ref genes file path')
    options, args = opt_parser.parse_args()
    cell = options.cell
    refgene_file = options.refgenes
    output_rpkm(cell, refgene_file)