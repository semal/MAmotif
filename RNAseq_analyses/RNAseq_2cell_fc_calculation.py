from math import log
from optparse import OptionParser
from RNA_seq_analysis import *


def output_rpkm(cell_folder, refgenes_file):
    """
    using multiple process to get rpkm, process number is 4.
    saving the result in a txt file named by cell_folder name
    """
    from joblib import Parallel, delayed
    fipath_list = [os.sep.join([cell_folder, fi]) for fi in os.listdir(cell_folder) if fi.endswith('.bed')]
    Parallel(n_jobs=4)(delayed(run1)(fi, refgenes_file) for fi in fipath_list)
    new_fipath_list = [fipath[:-4] + '_rpkm.txt' for fipath in fipath_list]
    # print new_fipath_list

    result = open(cell_folder + '_rpkms.txt', 'w')
    header = '\t'.join([''] + [fi for fi in os.listdir(cell_folder) if fi.endswith('.bed')]) + '\n'
    result.write(header)
    rpkms = [[] for _ in range(len(fipath_list))]
    refgenes = []
    for i, fipath in enumerate(new_fipath_list):
        fi = open(fipath)
        fi.readline()  # skip the header line
        for line in fi:
            rpkms[i].append(line.split()[1])
            if i == 0:
                refgenes.append(line.split()[0])
        fi.close()
    # print rpkms
    for i, refgene in enumerate(refgenes):
        line = '\t'.join([refgene] + [rpkm[i] for rpkm in rpkms]) + '\n'
        result.write(line)
    result.close()
    # delete the temp files
    for fipath in new_fipath_list:
        os.remove(fipath)


def calculate_fc(cell1_folder, cell2_folder, refgene_file):
    """
    calculate fold change between 2 cells, cell1 / cell2. The RNA-seq bed file should put in two folders.
    """
    if cell1_folder.endswith(os.sep):
        cell1_folder = cell1_folder[:-1]
    if cell2_folder.endswith(os.sep):
        cell2_folder = cell2_folder[:-1]
    if not os.path.exists(cell1_folder + '_rpkms.txt'):
        output_rpkm(cell1_folder, refgene_file)
    if not os.path.exists(cell2_folder + '_rpkms.txt'):
        output_rpkm(cell2_folder, refgene_file)
    # read the output rpkm result of the 2 folders. and
    cell1_lines = open(cell1_folder + '_rpkms.txt').readlines()
    cell2_lines = open(cell2_folder + '_rpkms.txt').readlines()
    cell1_name = os.path.basename(cell1_folder)
    cell2_name = os.path.basename(cell2_folder)
    name4save = '_'.join([cell1_name, cell2_name])
    fi = open(name4save + '_fold_change.txt', 'w')
    header = \
        '\t'.join([''] + ['_'.join([cell1_name, 'mean']), '_'.join([cell2_name, 'mean']),
                          'fold_change', 'log2_fold_change']) + '\n'
    fi.write(header)
    for l1, l2 in zip(cell1_lines[1:], cell2_lines[1:]):
        refgene = l1.split()[0]
        cell1_rpkm_mean = sum([float(e) for e in l1.split()[1:]]) / len(l1.split()[1:])
        cell2_rpkm_mean = sum([float(e) for e in l2.split()[1:]]) / len(l2.split()[1:])
        try:
            fold_change = cell1_rpkm_mean / cell2_rpkm_mean
            # print fold_change
            log2_fold_change = log(fold_change, 2)
        except:
            fold_change = None
            log2_fold_change = None
        line = '%s\t%f\t%f\t%s\t%s\n' % \
               (refgene, cell1_rpkm_mean, cell2_rpkm_mean, str(fold_change), str(log2_fold_change))
        fi.write(line)
    fi.close()


if __name__ == '__main__':
    opt_parser = OptionParser()
    opt_parser.add_option('--cell1', dest='cell1',
                          help='the folder including RNA-seq bed files of cell 1')
    opt_parser.add_option('--cell2', dest='cell2',
                          help='the folder including RNA-seq bed files of cell 2')
    opt_parser.add_option('--refgenes', dest='refgenes', help='the ref genes file path')
    options, args = opt_parser.parse_args()
    cell1 = options.cell1
    cell2 = options.cell2
    refgene_file = options.refgenes
    calculate_fc(cell1, cell2, refgene_file)