# coding=utf-8
# 调用DEseq程序。输入是两个样本的read count的文件。这两个read count文件由calculate_read_count.py生成。
import os
import sys


file_dir = os.path.dirname(os.path.realpath(__file__))


def run_DESeq(read_count_file1, read_count_file2):
    """
    call R script for running DESeq
    @param read_count_file1: using as treatment of DESeq
    @param read_count_file2: using as control of DESeq
    """
    # to get treatment and control sample's number
    treatment_num = get_sample_num(read_count_file1)
    control_num = get_sample_num(read_count_file2)

    # put 2 read count file together as input of DESeq R script input
    input_file_path = put_2together(read_count_file2, read_count_file1)

    # generate R script
    # print file_dir
    template_script_file_path = os.sep.join([file_dir, 'run_parrwise_deseq.r'])
    template_handle = open(template_script_file_path)
    cnt = template_handle.read()
    template_handle.close()
    new_script_name = 'temp.r'
    handle = open(new_script_name, 'w')
    new_cnt = \
        cnt.replace('file_path', input_file_path[:-4]) \
            .replace('control_num', str(control_num)) \
            .replace('treatment_num', str(treatment_num))
    # print new_cnt
    handle.write(new_cnt)
    handle.close()

    # run R script
    print 'Start running R script ...'
    os.system('Rscript %s' % new_script_name)
    os.remove(new_script_name)


def get_sample_num(read_count_file):
    f = open(read_count_file)
    sample_num = len(f.readline().split('\t')) - 1
    print sample_num
    f.close()
    return sample_num


def put_2together(file1, file2):
    """
    put 2 read count output file together with same format
    """
    file1_name = os.path.basename(file1)[:-4]
    file2_name = os.path.basename(file2)[:-4]
    file_new = '_'.join([file1_name, file2_name]) + '.txt'
    fi = open(file_new, 'w')

    gene_names = []
    fi1, fi2 = open(file1), open(file2)
    for i1, i2 in zip(fi1.readlines(), fi2.readlines()):
        i2_cnt = i2.split('\t')[1:]
        line = '\t'.join([i1.strip()] + i2_cnt)
        gene_name = i2.split('\t')[0].strip()
        if gene_name not in gene_names:
            fi.write(line)
            gene_names.append(gene_name)
    fi1.close(), fi2.close()

    return os.sep.join([os.path.abspath(os.curdir), file_new])


if __name__ == '__main__':
    if '-h' in sys.argv:
        print 'python DEseq.py treatment_read_count_file control_read_count_file'
        exit(1)
    treatment = sys.argv[1]
    control = sys.argv[2]
    run_DESeq(treatment, control)