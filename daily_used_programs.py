# coding=utf-8
# 将从macs到MAmotif的分析自动化，只适用于mm9的组蛋白修饰的数据
import os
import subprocess
import sys
import time
import re
from constant import mm9_refgenes_file as mm9  # 现在适用于mm9
from constant import jaspar3_file as jaspar


macs_fd, manorm_fd, motifscan_fd, mamotif_fd = '1.MACS', '2.MAnorm', '3.MotifScan', '4.MAmotif'
cut_top_peaks_cmd_fp = '/home/zhaohui/bin/macstools/get_top_num_peaks.py'
cut_p_peaks_cmd_fp = '/home/zhaohui/bin/macstools/cut_peaks_with_pvalue.py'
manorm_cmd_fp = 'MAnormFast'
motifscan_cmd_fp = '/mnt/MAmotif/motifscan_pkg/motifscan_pkg/motifscan'
mamotif_cmd_fp = '/home/zhaohui/MAmotif/MAmotif.py'


def make_dir(fd_name):
    try:
        os.mkdir(fd_name)
    except:
        pass


def generate_working_folders():
    """
    生成工作目录所需要的文件夹，程序的输出保存在这些文件夹中
    """
    make_dir(macs_fd)  # 用于保存MACS的输出结果
    make_dir(manorm_fd)  # 用于保存MAnorm的输出结果
    make_dir(motifscan_fd)  # 用于保存MotifScan的输出结果
    make_dir(mamotif_fd)  # 用于保存MAmotif的输出结果


def save_log(log_string, file_name):
    with open(file_name, 'w') as fo:
        fo.write(log_string)


def get_options(common_options):
    opts = ' '.join([' '.join([key, common_options[key]]) for key in common_options.keys()])
    return opts


def run_macs(treatment_fp, control_fp, name):
    """
    运行macs 1.3.7.1 (Oktoberfest, bug fixed #1)
    """
    os.chdir(working_fd)

    def get_tag_size():
        with open(treatment_fp) as fi:
            first_line = fi.readline().split('\t')
            tag_size = int(first_line[2]) - int(first_line[1])
            return tag_size

    common_options = \
        {'-t': treatment_fp,
         '-c': control_fp,
         '--name': name,
         '--nomodel': '',  # nomodel用于跑组蛋白修饰的peak，组蛋白的peak一般比较长
         # '--bw': '500',  # nomodel的时候，macs将bw的值的两倍用于scan window
         '--tsize': str(get_tag_size())}
    cmd = 'macs ' + get_options(common_options) + ' 2> %s_macs_output.log' % name
    os.chdir(macs_fd)
    make_dir(name)  # 将macs的结果输出到以macs的name选项命名的文件夹
    os.chdir(name)
    xls_peak_fn = '_'.join([name, 'peaks.xls'])
    if os.path.exists(xls_peak_fn):
        print '@MACS result already exist.'
        return os.path.abspath(xls_peak_fn)
    subprocess.check_output(cmd, shell=True)
    save_log(cmd, '%s_macs_cmd.log' % name)
    return os.path.abspath(xls_peak_fn)


def cut_macs_peaks(treatment_fp, control_fp, name):
    os.chdir(working_fd)
    cut_top = re.search(r'_Top\d+K', name)
    cut_p = re.search(r'_P\d+', name)
    if not cut_top and not cut_p:  # 没有找到cut peaks的需求, cut Top peaks和cut P peaks只有一种情况会存在
        return run_macs(treatment_fp, control_fp, name)
    elif cut_top:
        cut_top = cut_top.group()
        os.chdir(macs_fd)
        no_cut_fd = name.replace(cut_top, '')
        os.chdir(no_cut_fd)
        cut_peak_fn = '_'.join([name, 'peaks.xls'])
        if os.path.exists(cut_peak_fn):
            return os.path.abspath(cut_peak_fn)
        no_cut_peak_fn = '_'.join([no_cut_fd, 'peaks.xls'])
        cmd = 'python %s ' % cut_top_peaks_cmd_fp + \
              '-s %s ' % no_cut_peak_fn + \
              '-n %s ' % cut_top.replace('_Top', '').replace('K', '') + \
              '-t .'
        print '@cutting peaks ...'
        os.system(cmd)
        return os.path.abspath(cut_peak_fn)
    elif cut_p:
        cut_p = cut_p.group()
        os.chdir(macs_fd)
        no_cut_fd = name.replace(cut_p, '')
        os.chdir(no_cut_fd)
        cut_peak_fn = '_'.join([name, 'peaks.xls'])
        if os.path.exists(cut_peak_fn):
            return os.path.abspath(cut_peak_fn)
        no_cut_peak_fn = '_'.join([no_cut_fd, 'peaks.xls'])
        cmd = 'python %s ' % cut_p_peaks_cmd_fp + \
              '-p %s ' % no_cut_peak_fn + \
              '-v %s' % cut_p.replace('_P', '')
        print '@cutting peaks ...'
        os.system(cmd)
        return os.path.abspath(cut_peak_fn)


def run_manorm(xls_peak_fp1, read_fp1, xls_peak_fp2, read_fp2, comparison_name):
    os.chdir(working_fd)
    os.chdir(manorm_fd)
    common_options = {'--p1': xls_peak_fp1,
                      '--r1': read_fp1,
                      '--p2': xls_peak_fp2,
                      '--r2': read_fp2,
                      '-s': '',
                      '-o': comparison_name}
    cmd = '%s ' % manorm_cmd_fp + get_options(common_options) + \
          ' > %s_MAnorm_output.log' % comparison_name
    print cmd
    try:
        os.chdir(comparison_name)
        manorm_pk_fp1 = os.path.abspath(os.path.basename(xls_peak_fp1).replace('.xls', '_MAvalues.xls'))
        # print manorm_pk_fp1
        manorm_pk_fp2 = os.path.abspath(os.path.basename(xls_peak_fp2).replace('.xls', '_MAvalues.xls'))
        if os.path.exists(manorm_pk_fp1) and os.path.exists(manorm_pk_fp2):
            print '@MAnorm result already exist.'
            return manorm_pk_fp1, manorm_pk_fp2
    except:
        subprocess.call(cmd, shell=True)
        save_log(cmd, '%s_MAnorm_cmd.log' % comparison_name)
        os.chdir(comparison_name)
        manorm_pk_fp1 = os.path.abspath(os.path.basename(xls_peak_fp1).replace('.xls', '_MAvalues.xls'))
        print manorm_pk_fp1
        manorm_pk_fp2 = os.path.abspath(os.path.basename(xls_peak_fp2).replace('.xls', '_MAvalues.xls'))
        print manorm_pk_fp2
        return manorm_pk_fp1, manorm_pk_fp2


def run_motifscan(macs_xls_fp):
    cut_top = re.search(r'_Top\d+K', macs_xls_fp)
    cut_p = re.search(r'_P\d+', macs_xls_fp)
    if cut_top:
        macs_xls_fp = macs_xls_fp.replace(cut_top.group(), '')
    elif cut_p:
        macs_xls_fp = macs_xls_fp.replace(cut_p.group(), '')
    os.chdir(working_fd)
    os.chdir(motifscan_fd)
    common_options = {'-p': macs_xls_fp,
                      '-f': 'macs',
                      '-m': jaspar,
                      '-g': 'mm9',
                      '-t': mm9,
                      '-e': ''}
    cmd = 'python %s ' % motifscan_cmd_fp + get_options(common_options)
    motifscan_pkl_fp = \
        'motifscan_output_' + os.path.basename(macs_xls_fp).replace('.xls', '') + '/peak_result.pkl'
    if os.path.exists(motifscan_pkl_fp):
        print '@MotifScan pkl file already exist.'
        return os.path.abspath(motifscan_pkl_fp)
    os.system(cmd)
    print cmd
    save_log(cmd, macs_xls_fp.replace('.xls', '_motifscan_cmd.log'))
    return os.path.abspath(motifscan_pkl_fp)


def run_mamotif_pipeline(treatment_fp1, control_fp1, name1, treatment_fp2, control_fp2, name2, comparison_name):
    print '#1 running MACS ...'
    # macs_xls_peak_fp1 = run_macs(treatment_fp1, control_fp1, name1)
    # macs_xls_peak_fp2 = run_macs(treatment_fp2, control_fp2, name2)
    macs_xls_peak_fp1, macs_xls_peak_fp2 = \
        map(cut_macs_peaks, [treatment_fp1, treatment_fp2], [control_fp1, control_fp2], [name1, name2])
    print macs_xls_peak_fp1
    print macs_xls_peak_fp2
    print '#2 running MAnorm ...'
    manorm_pk_fp1, manorm_pk_fp2 = \
        run_manorm(macs_xls_peak_fp1, treatment_fp1, macs_xls_peak_fp2, treatment_fp2, comparison_name)
    print '#3 running MotifScan ...'
    # motifscan_pkl_fp1 = run_motifscan(macs_xls_peak_fp1)
    # motifscan_pkl_fp2 = run_motifscan(macs_xls_peak_fp2)
    motifscan_pkl_fp1, motifscan_pkl_fp2 = \
        map(run_motifscan, [macs_xls_peak_fp1, macs_xls_peak_fp2])
    time.sleep(2)
    print '#4 running MAmotif ...'
    os.chdir(working_fd)
    os.chdir(mamotif_fd)
    cmd1 = 'python %s ' % mamotif_cmd_fp + \
           '-p %s ' % manorm_pk_fp1 + \
           '-M %s ' % motifscan_pkl_fp1 + \
           '-r %s' % mm9
    os.system(cmd1)
    cmd2 = 'python %s ' % mamotif_cmd_fp + \
           '-p %s ' % manorm_pk_fp2 + \
           '-M %s ' % motifscan_pkl_fp2 + \
           '-r %s ' % mm9 + '-n'
    # 因为都在同一个文件夹下，避免同时创建文件夹时出现死锁的情况不使用多进程同时运行cmd1和cmd2
    os.system(cmd2)


def test_run_mamotif_pipeline():
    # os.chdir('/mnt/MAmotif/2.Processing/7.primed_vs_naive/New_version')
    # working_fd = '/mnt/MAmotif/2.Processing/7.primed_vs_naive/New_version'
    read_fp1 = '/mnt/MAmotif/1.RAWdata/Mixed_primed_naive_hESC_public_mm9/GSE59434_primed_vs_naive_MIT/' \
               'hESC_WIBR2_primed_H3K4me3_MIT.bed'
    input_fp1 = '/mnt/MAmotif/1.RAWdata/Mixed_primed_naive_hESC_public_mm9/GSE59434_primed_vs_naive_MIT/' \
                'hESC_WIBR2_primed_H3K4me3_MIT_input.bed'
    name1 = 'hESC_WIBR2_primed_H3K4me3_MIT_P100'
    read_fp2 = '/mnt/MAmotif/1.RAWdata/Mixed_primed_naive_hESC_public_mm9/GSE59434_primed_vs_naive_MIT/' \
               'hESC_WIBR2_6iLA_naive_H3K4me3_MIT.bed'
    input_fp2 = '/mnt/MAmotif/1.RAWdata/Mixed_primed_naive_hESC_public_mm9/GSE59434_primed_vs_naive_MIT/' \
                'hESC_WIBR2_6iLA_naive_H3K4me3_MIT_input.bed'
    name2 = 'hESC_WIBR2_naive_H3K4me3_MIT_P100'
    comparison_name = 'H3K4me3_hESC_WIBR2_primed_vs_naive_MIT_P100'
    run_mamotif_pipeline(read_fp1, input_fp1, name1, read_fp2, input_fp2, name2, comparison_name)


def cmd():
    if not os.path.exists(macs_fd):
        yn = raw_input('Current working directory do not exist requisite folders, '
                       'do you want to create them?(yes or no)')
        if yn == 'yes':
            generate_working_folders()
        else:
            sys.exit(0)
    from optparse import OptionParser
    optparser = OptionParser()
    optparser.add_option('--t1', dest='t1', help='treatment bed file1')
    optparser.add_option('--c1', dest='c1', help='control bed file1')
    optparser.add_option('--n1', dest='n1', help='name1')
    optparser.add_option('--t2', dest='t2', help='treatment bed file2')
    optparser.add_option('--c2', dest='c2', help='control bed file2')
    optparser.add_option('--n2', dest='n2', help='name2')
    optparser.add_option('--name', dest='name', help='comparison name')
    options, args = optparser.parse_args()
    t1, c1, n1 = os.path.abspath(options.t1), os.path.abspath(options.c1), options.n1
    t2, c2, n2 = os.path.abspath(options.t2), os.path.abspath(options.c2), options.n2
    name = options.name
    run_mamotif_pipeline(t1, c1, n1, t2, c2, n2, name)


if __name__ == '__main__':
    # jaspar = '/mnt/MAmotif/motifscan_pkg/input_data/motif_list'
    working_fd = os.path.abspath(os.getcwd())
    # working_fd = '/mnt/MAmotif/2.Processing/7.primed_vs_naive/New_version'
    cmd()
    # test_run_mamotif_pipeline()