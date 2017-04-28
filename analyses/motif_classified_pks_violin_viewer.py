# violin plot of M values distribution of motif exist and motif absent
import os
import numpy
from scipy.stats import gaussian_kde
import pylab as plt

from import_MAmotif import match_manorm_with_motifscan


def half_violin_plot(data, pos, left=False, **kwargs):
    # http://pyinsci.blogspot.it/2009/09/violin-plot-with-matplotlib.html
    # get the value of the parameters
    amplitude = kwargs.pop('amplitude', 0.33)
    # amplitude = kwargs.pop('amplitude', 0.4)
    ax = kwargs.pop('ax', plt.gca())
    # evaluate the violin plot
    x = numpy.linspace(min(data), max(data), 101)  # support for violin
    v = gaussian_kde(data).evaluate(x)  # violin profile (density curve)
    v = v / v.max() * amplitude * (1 if left else -1)  # set the length of the profile
    # kwargs.setdefault('facecolor', 'r')
    kwargs.setdefault('alpha', 0.7)
    return ax.fill_betweenx(x, pos, pos + v, **kwargs)


def violin_plot(data1, motifs, data2, **kwargs):
    ax = kwargs.get('ax', plt.gca())
    ax.spines['left'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    # positions = range(len(data1))
    # positions = [0., 0.8, 1.6, 2.4]
    positions = [i * 0.8 for i in range(len(motifs))]
    assert len(motifs) == len(data1) and len(motifs) == len(data2)
    for pos, key in zip(positions, motifs):
        try:
            d1, d2 = data1[key], data2[key]
        except TypeError:
            d1, d2 = data1[pos], data2[pos]
        half_violin_plot(d1, pos, False, facecolor='r')
        half_violin_plot(d2, pos, True, facecolor='g')
        # division line between the two half
        plt.plot([pos] * 2, [min(min(d1), min(d2)), max(max(d1), max(d2))], 'k-', linewidth=2)

    ax.set_xticks(positions)
    ax.set_xticklabels(motifs)
    ax.set_yticks(range(-5, 15, 5))
    ax.set_yticklabels(['-5', '0', '5', '10'])


# inputs: MAnorm_pk_file, motifscan_file, limit_pvalue
def violin_viewer(pk_file, motifscan_file, motif_pos, limit_pvalue=0.001, neg=False):
    """
    using violin plot to see the M-value distribution difference of motif-exist peaks and motif-absent peaks
    """
    data1 = {}
    data2 = {}

    pk_list, tarnum_dict, motifs = match_manorm_with_motifscan(pk_file, motifscan_file, neg)
    pk_array = numpy.array(pk_list)
    for motif in motifs:
        yes = pk_array[numpy.where(tarnum_dict[motif] > 0)[0]]
        no = pk_array[numpy.where(tarnum_dict[motif] == 0)[0]]
        if motif in motif_pos:
            print motif
            mvalue_yes, mvalue_no = [pk.mvalue for pk in yes], [pk.mvalue for pk in no]
            data1[motif] = mvalue_yes
            data2[motif] = mvalue_no
            print 'yes length: %d' % len(mvalue_yes)
            print 'no length: %d' % len(mvalue_no)

    plt.figure(figsize=(12, 4))
    plt.tick_params(width=2, labelsize=15)
    violin_plot(data1, motif_pos, data2)
    name = pk_file[:-4] + '_violin.png'
    i = 1
    while os.path.exists(name):
        name = name.replace('.png', '%d.png' % i)
        i += 1
    plt.savefig(name)
    plt.close()


def violin_plot2():
    os.chdir('F:\\2.MAmotif\\7. Violin_plot\\1.H1hESC_H3K27ac_Rep2_peaks')
    from statsmodels.graphics.boxplots import violinplot
    import numpy as np

    def read_mvalues(fp):
        m = []
        with open(fp) as fi:
            for line in fi:
                m.append(float(line))
        return m

    m11 = read_mvalues('H1hesc_H3K27ac_Broad_Rep2_peak_MAvalues_absent_Pou5f1.txt')
    m12 = read_mvalues('H1hesc_H3K27ac_Broad_Rep2_peak_MAvalues_exist_Pou5f1.txt')
    m11 = np.array(m11)
    m12 = np.array(m12)

    m21 = read_mvalues('H1hesc_H3K27ac_Broad_Rep2_peak_MAvalues_absent_Sox2.txt')
    m22 = read_mvalues('H1hesc_H3K27ac_Broad_Rep2_peak_MAvalues_exist_Sox2.txt')
    m21 = np.array(m21)
    m22 = np.array(m22)

    m31 = read_mvalues('H1hesc_H3K27ac_Broad_Rep2_peak_MAvalues_absent_Spz1.txt')
    m32 = read_mvalues('H1hesc_H3K27ac_Broad_Rep2_peak_MAvalues_exist_Spz1.txt')
    m31 = np.array(m31)
    m32 = np.array(m32)

    m41 = read_mvalues('H1hesc_H3K27ac_Broad_Rep2_peak_MAvalues_absent_HNF4G.txt')
    m42 = read_mvalues('H1hesc_H3K27ac_Broad_Rep2_peak_MAvalues_exist_HNF4G.txt')
    m41 = np.array(m41)
    m42 = np.array(m42)

    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax.set_ylim([-7, 12])
    ax.set_ylabel('M value')
    # ax.set_xlabel('H1hESC')
    # plt.grid(True)

    labels = ['Pou5f1', 'Sox2', 'Spz1', 'HNF4G']
    m_right = [m11, m21, m31, m41]
    m_left = [m12, m22, m32, m42]
    violinplot(m_right, labels=labels, side='right', ax=ax, show_boxplot=False,
               plot_opts={'violin_fc': 'g'})
    violinplot(m_left, labels=labels, side='left', ax=ax, show_boxplot=False,
               plot_opts={'violin_fc': 'r'})

    ax.plot([1, 1.15], [m11.mean(), m11.mean()], 'g', linewidth=2, alpha=.5)
    ax.plot([1+1, 1.15+1], [m21.mean(), m21.mean()], 'g', linewidth=2, alpha=.5)
    ax.plot([1+1+1, 1.15+1+1], [m31.mean(), m31.mean()], 'g', linewidth=2, alpha=.5)
    ax.plot([1+1+1+1, 1.15+1+1+1], [m41.mean(), m41.mean()], 'g', linewidth=2, alpha=.5)

    ax.plot([1-0.15, 1], [m12.mean(), m12.mean()], 'r', linewidth=2, alpha=.5)
    ax.plot([1-0.15+1, 1+1], [m22.mean(), m22.mean()], 'r', linewidth=2, alpha=.5)
    ax.plot([1-0.15+1+1, 1+1+1], [m32.mean(), m32.mean()], 'r', linewidth=2, alpha=.5)
    ax.plot([1-0.15+1+1+1, 1+1+1+1], [m42.mean(), m42.mean()], 'r', linewidth=2, alpha=.5)

    os.chdir('F:\\2.MAmotif\\7. Violin_plot\\1.K562_H3K27ac_Rep2_peaks')
    m11 = read_mvalues('K562_H3K27ac_Broad_Rep2_Top25K_peak_MAvalues_absent_TAL1__GATA1.txt')
    m12 = read_mvalues('K562_H3K27ac_Broad_Rep2_Top25K_peak_MAvalues_exist_TAL1__GATA1.txt')
    m11 = [-m for m in m11]
    m12 = [-m for m in m12]
    m11 = np.array(m11)
    m12 = np.array(m12)

    m21 = read_mvalues('K562_H3K27ac_Broad_Rep2_Top25K_peak_MAvalues_absent_Gata1.txt')
    m22 = read_mvalues('K562_H3K27ac_Broad_Rep2_Top25K_peak_MAvalues_exist_Gata1.txt')
    m21 = [-m for m in m21]
    m22 = [-m for m in m22]
    m21 = np.array(m21)
    m22 = np.array(m22)

    m31 = read_mvalues('K562_H3K27ac_Broad_Rep2_Top25K_peak_MAvalues_absent_Spz1.txt')
    m32 = read_mvalues('K562_H3K27ac_Broad_Rep2_Top25K_peak_MAvalues_exist_Spz1.txt')
    m31 = [-m for m in m31]
    m32 = [-m for m in m32]
    m31 = np.array(m31)
    m32 = np.array(m32)

    m41 = read_mvalues('K562_H3K27ac_Broad_Rep2_Top25K_peak_MAvalues_absent_HNF4G.txt')
    m42 = read_mvalues('K562_H3K27ac_Broad_Rep2_Top25K_peak_MAvalues_exist_HNF4G.txt')
    m41 = [-m for m in m41]
    m42 = [-m for m in m42]
    m41 = np.array(m41)
    m42 = np.array(m42)

    labels = ['TAL1::GATA1', 'Gata1', 'Spz1', 'HNF4G']
    m_right = [m11, m21, m31, m41]
    m_left = [m12, m22, m32, m42]
    ax = fig.add_subplot(212)
    ax.set_ylim([-7, 12])
    ax.set_ylabel('-M value')
    # ax.set_xlabel('K562')
    # plt.grid(True)
    violinplot(m_right, labels=labels, side='right', ax=ax, show_boxplot=False,
               plot_opts={'violin_fc': 'g'})
    violinplot(m_left, labels=labels, side='left', ax=ax, show_boxplot=False,
               plot_opts={'violin_fc': 'r'})

    ax.plot([1, 1.15], [m11.mean(), m11.mean()], 'g', linewidth=2, alpha=.5)
    ax.plot([1+1, 1.15+1], [m21.mean(), m21.mean()], 'g', linewidth=2, alpha=.5)
    ax.plot([1+1+1, 1.15+1+1], [m31.mean(), m31.mean()], 'g', linewidth=2, alpha=.5)
    ax.plot([1+1+1+1, 1.15+1+1+1], [m41.mean(), m41.mean()], 'g', linewidth=2, alpha=.5)

    ax.plot([1-0.15, 1], [m12.mean(), m12.mean()], 'r', linewidth=2, alpha=.5)
    ax.plot([1-0.15+1, 1+1], [m22.mean(), m22.mean()], 'r', linewidth=2, alpha=.5)
    ax.plot([1-0.15+1+1, 1+1+1], [m32.mean(), m32.mean()], 'r', linewidth=2, alpha=.5)
    ax.plot([1-0.15+1+1+1, 1+1+1+1], [m42.mean(), m42.mean()], 'r', linewidth=2, alpha=.5)

    plt.show()


def test():
    import os
    dir = 'F:\\MAmotif.py\\7. K562_peaks'
    os.chdir(dir)
    motif_pos = ['TAL1::GATA1', 'Gata1', 'Spz1', 'HNF4G']
    # violin_viewer('H1hesc_H3K27ac_Broad_Rep2_peak_MAvalues.xls', 'H3K27ac_H1hesc_Rep2_peak_result', motif_pos)
    violin_viewer('K562_H3K27ac_Broad_Rep2_Top25K_peak_MAvalues.xls', 'K562_H3K27ac_Rep2.pkl', motif_pos, neg=True)
    # violin_viewer('K562_H3K27ac_Broad_Rep2_Top25K_distal_peak_MAvalues.xls', 'K562_H3K27ac_Rep2.pkl', motif_pos,
    #               neg=True)

    # motifs = ['Zfp423', 'En1', 'Pax6', 'RXR::RAR_DR5', 'ESR2', 'Spz1', 'Spi1', 'USF2', 'IRF2', 'RXRA::VDR', 'Nkx3-2',
    #           'Ddit3::Cebpa', 'FOXF2', 'NFIL3', 'JUND_var.2', 'FOXC1', 'TP63', 'DUX4', 'ESR1', 'HNF4G', 'HSF1']
    # motif_pos = ['Pou5f1', 'Sox2', 'Meis1', 'EWSR1-FLI1']
    # i = 0
    # while i < len(motifs):
    #     motif_pos = motifs[i: i + 4]
    #     violin_viewer('H1hesc_H3K27ac_Broad_Rep2_peak_MAvalues.xls', 'H3K27ac_H1hesc_Rep2_peak_result', motif_pos)
    #     violin_viewer('H1hesc_H3K27ac_Broad_Rep2_promoter_peak_MAvalues.xls',
    #                   'H3K27ac_H1hesc_Rep2_peak_result',
    #                   motif_pos)
    #     violin_viewer('H1hesc_H3K27ac_Broad_Rep2_stringent_distal_peak_MAvalues.xls',
    #                   'H3K27ac_H1hesc_Rep2_peak_result', motif_pos)
    #     i += 4


if __name__ == '__main__':
    # test()
    violin_plot2()