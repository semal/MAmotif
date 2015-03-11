from import_MAmotif import motif_targenes_de_fisher_test


def command():
    from optparse import OptionParser
    optparser = OptionParser()
    optparser.add_option('-p', dest='pk', help='MAnorm peak file path')
    optparser.add_option('-M', dest='ms', help='Motifscan target number file path')
    optparser.add_option('-r', dest='ref', help='RefSeq gene file path')
    optparser.add_option('-d', dest='de', help='Different expression gene list file path')
    options, args = optparser.parse_args()
    pk = options.pk
    ms = options.ms
    ref = options.ref
    de = options.de
    motif_targenes_de_fisher_test(pk, ms, de, ref)


if __name__ == '__main__':
    command()