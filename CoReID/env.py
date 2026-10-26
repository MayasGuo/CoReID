ALIGNER = None
CONTIG_LEN = None
GENOME = None



def initializer(aligner, contig_len, genome):
    global ALIGNER
    global CONTIG_LEN
    global GENOME
    global STRAND_INDEX

    ALIGNER = aligner
    CONTIG_LEN = contig_len
    GENOME = genome

