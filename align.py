# 第一步：实现read回帖，直接根据回帖的位置得到eccDNA的全长？
# 第二步：实现trim的功能，走CCS这一条流程是否还需要此功能？如果要trim的话也是在CCS之前进行trim

import os
import sys
import logging
from collections import defaultdict
LOGGER = logging.getLogger('ECCTRY')

import pysam
# from ECCTRY import env
# from ECCTRY.utils import *

LOGGER = logging.getLogger('ECCTRY')

# 去掉N：3，最后做统一的替换 [todo]
# OPERATION = {
#     'M': 0,
#     'I': 1,
#     'D': 2,
#     'S': 4,
#     'H': 5,
#     'P': 6,
#     '=': 7,
#     'X': 8,
#     0: 'M',
#     1: 'I',
#     2: 'D',
#     4: 'S',
#     5: 'H',
#     6: 'P',
#     7: '=',
#     8: 'X',
# }


class Faidx(object):
    """
    API for unify pysam::Fastafile and mappy2::Aligner
    """

    def __init__(self,infile):
        if not os.path.exists(infile + '.fai'):
            LOGGER.warning('Index of reference genome not found, generating ...')
        try:
            faidx = pysam.FastaFile(infile)
        except ValueError:
            LOGGER.error('Cannot generate index of reference genome, index file is missing')
            sys.exit()
        except IOError:
            LOGGER.error('Cannot generate index of reference genome, file could not be opened')
            sys.exit()
    
        self.faidx = faidx
        self.contig_len = {contig: faidx.get_reference_length(contig) for contig in faidx.references}

    def seq(self, contig, start, end):
        return self.faidx.fetch(contig, start, end)

    def close(self):
        self.faidx.close()


# class Hit(object):

#     def __init__(self, aln):
#         self.ctg = aln.rname
#         self.strand = 1 if aln.orient == '+' else -1
#         self.cigar_string = aln.cigar
#         self.r_st = aln.pos
#         self.r_en, self.q_st, self.q_en, self.blen = self.__parse_cigar()
#         self.mlen = self.q_en - self.q_st
#         self.is_primary = aln.is_primary
#         self.mapq = aln.mapq

#     @property
#     def cigar(self):
#         return [(int(length), OPERATION[operation]) for length, operation in
#                 re.findall(r'(\d+)([MIDNSHP=X])', self.cigar_string)]

#     def __parse_cigar(self):
#         r_en = self.r_st
#         q_st, q_en = 0, 0
#         blen = 0
#         for length, operation in self.cigar:
#             if operation == 0:
#                 q_en += length
#                 r_en += length
#                 blen += length
#             elif operation == 1:
#                 q_en += length
#             elif operation in [2]:
#                 r_en += length
#                 blen += length
#             elif operation in [4, 5]:
#                 if q_st == 0:
#                     q_st = length
#                     q_en = length
#             else:
#                 pass
#         return r_en, q_st, q_en, blen

#     def __str__(self):
#         return '\t'.join([str(x) for x in [self.q_st, self.q_en, self.ctg, self.r_st, self.r_en, self.mlen, self.blen,
#                                            self.cigar_string]])


def get_primary_alignment(hits):
    if not hits:
        return []  
        
    primary_hits = []
    for hit in hits:
        if hit.is_primary:
            primary_hits.append(hit)
    
    primary_hits_sort = sorted(primary_hits, key=lambda x: x.q_st)

    return primary_hits_sort

