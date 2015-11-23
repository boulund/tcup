from ete3 import NCBITaxa
from collections import defaultdict, OrderedDict
import os
import sys

from ..proteotyping.ncbi import NCBITaxa_mod as NCBITaxa


def LCA_test():
    print("LCA (4):", NCBITaxa.lowest_common_ancestor(((590,543,91347,1236,1224,2,131567,1),
                                (314322,570,543,91347,1236,1224,2,131567,1),
                                (314322,530,543,91347,1236,1224,2,131567,1),
                                (314323,530,533,91347,1236,1224,2,131567,1))))
    print("LCA (3):", NCBITaxa.lowest_common_ancestor([[590,543,913447,12236,12244,25,1321567,14],
                                [314322,570,5433,913247,12336,12424,23,13153567,13],
                                [314322,530,5523,913447,1236,12244,21,131567,12]]))
    print("_common_lineage (3):", NCBITaxa._common_lineage(((590,543,91347,1236,1224,2,131567,1),
                                                        (314322,570,543,91347,1236,1224,2,131567,1),
                                                        (314322,530,553,91347,1236,1224,2,131567,1))))
    try:
        print(NCBITaxa.lowest_common_ancestor((314323,530,533,91347,1236,1224,2,131567,1)))
    except TypeError as e:
        print("Caught TypeError:", e)
