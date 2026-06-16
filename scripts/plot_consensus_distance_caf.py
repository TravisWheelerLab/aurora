import sys

import matplotlib.pyplot as plt
import numpy as np

caf_file = sys.argv[1]

gap_info = {}
prior = {}

with open(caf_file, "r") as f:
    for line in f:
        tokens = line.strip().split(",")
        name = tokens[8].strip()
        qstart = int(tokens[10])
        qend = int(tokens[11])
        tstart = int(tokens[5])
        is_neg_strand = int(tokens[13])

        if name in prior:
            pstart, pend, p_is_neg, p_tstart = prior[name]

            if tstart - p_tstart < 10000:
                if name not in gap_info:
                    gap_info[name] = []

                if not is_neg_strand and not p_is_neg:
                    gap = qstart - pend
                elif is_neg_strand and p_is_neg:
                    gap = qend - pstart
                else:
                    continue

                gap_info[name].append(gap)

        prior[name] = (qstart, qend, is_neg_strand, tstart)


for gap_name, gap_vals in sorted(
    gap_info.items(), key=lambda k: len(k[1]), reverse=True
):
    plt.title(gap_name)
    plt.hist(gap_vals, 500)
    plt.show()
