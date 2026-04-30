import sys

import matplotlib.pyplot as plt
import numpy as np

caf_file = sys.argv[1]

gap_info = {}
prior = {}

with open(caf_file, "r") as f:
    for line in f:
        shared_name = line.strip().split("\t")[3].split("#")[-1]
        seqs = line.strip().split("\t")[-1].split(",")

        if len(seqs) <= 1:
            continue

        prior = None

        for seq in seqs:
            tokens = seq.split(" ")
            name = f"{tokens[9]}#{tokens[10]}"
            is_pos = tokens[8] == "+"

            if is_pos:
                qstart = int(tokens[11])
                qend = int(tokens[12])
            else:
                qstart = int(tokens[13])
                qend = int(tokens[12])

            assert qend >= qstart

            if prior is not None:
                pname, ppos, pstart, pend = prior

                try:
                    if is_pos and ppos:
                        assert pstart <= qstart <= qend and pstart <= pend <= qend
                        gap = qstart - pend
                    elif not is_pos and not ppos:
                        assert qstart <= qend <= pstart and qstart <= pstart <= pend
                        gap = pstart - qend
                    else:
                        gap = None
                except AssertionError:
                    if is_pos and ppos:
                        print(f"Violation: {pstart} => {pend} => {qstart} => {qend}")
                    else:
                        print(f"Violation: {pstart} <= {pend} <= {qstart} <= {qend}")

                    print(f"Offending annotation: {' '.join(line.split('\t')[:-1])}")
                    print(f"Offending sequences:\n\t {'\n\t'.join(seqs)}")
                    gap = None

                if gap is not None:
                    if shared_name not in gap_info:
                        gap_info[shared_name] = []

                    gap_info[shared_name].append(gap)

            prior = (name, is_pos, qstart, qend)


for gap_name, gap_vals in sorted(
    gap_info.items(), key=lambda k: len(k[1]), reverse=True
):
    plt.title(gap_name)
    plt.hist(gap_vals, 100)
    plt.show()
