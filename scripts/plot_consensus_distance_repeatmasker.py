# This file plots histograms of join distance per tandem repeat family given a repeatmasker formatted bed file as an argument.

import sys

import matplotlib.pyplot as plt
import numpy as np

caf_file = sys.argv[1]

gap_info = {}
length_estimate = {}
prior = {}

with open(caf_file, "r") as f:
    for line in f:
        shared_name = line.strip().split("\t")[3]  # .split("#")[-1]
        seqs = line.strip().split("\t")[-1].split(",")

        if len(seqs) <= 1:
            continue

        # Sort by location on the target...
        seqs = sorted(seqs, key=lambda k: int(k.split(" ")[5]))

        prior = None
        length = 0

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
            length += qend - qstart

            if prior is not None:
                pname, ppos, pstart, pend = prior
                if pname != name:
                    gap = None
                else:
                    try:
                        if is_pos and ppos:
                            assert pstart <= qstart <= qend and pstart <= pend <= qend
                            gap = qstart - pend
                        elif not is_pos and not ppos:
                            assert qstart <= qend <= pend and qstart <= pstart <= pend
                            gap = pstart - qend
                        else:
                            gap = None
                    except AssertionError:
                        if is_pos and ppos:
                            print(
                                f"Violation: {pstart} => {pend} => {qstart} => {qend}"
                            )
                        else:
                            print(
                                f"Violation: {pend} <= {pstart} <= {qend} <= {qstart}"
                            )

                        print(
                            f"Offending annotation: {' '.join(line.split('\t')[:-1])}"
                        )
                        print(f"Offending sequences:\n\t {'\n\t'.join(seqs)}")
                        gap = None

                    if gap is not None:
                        if shared_name not in gap_info:
                            gap_info[shared_name] = []

                        gap_info[shared_name].append(gap)

            prior = (name, is_pos, qstart, qend)

        old_length_est = length_estimate.get(shared_name, (0, 0, 0))
        length_estimate[shared_name] = (
            old_length_est[0] + length,
            old_length_est[1] + 1,
            max(old_length_est[2], length),
        )


for gap_name, gap_vals in sorted(
    gap_info.items(), key=lambda k: len(k[1]), reverse=True
):
    gap_vals = np.array(gap_vals)
    # gap_vals = gap_vals[np.abs(gap_vals) > 1]

    avg_len = length_estimate[gap_name][0] / length_estimate[gap_name][1]
    max_len = length_estimate[gap_name][2]

    plt.title(f"{gap_name} (Avg Length: {avg_len:.02f}, Max Length: {max_len})")

    avg_pos_d = np.mean(gap_vals[gap_vals >= 0])
    avg_neg_d = np.mean(gap_vals[gap_vals <= 0])

    info = f"Positive Avg: {avg_pos_d:.02f}, Negative Avg: {avg_neg_d:.02f}, ({avg_len / avg_pos_d:.02f}, {max_len / avg_pos_d:.02f})"

    range = (0, np.max(gap_vals))
    x = np.arange(*range)
    lamb = 1 / avg_pos_d

    plt.hist(gap_vals, 100, label=info)
    # plt.plot(x, lamb * np.exp(-lamb * x))
    plt.legend()
    plt.show()
