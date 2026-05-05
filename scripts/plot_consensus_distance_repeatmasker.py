# This file plots histograms of join distance per tandem repeat family given a repeatmasker formatted bed file as an argument.

import sys

import matplotlib.pyplot as plt
import numpy as np

caf_file = sys.argv[1]

gap_info = {}
length_estimate = {}
sequences = []


def get_gap(pstart, pend, qstart, qend, ppos, is_pos, check_valid: bool = True):
    try:
        if is_pos and ppos:
            if check_valid:
                assert pstart <= qstart <= qend and pstart <= pend <= qend
            gap = qstart - pend
        elif not is_pos and not ppos:
            if check_valid:
                assert qstart <= qend <= pend and qstart <= pstart <= pend
            gap = pstart - qend
        else:
            gap = None
    except AssertionError:
        if is_pos and ppos:
            print(f"Violation: {pstart} => {pend} => {qstart} => {qend}")
        else:
            print(f"Violation: {pend} <= {pstart} <= {qend} <= {qstart}")

        print(f"Offending annotation: {' '.join(line.split('\t')[:-1])}")
        print(f"Offending sequences:\n\t {'\n\t'.join(seqs)}")
        gap = None

    return gap


with open(caf_file, "r") as f:
    for line in f:
        shared_name = line.strip().split("\t")[3]  # .split("#")[-1]
        seqs = line.strip().split("\t")[-1].split(",")
        sequences.extend(seqs)

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
            join_id = int(tokens[14])

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
                    gap = get_gap(pstart, pend, qstart, qend, ppos, is_pos)

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


sequences.sort(key=lambda seq: int(seq.split(" ")[5]))
gap_nojoin_info = {}
other_priors = {}
target_nojoin_info = {}

for seq in sequences:
    tokens = seq.split(" ")
    name = f"{tokens[9]}#{tokens[10]}"
    is_pos = tokens[8] == "+"
    join_id = int(tokens[14])
    tstart = int(tokens[5])

    if is_pos:
        qstart = int(tokens[11])
        qend = int(tokens[12])
    else:
        qstart = int(tokens[13])
        qend = int(tokens[12])

    random_prior = other_priors.get(name, None)
    if random_prior is not None:
        ppos, pstart, pend, pjoin_id, p_tstart = random_prior
        if pjoin_id == join_id:
            gap = None
        else:
            gap = get_gap(pstart, pend, qstart, qend, ppos, is_pos, False)

        target_gap = tstart - p_tstart

        if gap is not None:
            if name not in target_nojoin_info:
                target_nojoin_info[name] = []
            target_nojoin_info[name].append(target_gap)

            if name not in gap_nojoin_info:
                gap_nojoin_info[name] = []
            gap_nojoin_info[name].append(gap)

    other_priors[name] = (is_pos, qstart, qend, join_id, tstart)


for gap_name, gap_join_vals in sorted(
    gap_info.items(), key=lambda k: len(k[1]), reverse=True
):
    gap_join_vals = np.array(gap_join_vals)
    # gap_vals = gap_vals[np.abs(gap_vals) > 1]

    avg_len = length_estimate[gap_name][0] / length_estimate[gap_name][1]
    max_len = length_estimate[gap_name][2]

    avg_pos_d = np.mean(gap_join_vals[gap_join_vals >= 0])
    avg_neg_d = np.mean(gap_join_vals[gap_join_vals <= 0])

    info = f"Positive Avg: {avg_pos_d:.02f},\n Negative Avg: {avg_neg_d:.02f}, ({avg_len / avg_pos_d:.02f}, {max_len / avg_pos_d:.02f})"

    range = (0, np.max(gap_join_vals))
    x = np.arange(*range)
    lamb = 1 / avg_pos_d

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, squeeze=True)
    ax1.set_title(f"{gap_name} (Avg Length: {avg_len:.02f}, Max Length: {max_len})")

    ax1.hist(gap_join_vals, 50, label=info)
    # plt.plot(x, lamb * np.exp(-lamb * x))
    ax1.legend()

    ax2.set_title(f"{gap_name} Non-Join Consensus Gaps")
    ax2.hist(gap_nojoin_info[gap_name], 50)

    ax3.set_title(f"{gap_name} Scatterplot")
    ax3.scatter(target_nojoin_info[gap_name], gap_nojoin_info[gap_name])
    ax3.set_xlabel("Target Distance")
    ax3.set_ylabel("Consensus Distance")

    w, h = fig.get_size_inches()
    fig.set_size_inches(w * 3, h)
    fig.tight_layout()
    plt.show()
