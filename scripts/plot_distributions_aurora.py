import sys
import typing
from dataclasses import dataclass, fields
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgba
from scipy.optimize import curve_fit, minimize
from scipy.stats import (
    ecdf,
    expon,
    genpareto,
    gumbel_r,
    halfnorm,
    invweibull,
    laplace_asymmetric,
    norm,
    weibull_min,
)


@dataclass
class AuroraEntry:
    annotation_count: int
    target: str
    target_start: list[int]
    target_end: list[int]
    query: list[str]
    query_start: list[int]
    query_end: list[int]
    strand: list[str]
    score: float
    kimera80: list[float]
    join_id: int
    region: int

    @classmethod
    def from_line(cls, line: str) -> typing.Self:
        parts = line.strip().split()

        parsed = []
        max_list_len = 0

        for field, part in zip(fields(cls), parts):
            if isinstance(field.type, str):
                raise ValueError("Can't parse strings!")

            origin = typing.get_origin(field.type)
            if origin is not None and issubclass(origin, list):
                tp = typing.get_args(field.type)[0]
                parsed.append([tp(v) for v in part.split(",")])
                max_list_len = max(max_list_len, len(parsed[-1]))
            else:
                parsed.append(field.type(part))

        if max_list_len > 0:
            parsed_new = []

            for v in parsed:
                if isinstance(v, list):
                    if len(v) == max_list_len:
                        parsed_new.append(v)
                    elif len(v) == 1:
                        parsed_new.append(v * max_list_len)
                    else:
                        raise ValueError(
                            f"One of the elements has and invalid length: {v}"
                        )
                else:
                    parsed_new.append(v)

            parsed = parsed_new

        return cls(*parsed)

    def get_key(self) -> tuple[int, int]:
        return (self.region, self.join_id)

    def select(self, idx: int) -> typing.Self:
        return type(self)(
            self.annotation_count,
            self.target,
            [self.target_start[idx]],
            [self.target_end[idx]],
            [self.query[idx]],
            [self.query_start[idx]],
            [self.query_end[idx]],
            [self.strand[idx]],
            self.score,
            [self.kimera80[idx]],
            self.join_id,
            self.region,
        )


if len(sys.argv) not in [2, 3]:
    print("Usage:")
    print(f"\t{Path(sys.argv[0]).name} AURORA_OUTPUT_FILE [dist|scatter]")
    sys.exit(1)

aurora_file = sys.argv[1]
mode = sys.argv[2] if len(sys.argv) > 2 else "dist"

if mode not in ["scatter", "dist"]:
    print("Second argument (the mode) must be 'scatter' or 'dist'!")
    sys.exit(1)

if mode == "dist":
    print("Generating distributions plots...")
else:
    print("Generating scatter plots...")

joined_annots: dict[tuple[int, int], list[AuroraEntry]] = {}

with open(aurora_file, "r") as f:
    for line in f:
        entry = AuroraEntry.from_line(line)
        key = entry.get_key()
        if key not in joined_annots:
            joined_annots[key] = []
        joined_annots[key].append(entry)

for key, annots in joined_annots.items():
    annots.sort(key=lambda k: min(k.target_start))


random_stats = {}
join_stats = {}
prior_vals = {}
seq_size = {}


def consensus_dist(a: AuroraEntry, b: AuroraEntry, a_idx: int, b_idx: int):
    def best(*a):
        return min(a, key=lambda v: abs(v))

    match (a.strand[a_idx], b.strand[b_idx]):
        case ("+", "+"):
            return b.query_start[b_idx] - a.query_end[a_idx] - 1
        case ("-", "-"):
            return a.query_end[a_idx] - b.query_start[b_idx] - 1
        case ("-", "+"):
            return None
            return best(
                a.query_start[a_idx] - b.query_start[b_idx] - 1,
                b.query_end[b_idx] - a.query_end[a_idx] - 1,
            )
        case ("+", "-"):
            return None
            return best(
                b.query_start[b_idx] - a.query_start[a_idx] - 1,
                a.query_end[a_idx] - b.query_end[b_idx] - 1,
            )
        case _:
            raise ValueError("Unknown stand configuration.")


def _target_distance(a: AuroraEntry, b: AuroraEntry, ai: int, bi: int):
    assert b.target_start[bi] - a.target_end[ai] - 1 >= 0
    return b.target_start[bi] - a.target_end[ai] - 1


def _kimura_dist(a, b, ai, bi):
    return abs(b.kimera80[bi] - a.kimera80[ai])


def _relative_consensus_dist(a, b, ai, bi):
    d = consensus_dist(a, b, ai, bi)
    sum_seq_len = max(
        [
            abs(v)
            for v in (
                a.query_end[ai] - a.query_start[ai],
                b.query_end[bi] - b.query_start[bi],
            )
        ]
    )
    if d is None:
        return None
    try:
        return d / sum_seq_len
    except ZeroDivisionError:
        return None


stats_to_compute = {
    "Target Distance": _target_distance,
    "Divergence Change": _kimura_dist,
    "Relative Consensus Distance": _relative_consensus_dist,
}


class Distribution:
    def __init__(self, dist, defaults, exclude_location=True):
        self._dist = dist
        self.DEFAULTS = list(defaults)
        self._excl_loc = exclude_location
        self.NAME = getattr(
            dist, "__name__", getattr(type(dist), "__qualname__", repr(dist))
        )

    def pdf(self, x, *args):
        # print(*args)
        if self._excl_loc:
            return self._dist.pdf(x, *args[:-1], 0.0, args[-1])
        else:
            return self._dist.pdf(x, *args)

    def cdf(self, x, *args):
        if self._excl_loc:
            return self._dist.cdf(x, *args[:-1], 0.0, args[-1])
        else:
            return self._dist.cdf(x, *args)

    def logcdf(self, x, *args):
        print(*args)
        if self._excl_loc:
            return self._dist.logcdf(x, *args[:-1], 0.0, args[-1])
        else:
            return self._dist.logcdf(x, *args)


estimator = {
    "Relative Consensus Distance": Distribution(
        invweibull, (1.0, 0.0, 1.0), False
    ),  # Distribution(invweibull, (1.0, 0.0, 1.0), False), Distribution(laplace_asymmetric, (1.0, 0.0, 1.0), False)
    "Target Distance": Distribution(
        genpareto, (0.0, 1.0)
    ),  # Distribution(expon, (1.0,)),  # Distribution(genpareto, (0.0, 1.0)), Distribution(weibull_min, (1.0, 10000)
    "Divergence Change": Distribution(halfnorm, (1.0,)),
}


def fit_dist(data, dist):
    emp_cdf = ecdf(data).cdf

    return curve_fit(
        dist.cdf,
        emp_cdf.quantiles,
        emp_cdf.probabilities,
        dist.DEFAULTS,
        full_output=True,
    )[0]


random_idx_reference = {}
random_is_join = {}

for name in stats_to_compute:
    join_stats[name] = {}
    random_stats[name] = {}

all_anots_flat = [ann for annots in joined_annots.values() for ann in annots]
all_anots_flat.sort(key=lambda v: min(v.target_start))


for ann_i, ann in enumerate(all_anots_flat):
    for i in range(len(ann.query)):
        name = ann.query[i]
        (pann, j, pann_i) = prior_vals.get(name, (None, None, None))
        if pann is not None:
            # if pann.region == ann.region and pann.join_id == ann.join_id:
            #    continue

            (ann1, idx1), (ann2, idx2) = sorted(
                [(pann, j), (ann, i)], key=lambda v: v[0].target_start[v[1]]
            )

            stats = {
                stat_name: stats_to_compute[stat_name](ann1, ann2, idx1, idx2)
                for stat_name in stats_to_compute
            }

            if all([v is not None for v in stats.values()]):
                for stat_name, value in stats.items():
                    values_per_query = random_stats[stat_name]

                    if name not in values_per_query:
                        values_per_query[name] = []
                    values_per_query[name].append(value)

                if name not in random_idx_reference:
                    random_idx_reference[name] = []
                    random_is_join[name] = []

                random_idx_reference[name].append((ann_i, i, pann_i, j))
                random_is_join[name].append(
                    ann in joined_annots.get(pann.get_key(), [])
                )

        seq_size[name] = max(
            seq_size.get(name, 1),
            ann.query_start[i],
            ann.query_end[i],
        )
        prior_vals[name] = (ann, i, ann_i)

for k, annots in joined_annots.items():
    if len(annots) <= 1:
        continue

    for a, b in zip(annots[:-1], annots[1:]):
        for i in range(len(a.query)):
            name = a.query[i]

            stats = {
                stat_name: stats_to_compute[stat_name](a, b, i, i)
                for stat_name in stats_to_compute
            }

            if all([v is not None for v in stats.values()]):
                for stat_name, value in stats.items():
                    values_per_query = join_stats[stat_name]

                    if name not in values_per_query:
                        values_per_query[name] = []
                    values_per_query[name].append(value)


for query_name, _ in sorted(
    next(iter(join_stats.values())).items(), key=lambda k: -len(k[1])
):
    # if not query_name.startswith("sin"):
    #    continue
    join_indexes = np.flatnonzero(random_is_join[query_name])
    not_join_indexes = np.flatnonzero(~np.array(random_is_join[query_name]))

    if mode == "dist":
        fig, axs = plt.subplots(3, len(stats_to_compute))
        axs = axs.T

        fig.suptitle(f"{query_name} (Size: {seq_size.get(query_name, 0)})")

        for name, (ax1, ax2, ax3) in zip(stats_to_compute, axs):
            est = estimator[name]

            join_samples = np.array(join_stats[name][query_name])
            sx = np.linspace(join_samples.min(), join_samples.max(), 1000)
            fit = fit_dist(join_samples, est)
            ax1.set_title(f"Join {name}")
            ax1.hist(
                join_samples,
                200,
                density=True,
                label=f"Mean: {np.mean(join_samples):.02f}\nSTD: {np.std(join_samples):.02f}",
            )
            ax1.plot(
                sx,
                est.pdf(sx, *fit),
                label=f"Fit: {', '.join(f'{v:.02f}' for v in fit)}",
            )
            ax1.legend(fontsize="xx-small")

            random_samples = np.array(random_stats[name][query_name])
            sx2 = np.linspace(random_samples.min(), random_samples.max(), 1000)
            fit2 = fit_dist(random_samples[not_join_indexes], est)
            ax2.set_title(f"All {name}")
            ax2.plot(
                [0],
                [0],
                color="black",
                visible=False,
                label=f"Mean: {np.mean(random_samples):.02f}\nSTD: {np.std(random_samples):.02f}",
            )
            ax2.hist(
                [random_samples[not_join_indexes], random_samples[join_indexes]],
                200,
                label=["Not Joined", "Joined"],
                density=True,
                stacked=True,
                color=[to_rgba(c) for c in ["tab:orange", "tab:blue"]],
            )
            ax2.plot(
                sx2,
                est.pdf(sx2, *fit2),
                "black",
                label=f"Fit: {', '.join(f'{v:.02f}' for v in fit2)}",
            )
            ax2.legend(fontsize="xx-small")

            ax3.set_title("CDFs")
            ax3.ecdf(join_samples, label="Joins CDF")
            ax3.ecdf(random_samples[not_join_indexes], label="No Joins CDF")
            ax3.plot(sx, est.cdf(sx, *fit), label="Est. Join CDF")
            ax3.plot(sx2, est.cdf(sx2, *fit2), label="Est. No Join CDF")
            ax3.legend(fontsize="xx-small")

        fig.set_size_inches(16, 8)
        fig.tight_layout()
        plt.show()
    else:
        plt.title(f"{query_name} (Size: {seq_size.get(query_name, 0)})")

        join_art = plt.plot(
            np.array(random_stats["Target Distance"][query_name])[join_indexes],
            np.array(random_stats["Consensus Distance"][query_name])[join_indexes],
            "ro",
            picker=5,
            label="Joins",
        )
        no_join_art = plt.plot(
            np.array(random_stats["Target Distance"][query_name])[not_join_indexes],
            np.array(random_stats["Consensus Distance"][query_name])[not_join_indexes],
            "bo",
            picker=5,
            label="Not Joins",
        )
        plt.xlabel("Target Distance")
        plt.ylabel("Consensus Distance")
        plt.legend()
        fig = plt.gcf()

        def on_pick(evt):
            mask = join_indexes if evt.artist == join_art else not_join_indexes

            for idx in evt.ind:
                idx = mask[idx]
                annot_idx, sub_i, pann_idx, p_sub_i = random_idx_reference[query_name][
                    idx
                ]
                print(f"Index: {annot_idx}, Sub-Index: {sub_i}")
                print(
                    f"\tTarget Distance: {random_stats['Target Distance'][query_name][idx]}"
                )
                print(
                    f"\tConsensus Distance: {random_stats['Consensus Distance'][query_name][idx]}"
                )
                print(f"\tPrior: {all_anots_flat[pann_idx].select(p_sub_i)}")
                print(f"\tCurrent: {all_anots_flat[annot_idx].select(sub_i)}")
                print(f"\tIs Joined: {random_is_join[query_name][idx]}")

        fig.canvas.mpl_connect("pick_event", on_pick)
        plt.show()
