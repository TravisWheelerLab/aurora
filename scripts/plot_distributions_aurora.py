import sys
import typing
from dataclasses import dataclass, fields

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import ecdf, expon, genpareto, gumbel_r, invweibull, norm, weibull_min


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


aurora_file = sys.argv[1]

joined_annots: dict[tuple[int, int], list[AuroraEntry]] = {}

with open(aurora_file, "r") as f:
    for line in f:
        entry = AuroraEntry.from_line(line)
        key = (entry.region, entry.join_id)
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
            return best(
                a.query_start[a_idx] - b.query_start[b_idx] - 1,
                b.query_end[b_idx] - a.query_end[a_idx] - 1,
            )
        case ("+", "-"):
            return best(
                b.query_start[b_idx] - a.query_start[a_idx] - 1,
                a.query_end[a_idx] - b.query_end[b_idx] - 1,
            )
        case _:
            raise ValueError("Unknown stand configuration.")


def _target_distance(a: AuroraEntry, b: AuroraEntry, ai: int, bi: int):
    assert b.target_start[bi] - a.target_end[ai] - 1 >= 0
    return b.target_start[bi] - a.target_end[ai] - 1


stats_to_compute = {
    "Consensus Distance": consensus_dist,
    "Target Distance": _target_distance,
    "Divergence Change": lambda a, b, ai, bi: b.kimera80[bi] - a.kimera80[ai],
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
    "Consensus Distance": Distribution(invweibull, (1.0, 0.0, 1.0), False),
    "Target Distance": Distribution(
        weibull_min, (1.0, 10000.0)
    ),  # Distribution(genpareto, (0.0, 1.0)),
    "Divergence Change": Distribution(norm, (0.0, 1.0), False),
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


for name in stats_to_compute:
    join_stats[name] = {}
    random_stats[name] = {}

for k, annots in sorted(
    joined_annots.items(), key=lambda k: min(min(v.target_start) for v in k[1])
):
    for ann in annots[:1]:
        for i in range(len(ann.query)):
            name = ann.query[i]
            (pann, j) = prior_vals.get(name, (None, None))
            if pann is not None:
                for stat_name, values_per_query in random_stats.items():
                    (ann1, idx1), (ann2, idx2) = sorted(
                        [(pann, j), (ann, i)], key=lambda v: v[0].target_start[v[1]]
                    )

                    if name not in values_per_query:
                        values_per_query[name] = []
                    values_per_query[name].append(
                        stats_to_compute[stat_name](ann1, ann2, idx1, idx2)
                    )

            seq_size[name] = max(
                seq_size.get(name, 1),
                ann.query_start[i],
                ann.query_end[i],
            )
            prior_vals[name] = (ann, i)

    if len(annots) <= 1:
        continue

    for a, b in zip(annots[:-1], annots[1:]):
        for i in range(len(a.query)):
            name = a.query[i]
            for stat_name, values_per_query in join_stats.items():
                if name not in values_per_query:
                    values_per_query[name] = []
                values_per_query[name].append(stats_to_compute[stat_name](a, b, i, i))


for query_name, _ in sorted(
    join_stats["Consensus Distance"].items(), key=lambda k: -len(k[1])
):
    fig, axs = plt.subplots(3, len(stats_to_compute))
    axs = axs.T

    fig.suptitle(f"{query_name} (Size: {seq_size.get(query_name, 0)})")

    for name, (ax1, ax2, ax3) in zip(stats_to_compute, axs):
        est = estimator[name]

        sx = np.sort(join_stats[name][query_name])
        fit = fit_dist(sx, est)
        ax1.set_title(f"Join {name}")
        ax1.hist(
            join_stats[name][query_name],
            50,
            density=True,
            label=f"Mean: {np.mean(join_stats[name][query_name]):.02f}\nSTD: {np.std(join_stats[name][query_name]):.02f}",
        )
        ax1.plot(
            sx, est.pdf(sx, *fit), label=f"Fit: {', '.join(f'{v:.02f}' for v in fit)}"
        )
        ax1.legend()

        sx2 = np.sort(random_stats[name][query_name])
        fit2 = fit_dist(sx2, est)
        ax2.set_title(f"All {name}")
        ax2.hist(
            random_stats[name][query_name],
            50,
            density=True,
            label=f"Mean: {np.mean(random_stats[name][query_name]):.02f}\nSTD: {np.std(random_stats[name][query_name]):.02f}",
        )
        ax2.plot(
            sx2,
            est.pdf(sx2, *fit2),
            label=f"Fit: {', '.join(f'{v:.02f}' for v in fit2)}",
        )
        ax2.legend()

        ax3.set_title("CDFs")
        ax3.ecdf(join_stats[name][query_name], label="Joins CDF")
        ax3.ecdf(random_stats[name][query_name], label="All CDF")
        ax3.plot(sx, est.cdf(sx, *fit), label="Est. Join CDF")
        ax3.plot(sx2, est.cdf(sx2, *fit2), label="Est. All CDF")
        ax3.legend()

    fig.set_size_inches(12, 8)
    fig.tight_layout()
    plt.show()
