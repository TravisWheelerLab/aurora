# Plots distributions of distances between sequences of the same family in the Genome
# Takes a single bed file as an argument.
# Results vary between exponential and power-law depending on the family. Basically all look exponential with fatter tails for a good chunk of families.
#
# Most likely, based on results on understanding the genome (and shape!), the actual distributions are very likely Hyper-Exponential Distributions.
# Basically, it means each family is sampled from one of a weighted set of exponential distributions with different rates.
# This would make sense as TE's have different rates depending on what part of the genome your in (basically, there are functional domains or regions)
#
# Testing this though, would require an implementation of Prony's method to determine the fit. I wasn't able to find any implementations so this
# and it seems implementing such a method would take quite some time. It may be worth eventually adding if better fit's are needed.
#
# For now, an exponential distribution seems to provide a good enough approximation for use in aurora. It also is easy to fit well.
# Here's an interesting paper on the topic that seems to have landed in the same space I've been in: https://www.columbia.edu/~ww2040/FittingMixturesPerfEval98.pdf

bed_file = sys.argv[1]

seq_info = {}
cs = np.inf
ce = -np.inf

with open(bed_file, "r") as f:
    for line in f:
        tokens = line.strip().split()
        start = int(tokens[1])
        end = int(tokens[2])
        join_id = int(tokens[12])
        name = str(tokens[3])

        if name.endswith("Simple_repeat"):
            continue

        name = name.upper()

        if name not in seq_info:
            seq_info[name] = []
        seq_info[name].append((start, end, join_id))

        cs = min(cs, start)
        ce = max(ce, end)


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


distributions = [
    Distribution(expon, (1000.0,)),
    Distribution(weibull_min, (2.0, 1000.0)),
    # Distribution(invweibull, (1.0, 20000.0)),
    Distribution(lomax, (1.0, 1000.0)),
    # Distribution(burr12, (1.0, 1.0, 1000.0)),
    # Distribution(lognorm, (1.0, 1000.0)),
    # Distribution(fisk, (1.0, 1000.0)),
    Distribution(genpareto, (0.0, 1000.0)),
    Distribution(betaprime, (1.0, 2.0, 1000.0)),
]


for name, info in sorted(seq_info.items(), key=lambda k: -len(k[1])):
    # if not name.startswith("CHARLIE1#"):
    #    continue

    info = np.array(info)
    print(info.shape)
    info = info[np.argsort(info[:, 0])]
    starts = info[:, 0]
    ends = info[:, 1]

    # Get distances between consecutive sequences...
    _, indexes = np.unique(info[:, 2], return_index=True)
    indexes = np.sort(indexes)

    dists = np.diff(starts[indexes])

    counts, bins = np.histogram(dists, "auto", density=True)
    x = (bins[1:] + bins[:-1]) / 2

    # Get fit for exponential dist...
    beta = np.mean(dists)  # np.median(dists) / np.log(2)

    # Get fit for weibull dist... (and CDF)...

    valid_values = np.isfinite(wcdfx) & np.isfinite(wcdfy)
    fit_line = curve_fit(
        lambda x, m, b: m * x + b,
        wcdfx[valid_values],
        wcdfy[valid_values],
        sigma=wcdfy_err[valid_values],
    )[0]

    wk = fit_line[0]
    w_beta = np.exp(-fit_line[1] / wk)

    cdf_fits = [
        curve_fit(
            distrib.cdf,
            cdf.quantiles,
            cdf.probabilities,
            distrib.DEFAULTS,
            full_output=True,
        )[0]
        for distrib in distributions
    ]
    # print(cdf_fits)
    # raise ValueError

    # Printout results...
    print(f"Count: {starts.shape[0]}")
    print(f"Genome Size: {ce - cs}")
    est_scale = (ce - cs) / starts.shape[0]
    print(f"Avg Occurance Rate (in nucleotides): {est_scale}")
    print(f"Exponential Fit: Scale (beta): {beta}, Rate (lambda): {1 / beta}")
    print(
        f"Weibull Line Fit: Shape (k): {wk}, Scale (beta): {w_beta}, Rate (lambda): {1 / w_beta}"
    )
    with np.printoptions(suppress=True):
        for distrib, fit in zip(distributions, cdf_fits):
            print(f"{distrib.NAME} CDF Fit: {fit}")

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    fig.suptitle(name)

    ax1.set_title("Histogram and Predicted PDFs")
    ax1.set_ylabel("Density")
    ax1.set_xlabel("Distance (nucleotides)")
    ax1.stairs(counts, bins, fill=True, color="tan", label="Histogram")
    ax1.plot(
        x, weibull_min.pdf(x, wk, 0.0, w_beta), color="orange", label="Weibull Line Fit"
    )
    ax1.plot(x, expon.pdf(x, 0.0, est_scale), color="green", label="Naive Fit Exp")
    for distrib, fit in zip(distributions, cdf_fits):
        ax1.plot(x, distrib.pdf(x, *fit), label=f"{distrib.NAME} CDF Fit")
    ax1.legend()

    ax2.set_title("Weibull-Space of EDF")
    ax2.set_ylabel(r"$ \ln(x) $")
    ax2.set_xlabel(r"$ \ln(-\ln(1 - F(x))) $")
    ax2.scatter(wcdfx, wcdfy, label="Empirical Distribution Function")
    ax2.plot(
        wcdfx,
        fit_line[1] + fit_line[0] * wcdfx,
        label=f"Line: {fit_line[0]:.02f}x + {fit_line[1]:.02f}",
    )
    ax2.legend()

    ax3.set_title("Cumulative Distribution Functions")
    ax3.set_ylabel(r"$ P(X <= x) $")
    ax3.set_xlabel("Distance (nucleotides)")
    ax3.scatter(cdf.quantiles, cdf.probabilities, color="tan", label="EDF")
    ax3.plot(
        cdf.quantiles,
        weibull_min.cdf(cdf.quantiles, wk, 0.0, w_beta),
        color="orange",
        label="Weibull Line Fit",
    )
    ax3.plot(
        cdf.quantiles,
        expon.cdf(cdf.quantiles, 0.0, est_scale),
        color="green",
        label="Naive Fit Exp",
    )
    for distrib, fit in zip(distributions, cdf_fits):
        ax3.plot(
            cdf.quantiles,
            distrib.cdf(cdf.quantiles, *fit),
            label=f"{distrib.NAME} CDF Fit",
        )
    ax3.legend()

    ax4.plot(cdf.quantiles, cdf.probabilities)

    fig.set_size_inches(12, 12)
    fig.tight_layout()
    plt.show()
