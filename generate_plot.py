from functools import lru_cache

def _PoleGap(K, L, X):
    assert K % 2 == 0
    return 3 * K * L // 2 + K // 2 + 3 * X - 2

def PoleGap(K, L, X):
    if K % 2 == 0 and L % 2 == 0:
        return min(_PoleGap(K, L, X), _PoleGap(L, K, X))
    elif K % 2 == 0:
        return _PoleGap(K, L, X)
    elif L % 2 == 0:
        return _PoleGap(L, K, X)
    raise

@lru_cache
def GASP_r(K, L, X, r=None):
    if r is None:
        return [GASP_r(K, L, X, r) for r in range(1, min(K, X) + 1)]

    f_degrees_1 = list(range(K))
    f_degrees_2 = []
    i = 0
    while len(f_degrees_2) < X:
        f_degrees_2 += [K * L + i * K + j for j in range(r)]
        i += 1
    f_degrees_2 = f_degrees_2[:X]

    g_degrees_1 = list(range(0, K * L, K))
    g_degrees_2 = [K * L + i for i in range(X)]

    f_degrees = f_degrees_1 + f_degrees_2
    g_degrees = g_degrees_1 + g_degrees_2

    return len(set(f + g for f in f_degrees for g in g_degrees))

def _GASP_small(K, L, X):
    return GASP_r(K, L, X, 1)

def GASP_small(K, L, X):
    return min(_GASP_small(K, L, X), _GASP_small(L, K, X))

def _GASP_big(K, L, X):
    return GASP_r(K, L, X, min(K, X))

def GASP_big(K, L, X):
    return min(_GASP_big(K, L, X), _GASP_big(L, K, X))

def _GASP(K, L, X):
    return min(GASP_r(K, L, X, r) for r in range(1, min(K, X) + 1))

def GASP(K, L, X):
    return min(_GASP(K, L, X), _GASP(L, K, X))

def GASP_old(K, L, X):
    return min(GASP_small(K, L, X), GASP_big(K, L, X))

def _A3S(K, L, X):
    return (K + 1) * (L + X) - 1

def A3S(K, L, X):
    return min(_A3S(K, L, X), _A3S(L, K, X))

if __name__ == "__main__":

    import matplotlib as mpl
    mpl.use("pgf")

    import matplotlib.pyplot as plt
    import numpy as np

    plt.rcParams.update(
        {
            "font.family": "serif",  # use serif/main font for text elements
            "text.usetex": True,  # use inline math for ticks
            "pgf.rcfonts": False,  # don't setup fonts from rc parameters
            "pgf.preamble": r"\usepackage{amsmath, amssymb}",  # use ams packages
        }
    )

    fig, (ax1, ax2) = plt.subplots(
        1, 2, sharey=False, tight_layout=True, frameon=False, dpi=200.0, figsize=(7, 3)
    )

    # Plot of rate as a function of X, with fixed K = L = 14

    K = L = 14
    Xs = range(1, 81)

    rate_A3S        = [K * L / A3S(K, L, X)        for X in Xs]
    rate_GASP_big   = [K * L / GASP_big(K, L, X)   for X in Xs]
    rate_GASP_small = [K * L / GASP_small(K, L, X) for X in Xs]
    rate_GASP       = [K * L / GASP(K, L, X)       for X in Xs]
    rate_PoleGap    = [K * L / PoleGap(K, L, X)    for X in Xs]

    print("Largest increase in rate:", (np.array(rate_PoleGap) / np.array(rate_GASP)).max() - 1)

    # ax1.plot(Xs, rate_A3S,        linestyle=":",  color="tab:pink",   label="A3S")
    ax1.plot(Xs, rate_GASP_big,   linestyle="-.", color="tab:blue",   label=r"GASP\textsubscript{big}")
    ax1.plot(Xs, rate_GASP_small, linestyle="-.", color="tab:orange", label=r"GASP\textsubscript{small}")
    ax1.plot(Xs, rate_GASP,       linestyle="-.", color="tab:green",  label="GASP")
    ax1.plot(Xs, rate_PoleGap,    linestyle="-",  color="tab:red",    label="PoleGap (Thm. IV.2)")
    ax1.legend()
    ax1.grid()
    ax1.set_xlabel("$X$")
    ax1.set_ylabel(r"$\mathcal{R}$")

    del K, L, Xs, rate_A3S, rate_GASP_big, rate_GASP_small, rate_GASP, rate_PoleGap

    # Plot as a function of K, with K = L and fixed X = 50

    X = 50
    Ks = list(range(2, 28, 2))  # we only do even K

    rate_A3S        = [K * K / A3S(K, K, X)        for K in Ks]
    rate_GASP_big   = [K * K / GASP_big(K, K, X)   for K in Ks]
    rate_GASP_small = [K * K / GASP_small(K, K, X) for K in Ks]
    rate_GASP       = [K * K / GASP(K, K, X)       for K in Ks]
    rate_PoleGap    = [K * K / PoleGap(K, K, X)    for K in Ks]

    print("Largest increase in rate:", (np.array(rate_PoleGap) / np.array(rate_GASP)).max() - 1)

    # ax2.plot(Ks, rate_A3S,        linestyle=":",  color="tab:pink",   label="A3S")
    ax2.plot(Ks, rate_GASP_big,   linestyle="-.", color="tab:blue",   label=r"GASP\textsubscript{big}")
    ax2.plot(Ks, rate_GASP_small, linestyle="-.", color="tab:orange", label=r"GASP\textsubscript{small}")
    ax2.plot(Ks, rate_GASP,       linestyle="-.", color="tab:green",  label="GASP")
    ax2.plot(Ks, rate_PoleGap,    linestyle="-",  color="tab:red",    label="PoleGap (Thm. IV.2)")
    # ax2.legend()  # do not include legend for this one
    ax2.grid()
    ax2.set_xlabel(r"$K = L$")
    ax2.set_ylabel(r"$\mathcal{R}$")

    fig.savefig("plot.pdf")
