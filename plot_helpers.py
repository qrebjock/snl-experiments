import numpy as np

import math

import matplotlib.pyplot as plt

plt.rc("text", usetex=True)
plt.rc("text.latex")
plt.rc("font", family="serif")


def loghist(x, ax, bins=None, **kwargs):
    x = np.asarray(x)
    x = np.where(x < 1e-32, 1e-32, x)

    hist, bin_edges = np.histogram(x, bins=bins)
    logbins = np.logspace(
        np.log10(bin_edges[0]), np.log10(bin_edges[-1]), len(bin_edges)
    )

    ax.hist(x, bins=logbins, alpha=0.85, **kwargs)
    ax.set_xscale("log")

    _, labels = ax.get_legend_handles_labels()
    if labels:
        ax.legend()


def plot_statistics(mats):
    metrics = [
        ("costs", (0, 0), None, "loghist"),
        ("gradnorms", (0, 1), None, "loghist"),
        ("mineigvals", (1, 0), lambda x: np.maximum(-x, 1e-32), "loghist"),
        ("niters", (1, 1), None, "hist"),  # Use linear bins for integer counts
        ("normprocrustdiff", (0, 2), None, "loghist"),
    ]

    for mat in mats:
        print(f"Data size: {mat['costs'].shape}")
        fig, ax = plt.subplots(2, 3, figsize=(10, 5))

        for name, (i, j), transform, method in metrics:
            x = np.ravel(mat[name])
            x_reopt = np.ravel(mat[f"{name}_reopt"])
            if transform:
                x = transform(x)
                x_reopt = transform(x_reopt)

            if method == "loghist":
                loghist(x, ax[i, j], bins=100, label=name, log=True)
                loghist(x_reopt, ax[i, j], bins=100, label=f"{name}_reopt", log=True)
            else:
                ax[i, j].hist(x, bins=50, alpha=0.8, label=name, log=True)
                ax[i, j].hist(
                    x_reopt, bins=20, alpha=0.8, label=f"{name}_reopt", log=True
                )

            ax[i, j].legend()
            ax[i, j].set_title(name)

        plt.tight_layout()
        plt.show()


def compute_connectivity(n, ps):
    con = np.ones((n + 1, ps.size))
    for m in range(1, n + 1):
        for k in range(1, m):
            con[m, :] = con[m, :] - con[k, :] * math.comb(m - 1, k - 1) * (1 - ps) ** (
                k * (m - k)
            )
    return con


def plot_mats(
    mats,
    labels,
    metric1,
    metric2,
    threshold=1e-10,
    figsize=None,
    markers=None,
    markersize=None,
    markevery=None,
    colors=None,
    title=None,
    connectivity=True,
    linewidth=None,
):
    fig, ax = plt.subplots(figsize=figsize)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    for spine in ["bottom", "left"]:
        ax.spines[spine].set_position("zero")
        ax.spines[spine].set_zorder(-1)

    ax.set_ylim([0, 1.07])
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_xticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1"])
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_yticklabels(["", "0.2", "0.4", "0.6", "0.8", "1"])

    ax.plot(1, 0, ">k", markersize=4, clip_on=False, transform=ax.transAxes)
    ax.plot(0, ax.get_ylim()[1], "^k", markersize=4, clip_on=False)

    ax.set_xlabel("Erdős--Rényi density")
    ax.set_ylabel("success rate")

    ax.plot([0, 1], [1, 1], color="black", linewidth=0.75, linestyle="--")

    ps = np.linspace(0, 1, 301)
    if connectivity:
        connectedness = compute_connectivity(mats[0]["n"], ps)
        ax.plot(
            ps,
            connectedness[-1, :],
            linewidth=0.6,
            linestyle="-",
            color="black",
            alpha=0.8,
        )

    for i, mat in enumerate(mats):
        if metric2 is None:
            success = mat[metric1] < threshold * np.sqrt(mat["n"] * mat["dimgt"])
        else:
            success = mat[metric1] < (threshold + mat[metric2]) * (1 + threshold)

        success_rate = np.mean(success, axis=0)
        ax.plot(
            mat["ps"],
            success_rate,
            linewidth=linewidth,
            markevery=markevery,
            marker=markers[i] if markers else None,
            markersize=markersize,
            label=labels[i],
            color=colors[i] if colors else None,
        )
        print(f"Success rate at full density: {success_rate[-1]}")

    ax.legend(frameon=False)

    if title:
        ax.annotate(title, (0.45, 1.08), fontsize=12, annotation_clip=False)

    fig.tight_layout()
    return fig, ax
