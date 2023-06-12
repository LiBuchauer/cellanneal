import numpy as np
from matplotlib.pyplot import savefig, subplots, subplots_adjust
from matplotlib import rcParams, cycler
from seaborn import heatmap
from scipy.stats import spearmanr
from scipy.spatial.distance import correlation


rcParams["axes.prop_cycle"] = cycler(
    color=[
        "#e6194b",
        "#3cb44b",
        "#ffe119",
        "#4363d8",
        "#f58231",
        "#911eb4",
        "#46f0f0",
        "#f032e6",
        "#bcf60c",
        "#fabebe",
        "#008080",
        "#e6beff",
        "#9a6324",
        "#fffac8",
        "#800000",
        "#aaffc3",
        "#808000",
        "#ffd8b1",
        "#000075",
        "#808080",
    ]
)


def plot_pies(mix_df, save_path=None):
    # for each mixture, plot a pie chart
    plot_num = np.shape(mix_df)[0] + 1
    fig, axes = subplots(
        int(np.ceil(plot_num / 4)), 4, figsize=(12, 3 * np.ceil(plot_num / 4))
    )

    # if correlation values are included in the mix_df,
    # they should not be part of the pie
    corr_list = []
    if "rho_Spearman" in mix_df.columns:
        corr_list.append("rho_Spearman")
    if "rho_Pearson" in mix_df.columns:
        corr_list.append("rho_Pearson")
    if len(corr_list) > 0:
        plot_df = mix_df.drop(corr_list, axis=1)
    else:
        plot_df = mix_df

    # onto each axis, plot a mixture pie chart
    plotted_a_pie = 0
    for i, ax in enumerate(axes.flatten()):
        try:
            mixture = plot_df.iloc[i]
            pie = ax.pie(mixture)
            ax.set_title(plot_df.index.tolist()[i])
            plotted_a_pie = i
        except (IndexError, ValueError):
            ax.set_visible(False)
    # onto the next empty axis after the plotting stopped, plot a legend
    lax = axes.flatten()[plotted_a_pie + 1]
    lax.set_visible(True)
    lax.axis("off")
    lax.legend(
        pie[0], plot_df.columns, loc="lower center", ncol=2, bbox_to_anchor=(0.7, 0)
    )
    if save_path is not None:
        savefig(save_path, bbox_inches="tight")


def plot_mix_heatmap(mix_df, rownorm=False, save_path=None):
    fig, ax = subplots(
        figsize=(
            10 / np.shape(mix_df)[1] * np.shape(mix_df)[0],
            0.5 * np.shape(mix_df)[1],
        )
    )
    # if correlation values are included in the mix_df,
    # they should not be part of the pie
    corr_list = []
    if "rho_Spearman" in mix_df.columns:
        corr_list.append("rho_Spearman")
    if "rho_Pearson" in mix_df.columns:
        corr_list.append("rho_Pearson")
    if len(corr_list) > 0:
        plot_df = mix_df.drop(corr_list, axis=1)
    else:
        plot_df = mix_df

    if rownorm:
        tdf = plot_df.T.sort_index(axis=0)
        ax = heatmap(
            tdf.div(tdf.max(axis=1), axis=0),
            linewidths=0.5,
            square=True,
            cmap="viridis",
            cbar_kws={"shrink": 0.7},
        )
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    else:
        ax = heatmap(
            plot_df.T.sort_index(axis=0),
            linewidths=0.5,
            square=True,
            cmap="viridis",
            cbar_kws={"shrink": 0.7, "label": "fraction"},
        )
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

    bottom, top = ax.get_ylim()
    bottom, top = ax.get_ylim()
    ax.set_ylim(bottom + 0.5, top - 0.5)
    ax.invert_yaxis()

    if save_path is not None:
        savefig(save_path, bbox_inches="tight")


def plot_mix_heatmap_log(mix_df, rownorm=False, save_path=None):
    fig, ax = subplots(
        figsize=(
            10 / np.shape(mix_df)[1] * np.shape(mix_df)[0],
            0.5 * np.shape(mix_df)[1],
        )
    )
    # if correlation values are included in the mix_df,
    # they should not be part of the pie
    corr_list = []
    if "rho_Spearman" in mix_df.columns:
        corr_list.append("rho_Spearman")
    if "rho_Pearson" in mix_df.columns:
        corr_list.append("rho_Pearson")
    if len(corr_list) > 0:
        plot_df = mix_df.drop(corr_list, axis=1)
    else:
        plot_df = mix_df

    # for log version, add epsilon and take log of everything
    plot_df = plot_df.add(1e-4)
    for col in plot_df.columns:
        plot_df[col] = np.log10(plot_df[col])
    if rownorm:
        tdf = plot_df.T.sort_index(axis=0)
        ax = heatmap(
            tdf.div(tdf.max(axis=1), axis=0),
            linewidths=0.5,
            square=True,
            cmap="viridis",
            cbar_kws={"shrink": 0.7},
        )
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    else:
        ax = heatmap(
            plot_df.T.sort_index(axis=0),
            linewidths=0.5,
            square=True,
            cmap="viridis",
            cbar_kws={"shrink": 0.7, "label": "log10 of fraction (truncated at 1e-4)"},
        )
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

    bottom, top = ax.get_ylim()
    bottom, top = ax.get_ylim()
    ax.set_ylim(bottom + 0.5, top - 0.5)
    ax.invert_yaxis()

    if save_path is not None:
        savefig(save_path, bbox_inches="tight")


def plot_1D_lines(mix_df, save_path=None):
    # if correlation values are included in the mix_df,
    # they should not be part of the pie
    corr_list = []
    if "rho_Spearman" in mix_df.columns:
        corr_list.append("rho_Spearman")
        spear_list = np.round(mix_df["rho_Spearman"].values, 3)
    else:
        spear_list = ["-" for i in mix_df.index]
    if "rho_Pearson" in mix_df.columns:
        corr_list.append("rho_Pearson")
        pear_list = np.round(mix_df["rho_Pearson"].values, 3)
    else:
        pear_list = ["-" for i in mix_df.index]
    if len(corr_list) > 0:
        plot_df = mix_df.drop(corr_list, axis=1)
    else:
        plot_df = mix_df

    # preparations: count number of plots required, number of lines and columns
    plot_num = np.shape(plot_df)[1]
    xticks = np.arange(np.shape(plot_df)[0])
    col_num = 5
    row_num = int(np.ceil(plot_num / col_num))

    g = plot_df.plot(
        rot=90,
        subplots=True,
        layout=(row_num, col_num),
        figsize=(col_num * 3, row_num * 2.5),
        xticks=xticks,
        marker="o",
    )

    # all y_axes should start at 0
    for plot in g.flatten():
        plot.set_ylim(bottom=0)

    # adjust spacing
    subplots_adjust(hspace=0.2, wspace=0.4)
    # save
    if save_path is not None:
        savefig(save_path, bbox_inches="tight")


# function for pie plots from one lcm position set of results
def plot_scatter(mix_df, bulk_df, celltype_df, gene_dict, save_path=None):
    # if correlation values are included in the mix_df, remove
    corr_list = []
    if "rho_Spearman" in mix_df.columns:
        corr_list.append("rho_Spearman")
        spear_list = np.round(mix_df["rho_Spearman"].values, 3)
    else:
        spear_list = ["-" for i in mix_df.index]
    if "rho_Pearson" in mix_df.columns:
        corr_list.append("rho_Pearson")
        pear_list = np.round(mix_df["rho_Pearson"].values, 3)
    else:
        pear_list = ["-" for i in mix_df.index]
    if len(corr_list) > 0:
        plot_df = mix_df.drop(corr_list, axis=1)
    else:
        plot_df = mix_df

    # get gene restricted subsets of bulk and sc data
    sc_list = []
    bulk_comp_list = []  # compositional version of subset bulk data

    for b, bulk in enumerate(bulk_df.columns):
        # first, subset bulk data
        bulk_sub = bulk_df[bulk].loc[gene_dict[bulk]].values
        bulk_sub = bulk_sub / np.sum(bulk_sub)
        bulk_comp_list.append(bulk_sub)

        # next, subset sc data
        sc_sub = celltype_df.loc[gene_dict[bulk]].values
        sc_list.append(sc_sub)

    # for each mixture, plot a scatterplot of mixed vs real bulk
    plot_num = len(bulk_df.columns) + 1

    fig, axes = subplots(
        int(np.ceil(plot_num / 4)),
        4,
        figsize=(12, 3 * np.ceil(plot_num / 4)),
        sharex=False,
        sharey=True,
    )

    # onto each axis, plot a scatter after calculating mixture vector
    for i, ax in enumerate(axes.flatten()):
        try:
            mixed_counts = np.dot(plot_df.iloc[i], sc_list[i].T).T
            mixed_compositional = mixed_counts / mixed_counts.sum()
            comp_vec = bulk_comp_list[i]

            # make compositional
            mixed_compositional = mixed_counts / mixed_counts.sum(axis=0)
            ax.scatter(
                comp_vec,
                mixed_compositional,
                alpha=0.1,
                color="darkblue",
                rasterized=True,
            )

            # calc spearman and pearson
            spearman = spearmanr(mixed_compositional, comp_vec)[0]
            pearson = 1 - correlation(mixed_compositional, comp_vec)

            ax.set_xlabel("bulk data")
            ax.set_ylabel("best sc mix")
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim([5e-7, 0.02])
            ax.set_ylim([5e-7, 0.02])
            ax.set_title(
                bulk_df.columns.tolist()[i]
                + "\n P_corr={}, S_corr={}".format(
                    np.round(pearson, 2), np.round(spearman, 2)
                )
            )

        except (IndexError, ValueError):
            ax.set_visible(False)

    fig.tight_layout()
    if save_path is not None:
        savefig(save_path, bbox_inches="tight")
