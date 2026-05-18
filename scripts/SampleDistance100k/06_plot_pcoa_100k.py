from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
OUT_DIR = ROOT / "results/sample_distance_100k"
ANN_PATH = OUT_DIR / "subset_100k.h5ad"
EIG_PATH = OUT_DIR / "pcoa_eigenvalues.tsv"
FIG_DIR = OUT_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

COORD_FILES = {
    "bc": OUT_DIR / "pcoa_coords_bc.npy",
    "wunifrac": OUT_DIR / "pcoa_coords_wunifrac.npy",
}
METRIC_TITLES = {
    "bc": "Bray-Curtis",
    "wunifrac": "Weighted UniFrac",
}

BG_COLOR = "#d3d3d3"
BG_SIZE = 0.35
FG_SIZE = 1.2


def as_str(values):
    return values.astype(str).replace({"nan": "NA", "None": "NA", "<NA>": "NA"})


def axis_label(eig: pd.DataFrame, metric: str, axis: int) -> str:
    row = eig[(eig["metric"] == metric) & (eig["axis"] == axis)].iloc[0]
    return f"PC{axis} ({row['explained_trace_ratio'] * 100:.1f}%)"


def axis_limits(coords: np.ndarray) -> tuple[tuple[float, float], tuple[float, float]]:
    x = coords[:, 0]
    y = coords[:, 1]
    x_pad = (x.max() - x.min()) * 0.04
    y_pad = (y.max() - y.min()) * 0.04
    return (float(x.min() - x_pad), float(x.max() + x_pad)), (float(y.min() - y_pad), float(y.max() + y_pad))


def draw_background(ax, coords):
    ax.scatter(coords[:, 0], coords[:, 1], s=BG_SIZE, c=BG_COLOR, alpha=0.22,
               linewidths=0, rasterized=True)


def draw_groups(ax, coords, group_values, order, colors, labels=None):
    handles = []
    for group in order:
        mask = group_values == group
        if mask.sum() == 0:
            continue
        label = labels.get(group, group) if labels else group
        color = colors[group]
        ax.scatter(coords[mask, 0], coords[mask, 1], s=FG_SIZE, c=color,
                   alpha=0.55, linewidths=0, rasterized=True)
        handles.append(Line2D([0], [0], marker="o", color="none",
                              markerfacecolor=color, markeredgewidth=0,
                              markersize=5, label=f"{label} (n={int(mask.sum()):,})"))
    return handles


def finish_axis(ax, metric, eig, xlim, ylim, title):
    ax.set_title(title)
    ax.set_xlabel(axis_label(eig, metric, 1))
    ax.set_ylabel(axis_label(eig, metric, 2))
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)


def human_site_group(site: str) -> str:
    if site in {"gut", "oral", "skin", "urogenital", "lung", "nose", "gastric"}:
        return site
    return "other_or_na"


def rm_class(site: str) -> str:
    if site in {"Negative Control", "Positive Control"}:
        return "control"
    if site in {"NA", "Oral"}:
        return "other_or_na"
    return "respiratory"


print(f"Reading {ANN_PATH.name} ...")
adata = ad.read_h5ad(ANN_PATH, backed="r")
obs = adata.obs[[
    "Database", "sample_role", "overview_bucket", "ma_bucket_detail",
    "human_site", "rm_site", "paired_run_id",
]].copy()
obs = obs.reset_index().rename(columns={"index": "obs_name"})
obs["pos"] = np.arange(len(obs))
for col in ["Database", "sample_role", "overview_bucket", "ma_bucket_detail", "human_site", "rm_site", "paired_run_id"]:
    obs[col] = as_str(obs[col])

coords = {metric: np.load(path) for metric, path in COORD_FILES.items()}
for metric, arr in coords.items():
    assert arr.shape[0] == len(obs) and arr.shape[1] >= 2, f"Bad coords for {metric}: {arr.shape}"

eig = pd.read_csv(EIG_PATH, sep="\t")
limits = {metric: axis_limits(arr) for metric, arr in coords.items()}

ma_env = obs["overview_bucket"].to_numpy()
human_sites = obs["human_site"].map(human_site_group).to_numpy()
rm_groups = obs["rm_site"].map(rm_class).to_numpy()
paired_roles = obs["sample_role"].to_numpy()

ma_env_mask = obs["Database"].eq("MicrobeAtlas").to_numpy()
human_mask = obs["ma_bucket_detail"].eq("Human").to_numpy()
rm_mask = obs["Database"].eq("ResMicroDb").to_numpy()
paired_mask = obs["sample_role"].isin(["paired_ma", "paired_rm"]).to_numpy()

env_values = np.where(ma_env_mask, ma_env, "not_highlighted")
human_values = np.where(human_mask, human_sites, "not_highlighted")
rm_values = np.where(rm_mask, rm_groups, "not_highlighted")
paired_values = np.where(paired_mask, paired_roles, "not_highlighted")

ENV_ORDER = ["Animal", "Soil", "Aquatic", "Plant", "Unknown"]
ENV_COLORS = {
    "Animal": "#2ca25f",
    "Soil": "#8c510a",
    "Aquatic": "#2b8cbe",
    "Plant": "#41b6c4",
    "Unknown": "#636363",
}

HUMAN_ORDER = ["gut", "oral", "skin", "urogenital", "lung", "nose", "gastric", "other_or_na"]
HUMAN_COLORS = {
    "gut": "#d73027",
    "oral": "#1a9850",
    "skin": "#fee08b",
    "urogenital": "#4575b4",
    "lung": "#f46d43",
    "nose": "#984ea3",
    "gastric": "#4dd2d2",
    "other_or_na": "#969696",
}

RM_ORDER = ["respiratory", "control", "other_or_na"]
RM_COLORS = {
    "respiratory": "#fdae61",
    "control": "#542788",
    "other_or_na": "#969696",
}

PAIR_ORDER = ["paired_ma", "paired_rm"]
PAIR_COLORS = {
    "paired_ma": "#1b9e77",
    "paired_rm": "#d95f02",
}
PAIR_LABELS = {
    "paired_ma": "paired MA",
    "paired_rm": "paired RM",
}

ROWS = [
    ("MA environment overview", env_values, ENV_ORDER, ENV_COLORS, None),
    ("MA human body sites", human_values, HUMAN_ORDER, HUMAN_COLORS, None),
    ("ResMicroDb in corpus space", rm_values, RM_ORDER, RM_COLORS, None),
    ("Cross-database paired Runs", paired_values, PAIR_ORDER, PAIR_COLORS, PAIR_LABELS),
]

pair_table = obs.loc[paired_mask, ["paired_run_id", "sample_role", "pos"]].copy()
pair_pivot = pair_table.pivot(index="paired_run_id", columns="sample_role", values="pos")
pair_pivot = pair_pivot.dropna(subset=["paired_ma", "paired_rm"]).astype(int)
line_rng = np.random.default_rng(20260516)
line_n = min(1000, len(pair_pivot))
line_pairs = pair_pivot.iloc[np.sort(line_rng.choice(len(pair_pivot), size=line_n, replace=False))] if line_n else pair_pivot
print(f"paired complete pairs in plot: {len(pair_pivot):,}; drawing connector lines for {len(line_pairs):,}")

fig, axes = plt.subplots(4, 2, figsize=(14, 19), constrained_layout=False)
for row_i, (row_title, values, order, colors, labels) in enumerate(ROWS):
    for col_i, metric in enumerate(["bc", "wunifrac"]):
        ax = axes[row_i, col_i]
        arr = coords[metric]
        xlim, ylim = limits[metric]
        draw_background(ax, arr)
        if row_title.startswith("Cross-database") and len(line_pairs) > 0:
            ma_pos = line_pairs["paired_ma"].to_numpy()
            rm_pos = line_pairs["paired_rm"].to_numpy()
            for p_ma, p_rm in zip(ma_pos, rm_pos):
                ax.plot([arr[p_ma, 0], arr[p_rm, 0]], [arr[p_ma, 1], arr[p_rm, 1]],
                        color="#737373", linewidth=0.25, alpha=0.20, rasterized=True)
        handles = draw_groups(ax, arr, values, order, colors, labels)
        finish_axis(ax, metric, eig, xlim, ylim, f"{row_title} - {METRIC_TITLES[metric]}")
        if col_i == 1:
            ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5),
                      frameon=False, fontsize=8, markerscale=1.2)

fig.suptitle("PCoA of 100,000 representative corpus samples", y=0.995, fontsize=14)
fig.tight_layout(rect=[0, 0, 0.88, 0.985])
out = FIG_DIR / "pcoa_100k_4x2.png"
fig.savefig(out, dpi=180, bbox_inches="tight")
plt.close(fig)
print(f"Wrote {out}")

# Save each question as a compact 1x2 figure for slide/manuscript use.
for row_title, values, order, colors, labels in ROWS:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.6))
    for col_i, metric in enumerate(["bc", "wunifrac"]):
        ax = axes[col_i]
        arr = coords[metric]
        xlim, ylim = limits[metric]
        draw_background(ax, arr)
        if row_title.startswith("Cross-database") and len(line_pairs) > 0:
            for p_ma, p_rm in zip(line_pairs["paired_ma"].to_numpy(), line_pairs["paired_rm"].to_numpy()):
                ax.plot([arr[p_ma, 0], arr[p_rm, 0]], [arr[p_ma, 1], arr[p_rm, 1]],
                        color="#737373", linewidth=0.25, alpha=0.20, rasterized=True)
        handles = draw_groups(ax, arr, values, order, colors, labels)
        finish_axis(ax, metric, eig, xlim, ylim, METRIC_TITLES[metric])
    fig.legend(handles=handles, loc="upper center", ncol=min(len(handles), 4),
               bbox_to_anchor=(0.5, 1.02), frameon=False, fontsize=8, markerscale=1.3)
    fig.suptitle(row_title, y=1.08, fontsize=12)
    fig.tight_layout()
    slug = row_title.lower().replace(" ", "_").replace("-", "").replace("/", "_")
    out = FIG_DIR / f"pcoa_100k_{slug}.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")

# Supplement: RM sample-site detail.
rm_site_values = np.where(rm_mask, obs["rm_site"].to_numpy(), "not_highlighted")
rm_counts = pd.Series(rm_site_values[rm_site_values != "not_highlighted"]).value_counts()
rm_site_order = rm_counts.head(12).index.tolist()
rm_site_values = np.where(
    rm_site_values == "not_highlighted",
    "not_highlighted",
    np.where(np.isin(rm_site_values, rm_site_order), rm_site_values, "Other"),
)
rm_site_order = rm_site_order + ["Other"]
palette = list(plt.cm.tab20.colors)
rm_site_colors = {site: palette[i % len(palette)] for i, site in enumerate(rm_site_order)}
rm_site_colors["Other"] = "#bdbdbd"

fig, axes = plt.subplots(1, 2, figsize=(13, 5.8))
for col_i, metric in enumerate(["bc", "wunifrac"]):
    ax = axes[col_i]
    arr = coords[metric]
    xlim, ylim = limits[metric]
    draw_background(ax, arr)
    handles = draw_groups(ax, arr, rm_site_values, rm_site_order, rm_site_colors)
    finish_axis(ax, metric, eig, xlim, ylim, METRIC_TITLES[metric])
fig.legend(handles=handles, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.04),
           frameon=False, fontsize=7, markerscale=1.4)
fig.suptitle("ResMicroDb sample-site detail", y=1.09, fontsize=12)
fig.tight_layout()
out = FIG_DIR / "pcoa_100k_rm_sample_sites.png"
fig.savefig(out, dpi=180, bbox_inches="tight")
plt.close(fig)
print(f"Wrote {out}")
