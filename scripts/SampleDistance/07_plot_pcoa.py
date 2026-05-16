# ---
# jupyter:
#   jupytext:
#     formats: py:percent,ipynb
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: baseBio
#     language: python
#     name: python3
# ---

# %% [markdown]
# # 07: PCoA 散点可视化
#
# 输入：
# - `subset_50k.h5ad`：`obsm['X_pcoa_bc'/'X_pcoa_wunifrac']` + `obs['stratum_id']`
# - `pcoa_eigenvalues.tsv`：解释方差比
#
# 输出：
# - `results/sample_distance/figures/pcoa_bucket.png`    BC + wUniFrac 各 PC1×PC2 / PC1×PC3
# - `results/sample_distance/figures/pcoa_scree.png`     特征值 scree 图

# %%
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib as mpl

ROOT = Path("/hpcdisk1/limk_group/caiqy/project/260428_greengene2")
ANN_PATH = ROOT / "results/sample_distance/subset_50k.h5ad"
EIG_PATH = ROOT / "results/sample_distance/pcoa_eigenvalues.tsv"
OUT_DIR  = ROOT / "results/sample_distance/figures"
OUT_DIR.mkdir(exist_ok=True)

mpl.rcParams.update({
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "legend.fontsize": 8,
    "savefig.dpi": 150,
})

# %% [markdown]
# ## §1 读数据，把 stratum_id 收成大类着色

# %%
print(f"读 {ANN_PATH.name} ...")
adata = ad.read_h5ad(ANN_PATH, backed="r")
print(f"  shape: {adata.shape}")
pc_bc = np.asarray(adata.obsm["X_pcoa_bc"])
pc_wu = np.asarray(adata.obsm["X_pcoa_wunifrac"])
print(f"  pcoa_bc: {pc_bc.shape}  pcoa_wunifrac: {pc_wu.shape}")

stratum = adata.obs["stratum_id"].astype(str).values
print(f"\nstratum_id 取值（共 {len(set(stratum))} 个）")

# 收成 7 个大类
def collapse(s):
    if s.startswith("MA::"):
        return s.split("::", 1)[1]   # Human / Animal_other / Soil / Aquatic / Plant / Unknown
    return "RM_Respiratory"          # 所有 RM:: 折成一类（RM 几乎全是呼吸道）

bucket = np.array([collapse(s) for s in stratum])
print("\nbucket 分布:")
print(pd.Series(bucket).value_counts().to_string())

# %% [markdown]
# ## §2 颜色 & 绘图顺序（小类后画，盖在大类上）

# %%
# 大类按数量从大到小排序，小类后画
BUCKET_ORDER = ["Animal_other", "Human", "Soil", "Aquatic", "RM_Respiratory",
                "Plant", "Unknown"]
COLORS = {
    "Human":          "#e6194B",  # 红
    "Animal_other":   "#3cb44b",  # 绿
    "Soil":           "#9A6324",  # 棕
    "Aquatic":        "#4363d8",  # 蓝
    "Plant":          "#42d4f4",  # 青
    "RM_Respiratory": "#f58231",  # 橙
    "Unknown":        "#a9a9a9",  # 灰
}

# 读特征值用于坐标轴标签
eig = pd.read_csv(EIG_PATH, sep="\t").set_index("axis")
def axlabel(metric: str, axis: int) -> str:
    col = f"{metric}_explained_var"
    pct = eig.loc[axis, col] * 100
    return f"PC{axis} ({pct:.1f}%)"


def scatter_panel(ax, coords, axes_pair, title):
    i, j = axes_pair  # 1-based
    for b in BUCKET_ORDER:
        m = bucket == b
        if m.sum() == 0:
            continue
        ax.scatter(coords[m, i - 1], coords[m, j - 1],
                   s=1.5, alpha=0.4, c=COLORS[b], label=f"{b} (n={int(m.sum())})",
                   linewidths=0, rasterized=True)
    ax.set_title(title)
    ax.set_xlabel(axlabel(title.split(' ')[0].lower().replace("wunifrac", "wunifrac").replace("bc", "bc"), i))
    ax.set_ylabel(axlabel(title.split(' ')[0].lower(), j))


# %% [markdown]
# ## §3 主图：2×2  (BC PC12/PC13)  ×  (wUniFrac PC12/PC13)

# %%
fig, axes = plt.subplots(2, 2, figsize=(13, 11))

# BC
for b in BUCKET_ORDER:
    m = bucket == b
    if m.sum() == 0:
        continue
    axes[0, 0].scatter(pc_bc[m, 0], pc_bc[m, 1], s=1.5, alpha=0.4,
                       c=COLORS[b], label=f"{b} (n={int(m.sum()):,})",
                       linewidths=0, rasterized=True)
    axes[0, 1].scatter(pc_bc[m, 0], pc_bc[m, 2], s=1.5, alpha=0.4,
                       c=COLORS[b], linewidths=0, rasterized=True)
    axes[1, 0].scatter(pc_wu[m, 0], pc_wu[m, 1], s=1.5, alpha=0.4,
                       c=COLORS[b], linewidths=0, rasterized=True)
    axes[1, 1].scatter(pc_wu[m, 0], pc_wu[m, 2], s=1.5, alpha=0.4,
                       c=COLORS[b], linewidths=0, rasterized=True)

axes[0, 0].set_title("Bray-Curtis: PC1 vs PC2")
axes[0, 0].set_xlabel(axlabel("bc", 1)); axes[0, 0].set_ylabel(axlabel("bc", 2))
axes[0, 1].set_title("Bray-Curtis: PC1 vs PC3")
axes[0, 1].set_xlabel(axlabel("bc", 1)); axes[0, 1].set_ylabel(axlabel("bc", 3))

axes[1, 0].set_title("Weighted UniFrac: PC1 vs PC2")
axes[1, 0].set_xlabel(axlabel("wunifrac", 1)); axes[1, 0].set_ylabel(axlabel("wunifrac", 2))
axes[1, 1].set_title("Weighted UniFrac: PC1 vs PC3")
axes[1, 1].set_xlabel(axlabel("wunifrac", 1)); axes[1, 1].set_ylabel(axlabel("wunifrac", 3))

# 共用图例放在外面
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", ncol=4,
           bbox_to_anchor=(0.5, 0.99), frameon=False, markerscale=4)
fig.suptitle(f"PCoA on 50,000-sample subset (colored by environment bucket)", y=1.02)
fig.tight_layout(rect=[0, 0, 1, 0.95])

png = OUT_DIR / "pcoa_bucket.png"
fig.savefig(png, dpi=150, bbox_inches="tight")
print(f"已写出 {png} ({png.stat().st_size/1024:.0f} KB)")
plt.close(fig)

# %% [markdown]
# ## §4 Scree 图

# %%
fig, ax = plt.subplots(1, 1, figsize=(6, 4))
x = eig.index.values
ax.plot(x, eig["bc_explained_var"] * 100, "o-", label="Bray-Curtis", color="#e6194B")
ax.plot(x, eig["wunifrac_explained_var"] * 100, "s-", label="Weighted UniFrac", color="#4363d8")
ax.set_xlabel("PCoA axis")
ax.set_ylabel("Explained variance (%)")
ax.set_title("PCoA scree plot")
ax.set_xticks(x)
ax.grid(alpha=0.3)
ax.legend()
fig.tight_layout()
png = OUT_DIR / "pcoa_scree.png"
fig.savefig(png, dpi=150, bbox_inches="tight")
print(f"已写出 {png} ({png.stat().st_size/1024:.0f} KB)")
plt.close(fig)

# %% [markdown]
# ## §5 额外：MA Human 内部按 sample site 二级着色（独立小图）

# %%
ma_human_mask = (bucket == "Human")
print(f"\nMA Human 样本: {ma_human_mask.sum():,}")
sub_strat = adata.obs["sub_stratum"].astype(str).values[ma_human_mask]
print("MA Human sub_stratum:")
print(pd.Series(sub_strat).value_counts().head(10).to_string())

SITE_COLORS = {
    "gut":         "#e6194B",
    "oral":        "#3cb44b",
    "skin":        "#ffe119",
    "urogenital":  "#4363d8",
    "lung":        "#f58231",
    "nose":        "#911eb4",
    "gastric":     "#46f0f0",
    "bone":        "#f032e6",
    "NA":          "#a9a9a9",
    "nan":         "#a9a9a9",
}
fig, axes = plt.subplots(1, 2, figsize=(13, 6))
for site, col in SITE_COLORS.items():
    m = ma_human_mask & (adata.obs["sub_stratum"].astype(str).values == site)
    if m.sum() == 0:
        continue
    axes[0].scatter(pc_bc[m, 0], pc_bc[m, 1], s=1.5, alpha=0.5, c=col,
                    label=f"{site} (n={int(m.sum())})", linewidths=0, rasterized=True)
    axes[1].scatter(pc_wu[m, 0], pc_wu[m, 1], s=1.5, alpha=0.5, c=col,
                    linewidths=0, rasterized=True)
axes[0].set_title(f"MA Human only · Bray-Curtis PC1 vs PC2 (n={ma_human_mask.sum():,})")
axes[0].set_xlabel(axlabel("bc", 1)); axes[0].set_ylabel(axlabel("bc", 2))
axes[1].set_title(f"MA Human only · Weighted UniFrac PC1 vs PC2")
axes[1].set_xlabel(axlabel("wunifrac", 1)); axes[1].set_ylabel(axlabel("wunifrac", 2))
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", ncol=5, frameon=False,
           bbox_to_anchor=(0.5, 1.02), markerscale=4)
fig.tight_layout(rect=[0, 0, 1, 0.92])

png = OUT_DIR / "pcoa_human_sites.png"
fig.savefig(png, dpi=150, bbox_inches="tight")
print(f"已写出 {png} ({png.stat().st_size/1024:.0f} KB)")
plt.close(fig)

# %% [markdown]
# ## §6 额外：RM 内部按 sample site 二级着色

# %%
rm_mask = (bucket == "RM_Respiratory")
print(f"\nRM 样本: {rm_mask.sum():,}")
rm_sites = adata.obs["RM_Sample_Site"].astype(str).values
rm_top = pd.Series(rm_sites[rm_mask]).value_counts()
print("RM_Sample_Site top:")
print(rm_top.head(15).to_string())

# 取 top-12，其它合 Other
top12 = rm_top.head(12).index.tolist()
palette = plt.cm.tab20.colors
fig, axes = plt.subplots(1, 2, figsize=(13, 6))
for k, site in enumerate(top12 + ["Other"]):
    if site == "Other":
        m = rm_mask & (~np.isin(rm_sites, top12))
        col = "#cccccc"
    else:
        m = rm_mask & (rm_sites == site)
        col = palette[k % len(palette)]
    if m.sum() == 0:
        continue
    axes[0].scatter(pc_bc[m, 0], pc_bc[m, 1], s=1.5, alpha=0.5, c=col,
                    label=f"{site} (n={int(m.sum())})", linewidths=0, rasterized=True)
    axes[1].scatter(pc_wu[m, 0], pc_wu[m, 1], s=1.5, alpha=0.5, c=col,
                    linewidths=0, rasterized=True)
axes[0].set_title(f"RM only · Bray-Curtis PC1 vs PC2 (n={rm_mask.sum():,})")
axes[0].set_xlabel(axlabel("bc", 1)); axes[0].set_ylabel(axlabel("bc", 2))
axes[1].set_title(f"RM only · Weighted UniFrac PC1 vs PC2")
axes[1].set_xlabel(axlabel("wunifrac", 1)); axes[1].set_ylabel(axlabel("wunifrac", 2))
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", ncol=7, frameon=False,
           bbox_to_anchor=(0.5, 1.04), markerscale=4)
fig.tight_layout(rect=[0, 0, 1, 0.90])

png = OUT_DIR / "pcoa_rm_sites.png"
fig.savefig(png, dpi=150, bbox_inches="tight")
print(f"已写出 {png} ({png.stat().st_size/1024:.0f} KB)")
plt.close(fig)

print("\n全部图已生成。")
