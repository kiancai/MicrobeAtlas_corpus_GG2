# SampleDistance100k

Independent 100k sample-distance branch for corpus PCoA figures. This does not
modify the main corpus files or the existing 50k branch.

## Outputs

```text
results/sample_distance_100k/
├── subset_100k_index.tsv
├── subset_100k.h5ad
├── genus_tree.nwk
├── distance_bc.float16.npy
├── distance_wunifrac.float16.npy
├── pcoa_coords_bc.npy
├── pcoa_coords_wunifrac.npy
├── pcoa_eigenvalues.tsv
├── pcoa_negative_eigen_audit.tsv
└── figures/
```

## Sampling design

Total sample count is 100,000:

- 10,000 cross-database paired Runs, included as both MA and RM rows
- 70,000 additional MicrobeAtlas rows
- 10,000 additional ResMicroDb rows

Final database balance is 80,000 MicrobeAtlas and 20,000 ResMicroDb rows.
Additional MA quota is:

| group | quota |
|---|---:|
| Human | 22,000 |
| Animal_other | 12,000 |
| Soil | 13,000 |
| Aquatic | 12,000 |
| Plant | 6,000 |
| Unknown | 5,000 |

Paired Runs and additional RM rows are stratified by `RM_Sample_Site` and then
`Project_ID` with sqrt allocation. MA groups are stratified by body site or
environment subcategory with sqrt allocation.

## Run

Submit from the project root. This cluster wrapper requires explicit resources
on the `sbatch` command line even though the scripts include `#SBATCH` headers.

```bash
cd /hpcdisk1/limk_group/caiqy/project/260428_greengene2

sbatch -c 4  --mem=64G  -t 01:00:00 scripts/SampleDistance100k/01_stratified_sample_100k.sbatch
sbatch -c 8  --mem=128G -t 02:00:00 scripts/SampleDistance100k/02_build_subset_anndata_100k.sbatch
sbatch -c 64 --mem=450G -t 06:00:00 scripts/SampleDistance100k/03_compute_bc_100k.sbatch
sbatch -c 64 --mem=450G -t 04:00:00 scripts/SampleDistance100k/04_compute_wunifrac_100k.sbatch
sbatch -c 64 --mem=450G -t 06:00:00 scripts/SampleDistance100k/05_pcoa_100k.sbatch
sbatch -c 4  --mem=64G  -t 01:00:00 scripts/SampleDistance100k/06_plot_pcoa_100k.sbatch
```

## Figures

The main figure is `figures/pcoa_100k_4x2.png`.

Rows:

1. MA environment overview: Animal / Soil / Aquatic / Plant / Unknown
2. MA human body sites
3. ResMicroDb in corpus space
4. Cross-database paired Runs

Columns:

1. Bray-Curtis
2. Weighted normalized UniFrac

Every panel draws all 100,000 samples as a gray background and only changes the
highlight layer. Coordinates are therefore directly comparable within each
distance metric.
