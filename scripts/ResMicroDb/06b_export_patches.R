# 06b_export_patches.R
# 从 ps.16s.rds 和 ps.16s_0105_new7.rds 中导出 Run 级 patch 映射，
# 供 06b_fix_metadata_errors.py 消费。
#
# 输出：results/feature_table/metadata_patches_0105.tsv
#   只包含 5 个目标 project（PRJNA801796/822681/824137/914884/1058141）里
#   OLD vs NEW (0105) 在 {Body_Site, Phenotype, Case_Or_Control, Is_Healthy}
#   任一列有差异的 Run。预期 1,641 行。
#
# 改动规则与 jxt 漏改的处理详见：
# rawdata/ResMicroDb/supplement data/CHANGES_0105_new7.md
#
# 运行：conda activate baseR && Rscript scripts/ResMicroDb/06b_export_patches.R

suppressPackageStartupMessages({ library(phyloseq) })

ROOT  <- "/hpcdisk1/limk_group/caiqy/project/260428_greengene2"
F_OLD <- file.path(ROOT, "rawdata/ResMicroDb/supplement data/ps.16s.rds")
F_NEW <- file.path(ROOT, "rawdata/ResMicroDb/supplement data/ps.16s_0105_new7.rds")
OUT   <- file.path(ROOT, "results/feature_table/metadata_patches_0105.tsv")

PATCH_PROJECTS <- c("PRJNA801796","PRJNA822681","PRJNA824137",
                    "PRJNA914884","PRJNA1058141")

ps_old <- readRDS(F_OLD); ps_new <- readRDS(F_NEW)
sd_old <- as(sample_data(ps_old), "data.frame")
sd_new <- as(sample_data(ps_new), "data.frame")

shared <- intersect(rownames(sd_old), rownames(sd_new))
so <- sd_old[shared, , drop=FALSE]
sn <- sd_new[shared, , drop=FALSE]

mask_proj <- so$Project_ID %in% PATCH_PROJECTS
sub_runs  <- shared[mask_proj]
cat("5 个 patch project 在 OLD∩NEW 共 ", length(sub_runs), " 个样本\n", sep="")

eq_safe <- function(a, b) {
  a <- as.character(a); b <- as.character(b)
  a[is.na(a)] <- "<NA>"; b[is.na(b)] <- "<NA>"
  a == b
}

cols_to_check <- c("Body_Site","Phenotype","Case_Or_Control","Is_Healthy")
diff_mask <- rep(FALSE, length(sub_runs))
for (c in cols_to_check) {
  diff_mask <- diff_mask | !eq_safe(so[sub_runs, c], sn[sub_runs, c])
}
diff_runs <- sub_runs[diff_mask]
cat("OLD vs NEW 在 ", paste(cols_to_check, collapse=","),
    " 任一列有差异: ", length(diff_runs), " 行\n", sep="")

out <- data.frame(
  Run                 = diff_runs,
  Project_ID          = as.character(so[diff_runs, "Project_ID"]),
  Body_Site_old       = as.character(so[diff_runs, "Body_Site"]),
  Body_Site_new       = as.character(sn[diff_runs, "Body_Site"]),
  Phenotype_old       = as.character(so[diff_runs, "Phenotype"]),
  Phenotype_new       = as.character(sn[diff_runs, "Phenotype"]),
  Phenotype_ID_old    = as.character(so[diff_runs, "Phenotype_ID"]),
  Phenotype_ID_new    = as.character(sn[diff_runs, "Phenotype_ID"]),
  Case_Or_Control_old = as.character(so[diff_runs, "Case_Or_Control"]),
  Case_Or_Control_new = as.character(sn[diff_runs, "Case_Or_Control"]),
  Is_Healthy_old      = as.character(so[diff_runs, "Is_Healthy"]),
  Is_Healthy_new      = as.character(sn[diff_runs, "Is_Healthy"]),
  host_disease        = as.character(sn[diff_runs, "host.disease"]),
  real_group          = as.character(sn[diff_runs, "real_group"]),
  stringsAsFactors    = FALSE
)

# 各 patch 命中数
cat("\n各 project 命中行数:\n")
print(table(out$Project_ID))

cat("\nSample_Site (Body_Site_old→new):\n")
print(table(out$Body_Site_old, out$Body_Site_new))

cat("\nPhenotype (old→new) — PRJNA801796:\n")
ph_sub <- out[out$Project_ID == "PRJNA801796", c("Phenotype_old","Phenotype_new")]
print(table(ph_sub$Phenotype_old, ph_sub$Phenotype_new))

cat("\nPRJNA822681 host.disease 分布:\n")
print(table(out[out$Project_ID == "PRJNA822681", "host_disease"]))

cat("\nPRJNA824137 real_group 类型分布 (HC* vs TBZ*):\n")
rg <- out[out$Project_ID == "PRJNA824137", "real_group"]
print(table(ifelse(startsWith(rg, "HC"), "HC*",
            ifelse(startsWith(rg, "TBZ"), "TBZ*", "other"))))

dir.create(dirname(OUT), showWarnings=FALSE, recursive=TRUE)
write.table(out, OUT, sep="\t", row.names=FALSE, quote=FALSE, na="")
cat("\n写出: ", OUT, " (", nrow(out), " 行)\n", sep="")
