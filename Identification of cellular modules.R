meta <- data.table::fread("./metadata_total_8.14.csv")
meta <- data.frame(meta[,-1], row.names = meta$V1)
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(CoVarNet)
library(purrr)

get_prop_mat <- function(pattern, celltype_col = "celltype2") {
  meta %>%
    filter(str_detect(.data[[celltype_col]], pattern)) %>%
    count(sample, .data[[celltype_col]], name = "Freq") %>%
    pivot_wider(names_from = all_of(celltype_col), values_from = Freq, values_fill = 0) %>%
    column_to_rownames("sample") %>%
    {./rowSums(.)}
}

task_list <- list(
  TNK    = "CD4|CD8|NK|ILC",
  Myeloid= "Mo|Macro|RTM|DC",
  Bcell  = "B_",
  Plasma = "Plasma",
  Neu    = "Neu",
  Mast   = "Mast",
  Mes    = "Fib|SMC|Peri",
  EC     = "EC"
)

lst <- map(task_list, ~get_prop_mat(.x))
ref_samples <- rownames(lst[[1]])
lst_align <- map(lst, ~ {
  tmp <- .x[ref_samples, , drop = FALSE]
  tmp[is.na(tmp)] <- 0
  tmp
})

dat_bind <- bind_cols(lst_align)
dat_bind[is.na(dat_bind)] <- 0

mat_fq_raw <- t(dat_bind) %>% as.data.frame()

cluster_map <- c(
  "CD4|CD8|NK|ILC" = "TNK",
  "Mast"           = "Mast",
  "Neu"            = "Neu",
  "Macro|Mo|DC|RTM"= "Myeloid",
  "B_|Plasma"      = "B",
  "Fib|SMC|Peri"   = "Mes",
  "EC"             = "EC"
)

subtypes <- rownames(mat_fq_raw)
mat_fq_raw$majorCluster <- "NA"
for (pat in names(cluster_map)) {
  mat_fq_raw$majorCluster[str_detect(subtypes, pat)] <- cluster_map[pat]
}

mat_fq_raw <- mat_fq_raw %>%
  select(majorCluster, everything())
mat_fq_norm <- freq_normalize(mat_fq_raw, normalize = "minmax")

col_fun <- colorRamp2(seq(-1, 1, length.out = 11), brewer.pal(11, "RdBu")[11:1])
heat <- Heatmap(
  cor(t(mat_fq_norm)),
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  clustering_method_rows = "ward.D2",
  clustering_distance_rows = "pearson",
  clustering_distance_columns = "pearson",
  clustering_method_columns = "ward.D2",
  row_split = 11,
  column_split = 11,
  row_gap = unit(0, "mm"),
  column_gap = unit(0, "mm")
)
heat

row_split_order <- row_order(heat)
all_row_names <- rownames(mat_fq_norm)
cell_module_list <- lapply(row_split_order, function(idx) {
  all_row_names[idx]
})
names(cell_module_list) <- paste0("CM_", seq_along(cell_module_list))
cell_module_list

res <- nmf(mat_fq_norm, rank = 2:20, method = "nsNMF", seed = rep(123456, 6), .options = "vp")
plot(res)
K=13
NMF_K13<-nmf(mat_fq_norm, K, method = "nsNMF", seed = rep(77, 6), nrun = 30)
colnames(basis(NMF_K13))=c(paste0("CM_0",1:9),paste0("CM_",10:13))
rownames(coef(NMF_K13))=c(paste0("CM_0",1:9),paste0("CM_",10:13))
h<-coef(NMF_K13)
max_cm <- apply(h, 2, function(x) rownames(h)[which.max(x)])
max_cm<-gsub("CM", "CMT",max_cm)
max_cm[1:10]
top_features_list <- extractFeatures(NMF_K13, 0.3)
top_features_list <- lapply(top_features_list, function(x) rownames(mat_fq_norm)[x])
distance_matrix <- outer(top_features_list, cell_module_list, FUN = Vectorize(function(x, y) {
  intersection <- length(intersect(x, y))
  union <- length(union(x, y))
  if (union == 0) return(0)
  intersection / union
}))
rownames(distance_matrix) <- paste0("NMF",1:13)
colnames(distance_matrix) <- paste0("CM_",1:11)
pdf("CM_similarity.pdf",width = 7, height = 6)
Heatmap(distance_matrix
        ,col = colorRamp2(seq(0,0.4,length = 6),brewer.pal(11,"RdBu")[6:1]),
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = T,
        show_row_names = T,
        clustering_method_rows = "ward.D2",
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "pearson",
        clustering_method_columns = "ward.D2",
        column_gap = unit(0, "mm"),
        row_gap = unit(0, "mm"),
        name = "Similarity",
)
dev.off()
