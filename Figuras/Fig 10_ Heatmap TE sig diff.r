# ----------------------------
# 1. LEER E IDENTIFICAR MUESTRAS NO TRATADAS
# ----------------------------
data_noTreated <- read.table(
  "/home/vant/TFM/Browder2/TeTranscript/1month/TFM_JLR_YoungvsOld_NoTreated.cntTable",
  header = TRUE,
  row.names = 1
)

# Renombrar columnas para evitar duplicados
colnames(data_noTreated) <- paste0("NoT_", colnames(data_noTreated))

# Definir grupos para DESeq2 (TGroup vs CGroup)
groups_noTreated <- factor(c(rep("TGroup", 3), rep("CGroup", 4)))

# Crear objeto colData
sampleInfo_noTreated <- data.frame(
  groups = groups_noTreated,
  row.names = colnames(data_noTreated)
)

library(DESeq2)
dds_noTreated <- DESeqDataSetFromMatrix(
  countData = data_noTreated,
  colData   = sampleInfo_noTreated,
  design    = ~ groups
)
dds_noTreated <- DESeq(dds_noTreated)
res_noTreated <- results(dds_noTreated)

# Filtramos TEs significativamente regulados con p-adj < 0.05 y |log2FC| > 0.58
resSig_noTreated <- res_noTreated[
  (!is.na(res_noTreated$padj)) & (res_noTreated$padj < 0.05) &
  (abs(res_noTreated$log2FoldChange) > 0.58),
]

# Excluir los "genes" que empiecen por ENS (quedarnos con TEs)
teSig_noTreated <- resSig_noTreated[!grepl("^ENS", rownames(resSig_noTreated)), ]

# Obtener recuentos normalizados
normCounts_noTreated <- counts(dds_noTreated, normalized = TRUE)


# ----------------------------
# 2. LEER E IDENTIFICAR MUESTRAS TRATADAS
# ----------------------------
data_treated <- read.table(
  "/home/vant/TFM/Browder2/TeTranscript/1month/TFM_JLR_YoungvsOldTreated.cntTable",
  header = TRUE,
  row.names = 1
)
# Renombrar columnas
colnames(data_treated) <- paste0("Tr_", colnames(data_treated))

# Definir grupos para DESeq2 (TGroup vs CGroup)
groups_treated <- factor(c(rep("TGroup", 3), rep("CGroup", 4)))

sampleInfo_treated <- data.frame(
  groups = groups_treated,
  row.names = colnames(data_treated)
)

dds_treated <- DESeqDataSetFromMatrix(
  countData = data_treated,
  colData   = sampleInfo_treated,
  design    = ~ groups
)
dds_treated <- DESeq(dds_treated)
normCounts_treated <- counts(dds_treated, normalized = TRUE)


# ----------------------------
# 3. COMBINAR Y GENERAR HEATMAP
# ----------------------------

# 3.1. Intersección de features
commonFeatures <- intersect(rownames(normCounts_noTreated), rownames(normCounts_treated))
all_normCounts <- cbind(
  normCounts_noTreated[commonFeatures, ],
  normCounts_treated[commonFeatures, ]
)

# 3.2. Subset con los TE significativos en 'noTreated' (y que existan en treated)
teSig_ids <- rownames(teSig_noTreated)
teSig_ids <- intersect(teSig_ids, commonFeatures) 

cat("Total TE sig en intersección:", length(teSig_ids), "\n")
if (length(teSig_ids) == 0) {
  stop("No hay TE significativos en la intersección para generar heatmap.")
}

mat_heatmap <- all_normCounts[teSig_ids, ]


# 3.3. Anotación unificada (3 tipos de muestra)
# - Para las muestras "noTreated": TGroup -> "25m_NoTreated", CGroup -> "Control"
# - Para las muestras "treated":   TGroup -> "25m_Treated",   CGroup -> "Control"
sampleTypes_noTreated <- c(rep("25 months -doxycycline", 3), rep("3 months -doxycycline", 4))
sampleTypes_treated   <- c(rep("25 months +doxycycline", 3), rep("3 months -doxycycline", 4))

# Unificamos
anno_sampleTypes <- c(sampleTypes_noTreated, sampleTypes_treated)
annotation_col <- data.frame(SampleType = anno_sampleTypes)
rownames(annotation_col) <- c(colnames(data_noTreated), colnames(data_treated))


# 3.4. Eliminar las 4 últimas columnas (y sus anotaciones) tal como solicitas
mat_heatmap <- mat_heatmap[, -c(11:14)]
annotation_col <- annotation_col[-c(11:14), , drop = FALSE]


# ----------------------------
# 4. HEATMAP con distancia de Pearson en filas y columnas
# ----------------------------
library(pheatmap)

# Distancia de Pearson en columnas
correlation_matrix_cols <- cor(mat_heatmap, method = "pearson")
distance_matrix_cols <- as.dist(1 - correlation_matrix_cols)

# Distancia de Pearson en filas (correlación entre TEs)
# OJO: la correlación en filas se obtiene transponiendo mat_heatmap
correlation_matrix_rows <- cor(t(mat_heatmap), method = "pearson")
distance_matrix_rows <- as.dist(1 - correlation_matrix_rows)

pdf("/home/vant/TFM/Browder2/Figuras/Fig 10 : Heatmap TE sig diff.pdf", width = 10, height = 10)
pheatmap(
  mat_heatmap,
  scale                   = "row",
  annotation_col         = annotation_col,
  show_rownames          = TRUE,
  show_colnames          = TRUE,
  main = "Fig 10: Heatmap TE sig diff (25 months -dox/3 months -dox)(1month treatment)",
  clustering_distance_cols = distance_matrix_cols,
  clustering_distance_rows = distance_matrix_rows,
  clustering_method       = "complete",
  fontsize                = 8
)
dev.off()
