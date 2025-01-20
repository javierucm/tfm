############################################################
# 0. Importante -> Necesario tener cargado previamente la Fig 10.  (teSig_ids)
############################################################

library(DESeq2)
library(pheatmap)


# teSig_ids <- scan("/ruta/a/TE_significativos_del_estudio_1mes.txt", what = character())

############################################################
# 1. LEER E IDENTIFICAR MUESTRAS NO TRATADAS (10 meses)
############################################################

# Lee la tabla de recuentos para el estudio de 10 meses NO TRATADAS
data_noTreated_10 <- read.table(
  "/home/vant/TFM/Browder2/TeTranscript/10month/TFM_JLR_Long10MonthYoungvsOldNoTreated.cntTable",
  header = TRUE,
  row.names = 1
)

# Renombrar columnas para evitar duplicados
colnames(data_noTreated_10) <- paste0("NoT_10m_", colnames(data_noTreated_10))

# Definir grupos para DESeq2 (TGroup vs CGroup)
groups_noTreated_10 <- factor(c(rep("TGroup", 3), rep("CGroup", 4)))

# Crear objeto colData
sampleInfo_noTreated_10 <- data.frame(
  groups = groups_noTreated_10,
  row.names = colnames(data_noTreated_10)
)

# Crear y correr DESeqDataSet
dds_noTreated_10 <- DESeqDataSetFromMatrix(
  countData = data_noTreated_10,
  colData   = sampleInfo_noTreated_10,
  design    = ~ groups
)
dds_noTreated_10 <- DESeq(dds_noTreated_10)

# Obtener recuentos normalizados
normCounts_noTreated_10 <- counts(dds_noTreated_10, normalized = TRUE)

############################################################
# 2. LEER E IDENTIFICAR MUESTRAS TRATADAS (10 meses)
############################################################

# Lee la tabla de recuentos para el estudio de 10 meses TRATADAS
data_treated_10 <- read.table(
  "/home/vant/TFM/Browder2/TeTranscript/10month/TFM_JLR_Long10MonthYoungvsOldTreated.cntTable",
  header = TRUE,
  row.names = 1
)
# Renombrar columnas
colnames(data_treated_10) <- paste0("Tr_10m_", colnames(data_treated_10))

# Definir grupos para DESeq2 (TGroup vs CGroup)
groups_treated_10 <- factor(c(rep("TGroup", 3), rep("CGroup", 4)))

sampleInfo_treated_10 <- data.frame(
  groups = groups_treated_10,
  row.names = colnames(data_treated_10)
)

dds_treated_10 <- DESeqDataSetFromMatrix(
  countData = data_treated_10,
  colData   = sampleInfo_treated_10,
  design    = ~ groups
)
dds_treated_10 <- DESeq(dds_treated_10)

# Obtener recuentos normalizados
normCounts_treated_10 <- counts(dds_treated_10, normalized = TRUE)

############################################################
# 3. COMBINAR Y GENERAR LA MATRIZ PARA EL HEATMAP
############################################################

# 3.1. Intersección de features (filas)
commonFeatures_10 <- intersect(rownames(normCounts_noTreated_10), 
                               rownames(normCounts_treated_10))

# Unir los recuentos normalizados en una sola matriz
all_normCounts_10 <- cbind(
  normCounts_noTreated_10[commonFeatures_10, ],
  normCounts_treated_10[commonFeatures_10, ]
)

# 3.2. Usar los TEs del estudio anterior (teSig_ids)
#     --> Nos quedamos con aquellos que también existan en la intersección
teSig_ids_10 <- intersect(teSig_ids, commonFeatures_10)

cat("Total TE del primer estudio presentes en el de 10 meses:", length(teSig_ids_10), "\n")
if (length(teSig_ids_10) == 0) {
  stop("No hay TE significativos en la intersección para el estudio de 10 meses.")
}

mat_heatmap_10 <- all_normCounts_10[teSig_ids_10, ]

############################################################
# 4. ANOTACIÓN DE COLUMNAS (Tipos de muestra)
############################################################

# Aquí puedes ajustar la anotación a lo que quieras reflejar para tu estudio
# Por ejemplo, usando la misma idea que antes:
# - TGroup (noTreated_10) -> "25 months -doxycycline (10m)"
# - CGroup (noTreated_10) -> "3 months -doxycycline (10m)"
# - TGroup (treated_10)   -> "25 months +doxycycline (10m)"
# - CGroup (treated_10)   -> "3 months -doxycycline (10m)"

sampleTypes_noTreated_10 <- c(
  rep("22 months -doxycycline (10m)", 3),
  rep("3 months -doxycycline ", 4)
)
sampleTypes_treated_10 <- c(
  rep("22 months +doxycycline (10m)", 3),
  rep("3 months -doxycycline ", 4)
)

# Unificamos
anno_sampleTypes_10 <- c(sampleTypes_noTreated_10, sampleTypes_treated_10)
annotation_col_10 <- data.frame(SampleType = anno_sampleTypes_10)
rownames(annotation_col_10) <- c(colnames(data_noTreated_10), colnames(data_treated_10))

############################################################
# 5. (Opcional) Eliminar las últimas 4 columnas
#    Tal como se hizo en el código original
############################################################

# Asegúrate de que tienes al menos 14 columnas totales como en el ejemplo
if (ncol(mat_heatmap_10) >= 14) {
  mat_heatmap_10 <- mat_heatmap_10[, -c(11:14)]
  annotation_col_10 <- annotation_col_10[-c(11:14), , drop = FALSE]
} else {
  message("Advertencia: No hay suficientes columnas para eliminar las posiciones 11 a 14.")
}

############################################################
# 6. HEATMAP con distancia de Pearson en filas y columnas
############################################################

# Distancia de Pearson en columnas
correlation_matrix_cols_10 <- cor(mat_heatmap_10, method = "pearson")
distance_matrix_cols_10 <- as.dist(1 - correlation_matrix_cols_10)

# Distancia de Pearson en filas (correlación entre TEs, se calcula en la transpuesta)
correlation_matrix_rows_10 <- cor(t(mat_heatmap_10), method = "pearson")
distance_matrix_rows_10 <- as.dist(1 - correlation_matrix_rows_10)

# Generar PDF
pdf("/home/vant/TFM/Browder2/Figuras/Fig 10b: Heatmap_TE_sig_diff_10m.pdf", width = 10, height = 10)

pheatmap(
  mat_heatmap_10,
  scale                   = "row",
  annotation_col         = annotation_col_10,
  show_rownames          = TRUE,
  show_colnames          = TRUE,
  main = "Fig 10b: Heatmap TE (10 month) using TEs from 1 month treatment",
  clustering_distance_cols = distance_matrix_cols_10,
  clustering_distance_rows = distance_matrix_rows_10,
  clustering_method       = "complete",
  fontsize                = 8
)

dev.off()
