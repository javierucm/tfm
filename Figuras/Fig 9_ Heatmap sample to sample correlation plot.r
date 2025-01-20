#####################################
# Obtener el directorio de trabajo actual

current_dir='/home/vant/TFM/Browder2'

# Leer los dos archivos de conteo de TE
#data1 <- read.table(file.path(current_dir, "TFM_JLR_Long10MonthYoungvsOldNoTreated.cntTable"), header=T, row.names=1)
#data2 <- read.table(file.path(current_dir, "TFM_JLR_Long10MonthYoungvsOldTreated.cntTable"), header=T, row.names=1)
data1 <- read.table(file.path(current_dir, "TeTranscript/1month/TFM_JLR_YoungvsOld_NoTreated.cntTable"), header=T, row.names=1)
data2 <- read.table(file.path(current_dir, "TeTranscript/1month/TFM_JLR_YoungvsOldTreated.cntTable"), header=T, row.names=1)



data1 <- data1[!grepl("^ENS", rownames(data1)), ]
data2 <- data2[!grepl("^ENS", rownames(data2)), ]


# Asegurarse de que los genes (filas) sean los mismos
common_genes <- intersect(rownames(data1), rownames(data2))
data1 <- data1[common_genes, ]
data2 <- data2[common_genes, ]

# Verificar si los datos de los controles son idénticos
if(!identical(data1[,4:7], data2[,4:7])) {
    stop("Los datos de los controles no son idénticos entre los conjuntos de datos.")
}


# Combinar los datos (usando los controles de un solo archivo)
data_combined <- cbind(data1[,4:7], data1[,1:3], data2[,1:3])

# Definir los grupos
groups <- factor(c(rep("Control", 4), rep("Old_NoTreated", 3), rep("Old_Treated", 3)))

# Crear la información de muestra
sampleInfo <- data.frame(groups, row.names=colnames(data_combined))

# Asignar los datos combinados
data <- data_combined

# Filtrar genes con baja expresión
min_read <- 1
data <- data[apply(data, 1, function(x){max(x)}) > min_read, ]

# Cargar DESeq2 y crear el objeto DESeqDataSet
suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSetFromMatrix(countData = data, colData = sampleInfo, design = ~ groups)

# Reajustar el nivel de referencia del grupo
dds$groups <- relevel(dds$groups, ref="Control")

# Ejecutar el análisis de DESeq2
dds <- DESeq(dds)

# Extraer los resultados para las comparaciones de interés
res_NoTreated_vs_Control <- results(dds, contrast = c("groups", "Old_NoTreated", "Control"))
res_Treated_vs_Control <- results(dds, contrast = c("groups", "Old_Treated", "Control"))

#write.table(res_NoTreated_vs_Control, file="a2.txt",sep="\t", quote=F)
#write.table(res_Treated_vs_Control, file="b2.txt",sep="\t", quote=F)




# Filtrar y combinar los genes significativos
#resSig_NoTreated <- res_NoTreated_vs_Control[(!is.na(res_NoTreated_vs_Control$padj) & (res_NoTreated_vs_Control$padj < 0.05)), ]
#resSig_Treated <- res_Treated_vs_Control[(!is.na(res_Treated_vs_Control$padj) & (res_Treated_vs_Control$padj < 0.05)), ]


genes_sig <- unique(c(rownames(res_NoTreated_vs_Control), rownames(res_Treated_vs_Control)))
#genes_sig <- intersect(rownames(resSig_NoTreated), rownames(resSig_Treated))




# Extraer y escalar los conteos normalizados
normalized_counts <- counts(dds, normalized=TRUE)
norm_counts_sig <- normalized_counts[genes_sig, ]
norm_counts_scaled <- t(scale(t(norm_counts_sig)))




library(pheatmap)

# Suponiendo que tienes 'normalized_counts' (genes en filas, muestras en columnas)
# 1. Calculamos la correlación (Pearson o Spearman) entre columnas (muestras)
cor_matrix <- cor(normalized_counts, method = "pearson") 
# O 'spearman' si prefieres

# 2. Creamos un heatmap de la matriz de correlación
pheatmap(cor_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Fig 9:Sample-to-Sample Correlations TE (1month treatment)",
         annotation_col = sampleInfo,  # si quieres anotar el grupo
         fontsize = 10)

annotation_col <- data.frame(Group = sampleInfo$groups)
rownames(annotation_col) <- colnames(normalized_counts)
annotation_colors <- list(Group = c("Control"="#FFCC00", "Old_NoTreated"="#FF0000", "Old_Treated"="#0000FF"))
# Abre un dispositivo JPG:
png(file.path(current_dir,"Figuras/Fig 9: Heatmap Sample-to-Sample Correlations TE.jpg"), width = 800, height = 600)

# Llamas a pheatmap normalmente
pheatmap(
  cor_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Fig 9:Sample-to-Sample Correlations TE (1month treatment)",
  annotation_col = annotation_col,
  fontsize = 10
)

# Cierras el dispositivo
dev.off()