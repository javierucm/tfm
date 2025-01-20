#####################################
# Obtener el directorio de trabajo actual

current_dir='/home/vant/TFM/Browder2'

# Leer los dos archivos de conteo de TE
#data1 <- read.table(file.path(current_dir, "TFM_JLR_Long10MonthYoungvsOldNoTreated.cntTable"), header=T, row.names=1)
#data2 <- read.table(file.path(current_dir, "TFM_JLR_Long10MonthYoungvsOldTreated.cntTable"), header=T, row.names=1)
data1 <- read.table(file.path(current_dir, "TeTranscript/1month/TFM_JLR_YoungvsOld_NoTreated.cntTable"), header=T, row.names=1)
data2 <- read.table(file.path(current_dir, "TeTranscript/1month/TFM_JLR_YoungvsOldTreated.cntTable"), header=T, row.names=1)



data1 <- data1[grepl("^ENS", rownames(data1)), ]
data2 <- data2[grepl("^ENS", rownames(data2)), ]


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


# Usar los conteos normalizados de los genes significativos
pca_data <- t(norm_counts_scaled)  # Transponer: muestras como filas, genes como columnas

# Realizar el PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# Extraer las coordenadas del PCA para las muestras
pca_coordinates <- as.data.frame(pca_result$x)

# Añadir la información de los grupos
pca_coordinates$Group <- sampleInfo$groups

# Visualizar la varianza explicada por cada componente
explained_variance <- summary(pca_result)$importance[2, ]


# Invertir los ejes si es necesario (para que se asemeje al del Article Research)
#pca_coordinates$PC1 <- -pca_coordinates$PC1  # Invertir PC1
pca_coordinates$PC2 <- -pca_coordinates$PC2  # Invertir PC2 

# Generar el gráfico de PCA
library(ggplot2)
png(file.path(current_dir,"Figuras/Fig 1: PCA of Gene Normalized Counts (1 month treatment).png"), width = 800, height = 600)
ggplot(pca_coordinates, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  labs(title = "Fig 1: PCA of Gene Normalized Counts (1 month treatment)",
       x = paste0("PC1 (", round(explained_variance[1] * 100, 2), "%)"),
       y = paste0("PC2 (", round(explained_variance[2] * 100, 2), "%)")) +
  theme_minimal() +
  scale_color_manual(values = c("Control" = "orange", "Old_NoTreated" = "red", "Old_Treated" = "blue"),
   labels = c("3 months -doxycycline", "26 months -doxycycline", "26 months +doxycycline (1 month)")
  ) +
  theme(
    legend.position = "top",
    plot.title = element_text(size = 20, face = "bold"), # Título del gráfico
    axis.title = element_text(size = 16, face = "bold"), # Títulos de los ejes
    axis.text = element_text(size = 16),                # Texto de los ejes
    legend.text = element_text(size = 16),              # Texto de la leyenda
    legend.title = element_text(size = 16, face = "bold") # Título de la leyenda
  )

dev.off()

# Mostrar los genes que más contribuyen a PC1 y PC2
print("Top genes contributing to PC1:")
print(sort(abs(pca_result$rotation[, "PC1"]), decreasing = TRUE)[1:10])

print("Top genes contributing to PC2:")
print(sort(abs(pca_result$rotation[, "PC2"]), decreasing = TRUE)[1:10])