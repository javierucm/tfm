library(reshape)
library(ggplot2)

current_dir <- '/home/vant/TFM/Browder2'

# Cargar los datos
data_no_treated <- read.table(file.path(current_dir, "TeTranscript/10month/TFM_JLR_Long10MonthYoungvsOldNoTreated_gene_TE_analysis.txt"), header = TRUE, sep = "\t", row.names = 1)
data_treated <- read.table(file.path(current_dir, "TeTranscript/10month/TFM_JLR_Long10MonthYoungvsOldTreated_gene_TE_analysis.txt"), header = TRUE, sep = "\t", row.names = 1)

# Filtrar TEs (no comienzan con "ENS")
TE_no_treated <- data_no_treated[!grepl("^ENS", rownames(data_no_treated)), ]
TE_treated <- data_treated[!grepl("^ENS", rownames(data_treated)), ]

# Crear un conjunto combinado de TEs (UNIÓN)
all_TEs <- unique(c(rownames(TE_no_treated), rownames(TE_treated)))

# Crear una matriz combinada de log2FoldChange
combined_data <- data.frame(
  No_Treatment = ifelse(all_TEs %in% rownames(TE_no_treated), TE_no_treated$log2FoldChange[match(all_TEs, rownames(TE_no_treated))], NA),
  Treatment = ifelse(all_TEs %in% rownames(TE_treated), TE_treated$log2FoldChange[match(all_TEs, rownames(TE_treated))], NA),
  row.names = all_TEs
)

# Asegurarse de que las filas tengan nombres únicos (evitar duplicados)
combined_data <- combined_data[all_TEs, ]

combined_data_no_na <- combined_data
combined_data_no_na[is.na(combined_data_no_na)] <- 0

# Convertir los datos a formato largo
long_data <- melt(as.matrix(combined_data), varnames = c("TE", "Condition"), value.name = "log2FoldChange")
# Crear el boxplot
png(file.path(current_dir, "Figuras/Fig 11b: Boxplot 10 month treatment.png"), width = 800, height = 600)
ggplot(long_data, aes(x = Condition, y = value, fill = Condition)) +
  geom_boxplot() +
  labs(
    title = "Fig 11b: Distribución de log2FoldChange por condición (10 month treatment)",
    x = "Condición",
    y = "log2FoldChange",
    fill = "Condición"
  ) +
  theme_minimal() +
  scale_fill_manual(
    values = c("blue", "red"),
    labels = c("22 months -doxycycline / 3 months -doxycycline", 
               "22 months +doxycycline / 3 months -doxycycline")
  ) +
  theme(
    plot.title = element_text(size = 20, face = "bold"), # Título más grande
    axis.title = element_text(size = 16, face = "bold"), # Etiquetas de los ejes más grandes
    axis.text = element_text(size = 14),                # Texto de los ejes más grande
    legend.title = element_text(size = 16, face = "bold"), # Título de la leyenda más grande
    legend.text = element_text(size = 14)               # Texto de la leyenda más grande
  )
dev.off()


# Crear un MA plot para la condición No Tratada
png(file.path(current_dir,"Figuras/Fig 12b: MA Plot no treated Condition (22 month old).png"), width = 800, height = 600)
plot(TE_no_treated$baseMean,
     TE_no_treated$log2FoldChange,
     log = "x",
     pch = 19,
     col = ifelse(TE_no_treated$padj < 0.05, "red", "black"),
     xlab = "Mean of normalized counts (log scale)",
     ylab = "log2 Fold Change",
     main = "Fig 12b: MA Plot - No Treatment Condition (22 month old)",
     cex.lab = 1.3,    # Tamaño de las etiquetas de los ejes
     cex.main = 1.5,     # Tamaño del título
     cex.axis = 1.3    # Tamaño de los números en los ejes
)
abline(h = c(-1, 1), col = "blue", lty = 2) # Líneas de referencia
dev.off()


# Crear un MA plot para la condición Tratada
png(file.path(current_dir,"Figuras/Fig 13b: MA Plot treated Condition (10 month treatment - 22 month old).png"), width = 800, height = 600)
plot(TE_treated$baseMean,
     TE_treated$log2FoldChange,
     log = "x",
     pch = 19,
     col = ifelse(TE_treated$padj < 0.05, "red", "black"),
     xlab = "Mean of normalized counts (log scale)",
     ylab = "log2 Fold Change",
     main = "Fig 13b: MA Plot - Treatment Condition (10 month treatment - 22 month old)",
     cex.lab = 1.3,    # Tamaño de las etiquetas de los ejes
     cex.main = 1.5,     # Tamaño del título
     cex.axis = 1.3    # Tamaño de los números en los ejes
)
abline(h = c(-1, 1), col = "blue", lty = 2) # Líneas de referencia
dev.off()
