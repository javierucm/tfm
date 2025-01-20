##########################################################################3
     ############Version  Volcano Plot Paper:

     library(ggplot2)

# Cargar datos
current_dir='/home/vant/TFM/Browder2'

# Leer los dos archivos de conteo de TE
data <- read.table(file.path(current_dir, "/TeTranscript/1month/TFM_JLR_YoungvsOldTreated_gene_TE_analysis.txt"), header=T, row.names=1)


# Filtrar datos relevantes
data <- data[!grepl("^ENS", rownames(data)), ]

# Clasificar genes
data$category <- "No significativos"
data$category[data$padj < 0.05 & data$log2FoldChange >= 0.585] <- "Upregulados"
data$category[data$padj < 0.05 & data$log2FoldChange <= -0.585] <- "Downregulados"
data$category[data$padj < 0.05 & abs(data$log2FoldChange) < 0.585] <- "Cambios no marcados"

# Crear el volcano plot ajustado
png(file.path(current_dir,"Figuras/Fig 8: Volcanoplot_All_TE_Paper_YoungVsOldTreated_1month.png"), width = 800, height = 600)
ggplot(data, aes(x = log2FoldChange, y = -log10(padj), color = category)) +
  geom_point(alpha = 0.9, size = 1) +
  scale_color_manual(values = c("No significativos" = "grey", 
                                "Upregulados" = "red", 
                                "Downregulados" = "blue", 
                                "Cambios no marcados" = "black")) +
  # Añadir líneas verticales en LFC = ±1
  geom_vline(xintercept = c(-1, 1), linetype = "solid", color = "black") +
  # Añadir línea horizontal en -log10(padj) = 1.301 (equivalente a padj=0.05)
  geom_hline(yintercept = -log10(0.05), linetype = "solid", color = "black") +
  labs(title = "Fig 8: TE: 26 months +doxycycline(1month) / 3 months -doxycycline",
       x = "Log2 Fold Change",
       y = "-log10(p-value adjusted)") +
  xlim(c(-10, 10)) +
  ylim(c(0, 8)) +
  annotate("text", x = 5, y = 7.5, label = paste(nrow(data[data$category == "Upregulados", ]), "Upregulados"), 
           color = "red", size = 6, fontface = "bold") + # Ajusta tamaño y estilo
  annotate("text", x = -5, y = 7.5, label = paste(nrow(data[data$category == "Downregulados", ]), "Downregulados"), 
           color = "blue", size = 6, fontface = "bold") + # Ajusta tamaño y estilo
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"), # Título del gráfico
    axis.title = element_text(size = 16, face = "bold"), # Títulos de los ejes
    axis.text = element_text(size = 16),                # Texto de los ejes
    legend.text = element_text(size = 16),              # Texto de la leyenda (si fuera visible)
    legend.title = element_text(size = 16, face = "bold"), # Título de la leyenda (si fuera visible)
    legend.position = "none"                             # Mantener sin leyenda como en el código original
  )

dev.off()

