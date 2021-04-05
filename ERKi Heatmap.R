library(ComplexHeatmap)
library(openxlsx)
setwd("C:/Lab/R studio")
ktr<-read.xlsx('Combined H1299 bins.xlsx')
mat= as.matrix(ktr[-c(1:2),-1])
class(mat)<-"numeric"
rownames(mat) <- ktr[-c(1:2),1]
col_heat = colorRamp2(c(0, 7.5, 15), c("white","#9bdb24","#4e6e12"))
ha = HeatmapAnnotation(
  Drug= as.character(ktr[1,-1]),
  Hours= as.character(ktr[2,-1]),
  col = list(Drug = c("DMSO" = "#adadad", "LY" = "#90bff9","PD" = "#323232","Rav" = "#800080","SCH" = "#f71480","Ulix" = "#c0c000","VX" = "#00c000"),
             Hours = c('24' = '#d3e5fa', '48' = '#96a7bb', '72' = '#5b6c81', '96' = '#2c323a')),
  gp = gpar(col = "black")
  )
Heatmap(mat, col=col_heat, name = "Percentage", 
        clustering_distance_rows = 'euclidean', clustering_method_columns= 'ward.D', cluster_rows = FALSE, 
        top_annotation = ha, row_names_side = "left", show_column_names = FALSE, column_dend_height = unit(1, "cm"),
        row_title = "ERK activity bins")
