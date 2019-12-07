# SECTIOn -2 Scatter Plot
data <- read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)
mouse_data = read.csv("bimm143_05_rstats/feature_counts.txt",header = TRUE, sep = "\t")
barplot(mouse_data$Count, horiz = TRUE,names.arg =mouse_data$Feature)


