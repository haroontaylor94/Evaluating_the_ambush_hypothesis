setwd("/Users/haroontaylor/Documents/University/Final_Year/Project/")

create_box_plot <- function(x, x2, x_names, xlab, ylab, scatterplot = TRUE){
    continuous_vector = c(x, x2)
    categorical_vector = c(rep(x_names[1], length(x)), 
                           rep(x_names[2], length(x2)))
    boxplot(continuous_vector ~ categorical_vector, outline = FALSE, staplelty = 2,
            ylab = ylab, xlab = xlab, border = rgb(0, 0.6, 0.8))
    if (scatterplot == TRUE){
        stripchart(continuous_vector ~ categorical_vector, vertical = TRUE, 
                   method = "jitter", add = TRUE, pch = 20, col = rgb(0.2, 0.6,
                                                                      0.8, 0.3))
    }
} 

create_histogram <- function(data, xlab){
    hist(data, main=NULL, ylab = 'Frequency', 
         xlab = xlab, col = rgb(0.2, 0.6, 0.8))
    shapiro_p = signif(shapiro.test(data)$p.value, 3)
    legend('topright', legend = paste('shapiro p = \n', shapiro_p), bty = "n", 
           cex = 0.7)
}

create_scatterplot <- function(x, y, opacity = 1){
    xlab = gsub("_", " ", deparse(substitute(x)))
    ylab = gsub("_", " ", deparse(substitute(y)))
    plot(x, y, xlab = xlab, ylab = ylab, pch = 20, col = rgb(0.2, 0.6, 0.8, 
                                                             opacity))
    abline(lm(y ~ x))
    spearman_rank = cor.test(x, y, method = "spearman")
    print(spearman_rank)
    spearman_rho = signif(spearman_rank$estimate, 3)
    if (spearman_rho <= 0){
        legend_pos = "topright"
    }
    else {
        legend_pos = "bottomright"
    }
    legend(legend_pos, legend = paste('rho = ', spearman_rho), 
           bty = "n")
}


calculate_percentage_excess <- function(excess, stdev, frameshift){
    percent_excess = sum(excess > 0)/length(excess)
    significant_percent_excess = sum(excess/stdev > 2)/length(excess)
    cat(paste("Percent of genomes with", frameshift, "OSC excess =", 
              percent_excess, "\n\n"))
    cat(paste("Percent of genomes with significant", frameshift, "OSC excess =", 
              significant_percent_excess, "\n\n"))
}


make_all_of_the_plots <- function(input_file, output_directory, opacity){
    df <- read.csv(input_file)
    setwd(output_directory)
    attach(df)
    
    pdf("sequence_length_histogram.pdf")
    create_histogram(Sequence_length, "Sequence length")
    dev.off()
    
    pdf("plus1_vs_plus2_boxplot.pdf")
    create_box_plot(Plus_1_OSC_percentage, Plus_2_OSC_percentage, c("+1", "+2"), 
                    "Frame shift", "OSC percentage")
    dev.off()
    
    pdf("plus1_plus2_histograms.pdf")
    par(mfrow = c(2,2))
    create_histogram(Plus_1_OSC_percentage, "+1 OSC percentage")
    create_histogram(Plus_2_OSC_percentage, "+2 OSC percentage")
    dev.off()
    
    wilcox_test <- wilcox.test(Plus_1_OSC_percentage, Plus_2_OSC_percentage, 
                               paired = TRUE)
    cat(paste('Wilcoxon p = ', wilcox_test$p.value, "\n\n"))
    
    OSC_excess = Total_OSC_percentage - Total_OSC_simulation_mean
    calculate_percentage_excess(OSC_excess, Total_OSC_simulation_stdev, "Total")
    
    Plus_1_excess = Plus_1_OSC_percentage - Plus_1_OSC_simulation_mean
    calculate_percentage_excess(Plus_1_excess, Plus_1_OSC_simulation_stdev, "+1")
    
    Plus_2_excess = Plus_2_OSC_percentage - Plus_2_OSC_simulation_mean
    calculate_percentage_excess(Plus_2_excess, Plus_2_OSC_simulation_stdev, "+2")
    
    OSC_percentage_vs_GC_content_lm <- lm(Total_OSC_percentage ~ GC_content + Sequence_length)
    cat("***OSC percentage vs GC content lm***\n")
    print(summary(OSC_percentage_vs_GC_content_lm))
    
    OSC_percentage_vs_GC3_content_lm <- lm(Total_OSC_percentage ~ GC3_content + Sequence_length)
    cat("***OSC percentage vs GC3 content lm***\n")
    print(summary(OSC_percentage_vs_GC3_content_lm))
    
    OSC_excess_vs_GC_content_lm <- lm(OSC_excess ~ GC_content + Sequence_length)
    cat("***OSC excess vs GC content lm***\n")
    print(summary(OSC_excess_vs_GC_content_lm))
    
    OSC_excess_vs_GC3_content_lm <- lm(OSC_excess ~ GC3_content + Sequence_length)
    cat("***OSC excess vs GC3 content lm***\n")
    print(summary(OSC_excess_vs_GC3_content_lm))
    
    pdf("OSC vs GC scatter.pdf")
    par(mfrow = c(2, 2))
    cat("***OSC percentage vs GC content***\n")
    create_scatterplot(GC_content, Total_OSC_percentage, opacity)
    cat("***OSC percentage vs GC3 content***\n")   
    create_scatterplot(GC3_content, Total_OSC_percentage, opacity)
    dev.off()
    
    OSC_excess_z_score <- OSC_excess / Total_OSC_simulation_stdev
    pdf("OSC excess vs GC scatter.pdf")
    par(mfrow = c(2, 2))
    cat("***OSC excess vs GC content***\n")
    create_scatterplot(GC_content, OSC_excess_z_score)
    cat("***OSC excess vs GC3 content***\n")
    create_scatterplot(GC3_content, OSC_excess_z_score)
    dev.off()
    
    pdf("OSC intervals vs GC scatter.pdf")
    par(mfrow = c(2, 2))
    cat("***OSC intervals vs GC content***")
    create_scatterplot(GC_content, OSC_gap_mean)
    cat("***OSC intervals vs GC3 content***")
    create_scatterplot(GC3_content, OSC_gap_mean)
    dev.off()
    
    detach(df)
    setwd("/Users/haroontaylor/Documents/University/Final_Year/Project/")
}

#Make plots
cat("******Table 4 data******\n\n")
make_all_of_the_plots("Results/Table_4_results.csv", "Graphs/Table_4/", 0.7)
cat("******Table 11 data******\n\n")
make_all_of_the_plots("Results/Table_11_results_2.csv", "Graphs/Table_11/", 0.3)

#Compare sequence length to OSC percentage
table_11_data <- read.csv("Results/table_11_results.csv")
table_4_data <- read.csv("Results/table_4_results.csv")

#Create boxplot of table 4 vs table 11 sequence length
pdf("Graphs/Boxplot_table_11_vs_table_4_seq_length.pdf")
create_box_plot(table_4_data$Sequence_length, table_11_data$Sequence_length,
                c("4", "11"), "Translation table", "Sequence length", 
                scatterplot = FALSE)
dev.off()
shapiro.test(table_4_data$Sequence_length)
shapiro.test(table_11_data$Sequence_length)
wilcox.test(table_4_data$Sequence_length, table_11_data$Sequence_length)

#Create scatterplots for gap mean vs gap stdev
pdf("Graphs/Gap_mean_vs_stdev_scatter.pdf")
par(mfrow = c(2, 2))
attach(table_4_data)
cat("******Table 4 interval mean vs stdev******")
create_scatterplot(OSC_gap_stdev, OSC_gap_mean)
detach(table_4_data)
attach(table_11_data)
cat("******Table 11 interval mean vs stdev******")
create_scatterplot(OSC_gap_stdev, OSC_gap_mean)
detach(table_11_data)
dev.off()


