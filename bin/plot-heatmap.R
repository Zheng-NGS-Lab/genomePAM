options(width=204)
library(plyr)
library(ggseqlogo)
library(ggplot2)
# options(repos = c(CRAN = "https://cloud.r-project.org"))
# install.packages("reshape2")
# library(reshape)
library(patchwork)

# Command line arguments
IdentifiedOfftargetFile = commandArgs(TRUE)[1]
spacer = commandArgs(TRUE)[2]
target = commandArgs(TRUE)[3]
LibID = commandArgs(TRUE)[4]
umitagged_count = commandArgs(TRUE)[5]
consolidated_count = commandArgs(TRUE)[6]

# Load and process data
d1 = read.table(IdentifiedOfftargetFile, sep='\t', stringsAsFactors=F, colClasses=c("character"), comment.char="", header=T, fill=TRUE)
head(d1,3)
lenSpacer = nchar(spacer)
lenTarget = nchar(target)
if (substr(target, 1, lenSpacer) == spacer) {
    d1$pam = substr(d1$Site_SubstitutionsOnly.Sequence, 1+lenSpacer, lenTarget)  # 3' PAM
} else {
    d1$pam = substr(d1$Site_SubstitutionsOnly.Sequence, 1, lenTarget - lenSpacer)  # 5' PAM
}

pd_unsorted = subset(d1, pam !="")
pd = pd_unsorted[order(as.numeric(pd_unsorted$bi.sum.mi), decreasing=TRUE), ]
# Perfect Match
pd_PM = subset(pd, Site_SubstitutionsOnly.NumSubstitutions==0)
# Mismatch
pd_MM = subset(pd, Site_SubstitutionsOnly.NumSubstitutions>0)

#####Generate plots
# Initialize cumulative alignment read counts
colors <- c("PM" = "#742554", "MM" = "darkgrey")
p_cum <- ggplot() +
    labs(title=paste0(LibID, ": Cumulative alignments' read counts"), x="Site #", y="Aligment count (cumulative)", color="Legend", subtitle = paste0("Umitagged read count: ", umitagged_count, "; Unique read count: ", consolidated_count)) + 
    scale_color_manual(values = colors)

if (nrow(pd_MM)==0 & nrow(pd_PM)==0){
    p_cum <- p_cum + labs(caption="No PM or MM reads.")
}

# Mismatch
if (nrow(pd_MM) !=0){
    rownames(pd_MM) <- 1:nrow(pd_MM)
    # Cumulative alignment read counts
    p_cum <- p_cum + 
        geom_point(data=pd_MM, aes(x=as.numeric(rownames(pd_MM)), y=cumsum(bi.sum.mi), color="MM"), size = 1.8) + 
        geom_line(data=pd_MM, aes(x=as.numeric(rownames(pd_MM)), y=cumsum(bi.sum.mi), color="MM", group = 1), size = 1.2)    
    
    # Sequence logo for mismatch
    n.reads_MM = sum(as.numeric(pd_MM$bi.sum.mi))
    n.sites_MM = length(pd_MM$bi.sum.mi)
    pam_MM = c(rep(pd_MM$pam, pd_MM$bi.sum.mi))
    p_MM = ggplot() + 
        geom_logo(pam_MM) + 
        labs(title=paste0(LibID, ": Mismatch"), subtitle=paste0("Num of reads: ", n.reads_MM, "; Num of sites: ", n.sites_MM)) + 
        theme_logo() + ylim(0,2)
} else {
    # Set n.reads & n.sites to 0
    n.reads_MM = 0
    n.sites_MM = 0
}

# Perfect match
if (nrow(pd_PM) !=0){
    rownames(pd_PM) <- 1:nrow(pd_PM)
    # Cumulative alignment read counts
    p_cum <- p_cum + 
        geom_point(data=pd_PM, aes(x=as.numeric(rownames(pd_PM)), y=cumsum(bi.sum.mi), color="PM"), size = 1.8) + 
        geom_line(data=pd_PM, aes(x=as.numeric(rownames(pd_PM)), y=cumsum(bi.sum.mi), color="PM", group = 1), size = 1.2) 
    
    # Sequence logo for perfect match
    n.reads_PM = sum(as.numeric(pd_PM$bi.sum.mi))
    n.sites_PM = length(pd_PM$bi.sum.mi)
    pam_PM = c(rep(pd_PM$pam, pd_PM$bi.sum.mi))
    p_PM = ggplot() + 
        geom_logo(pam_PM) + 
        labs(title=paste0(LibID, ": Perfect Match"), subtitle=paste0("Num of reads: ", n.reads_PM, "; Num of sites: ", n.sites_PM)) + 
        theme_logo() + ylim(0,2)
} else {
    # Set n.reads & n.sites to 0
    n.reads_PM = 0
    n.sites_PM = 0
}

#####Generate Heatmap
# Aggregate data for heatmap
heatmap_data <- ddply(pd, .(pam), summarize, count = sum(bi.sum.mi))

# Convert to matrix format suitable for heatmap plotting
heatmap_matrix <- acast(heatmap_data, pam ~ ., value.var="count")

# Transform data to long format for heatmap
heatmap_m <- melt(heatmap_matrix)

# Plot heatmap
g <- ggplot(heatmap_m, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = topo.colors(100)) + 
  labs(fill = "Count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("PAM Sequences") +
  ylab("Count")

#####Patch subplots and save as PDF
if (exists("p_PM") & exists("p_MM")){
    p_all <- p_cum / p_PM / p_MM / g
    plot_height <- 12
} else if (exists("p_PM")) {
    p_all <- p_cum / p_PM / g
    plot_height <- 9
} else if (exists("p_MM")) {
    p_all <- p_cum / p_MM / g
    plot_height <- 9
} else {
    p_all <- p_cum / g
    plot_height <- 6
}
ggsave(paste0(LibID, '_visualize.pdf'), width=8, height=plot_height)

#####Export read stats in CSV
stat_csv <- paste0('.', LibID, '_stats.csv')
# header <- c("LibID", "Umitagged_read_count", "Unique_read_count", "PM_n_reads", "PM_n_sites", "MM_n_reads", "MM_n_sites")
stats <- c(LibID, umitagged_count, consolidated_count, n.reads_PM, n.sites_PM, n.reads_MM, n.sites_MM)
# write(paste(header, collapse = ","), file=stat_csv, append=F)
write(paste(stats, collapse = ","), file=stat_csv, append=T)