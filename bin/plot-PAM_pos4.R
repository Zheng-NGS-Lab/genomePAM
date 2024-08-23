
# An R script to visualize the dominat 4-base PAM base combinations.
# - Only the identifiedOfftargetsFile.txt file from the GUIDE-Seq output is needed. 
# - Optional: add 2nd argument for 'runDate', which will be used for output PDF names. Default: 000000
# - The 4-base position is automatically detected.
# example:
#	Rscript plot-PAM_pos4.R Lib091_identifiedOfftargets.txt

options(width=204)
library(ggseqlogo)
library(ggplot2)
OTfile = commandArgs(TRUE)[1]
runDate = commandArgs(TRUE)[2]
if(!exists('runDate')){runDate='000000'} 
cat('runDate: ', runDate)

    plot_4_PAM_pos = function(IdentifiedOfftargetsFile, runDate){
	LibID=sub("_identifiedOfftargets.*", "", sub(".*/","",IdentifiedOfftargetsFile)); cat("LibID: ", LibID, "\n")
	d1 = read.table(IdentifiedOfftargetsFile, sep='\t', stringsAsFactors=F, colClasses = c("character"), comment.char = "", header=T)
	target=d1$TargetSequence[1]
	spacer=sub('N{2,}', '', target)
	PAMcand=sub(spacer,'', target)
	lenPAM = nchar(PAMcand)
	if (nrow(d1) !=0){
		if (substr(target, 1, nchar(spacer)) == spacer){
			pamDir=3; d1$pam = substr(d1$Site_SubstitutionsOnly.Sequence, 1+nchar(spacer), nchar(target)) # 3' PAM
		}else{	pamDir=5; d1$pam = substr(d1$Site_SubstitutionsOnly.Sequence, 1, nchar(target)-nchar(spacer))} # 5' PAM
		cat("PAM Direction:", pamDir, "\n")
	d2 = subset(d1, pam !="")
	d2$PMMM = ifelse(d2$Site_SubstitutionsOnly.NumSubstitutions==0, 'PM', 
			ifelse(d2$Site_SubstitutionsOnly.NumSubstitutions>0 & d2$Site_SubstitutionsOnly.NumSubstitutions<7, 'MM', ''))
	for (mType in c('PM', 'MM')){
		dataPMMM = subset(d2, PMMM==mType)
		nsites = nrow(dataPMMM)
		if (nsites !=0){
		    pos_start = 1; pos_pam = NA; pam_d=NA;  peak=1;
		    for (i in 1:(lenPAM-3)){
			    pam_i = c(rep(substr(dataPMMM$pam, i, i+3), dataPMMM$bi.sum.mi))
			    pam_tabi = table(pam_i); pam_maxi = max(pam_tabi); sum_im = sum(table(pam_i)); 
			    cat('i:', i, 'max cnt: ', pam_maxi, '; sum: ', sum_im, 'pam_table: ', pam_tabi, '\n');
			    if(pam_maxi > peak){peak = pam_maxi; pos_start=i; pam_d = data.frame(pam_tabi)}
		    }
		    pam_d$pam_i = as.character(pam_d$pam_i)
		    pam_d$x1 = substr(pam_d$pam_i,1,1)
		    pam_d$x2 = substr(pam_d$pam_i,2,2)
		    pam_d$x3 = substr(pam_d$pam_i,3,3)
		    pam_d$x4 = substr(pam_d$pam_i,4,4)
		    for (x1 in c('A', 'G', 'C', 'T')){for (x2 in c('A', 'G', 'C', 'T')){for (x3 in c('A', 'G', 'C', 'T')){for (x4 in c('A', 'G', 'C', 'T')){
			    pam_d=rbind(pam_d, c(paste0(x1, x2, x3, x4), 0, x1, x2, x3, x4))}}}}
		    pam_d = pam_d[!duplicated(pam_d$pam_i),]
		    pam_d$Fraction = as.numeric(pam_d$Freq)/sum(as.numeric(pam_d$Freq))
		    pos_str = ifelse(pamDir==5, pos_start - lenPAM -1, pos_start)
		    pam_d$x1 = sub('^', paste0('Position ', pos_str, ': '), pam_d$x1)
		    pam_d$x2 = sub('^', paste0('Position ', pos_str+1, ': '), pam_d$x2)
		pdf(file=paste0(runDate, "_", LibID, "_PAM_", mType, "_", pos_str, "-", pos_str+3, '.pdf'))
			myplot = ggplot(pam_d, aes(x=x3,y=x4, fill=Fraction)) +
			    geom_tile(color = "grey")+
			    facet_grid(x2 ~ x1)+
			    scale_fill_gradientn(colours=rev(heat.colors(100)))+
			    scale_x_discrete(expand=c(0,0))+
			    scale_y_discrete(expand=c(0,0))+
			    coord_fixed()+
			    labs(title=paste0(runDate,"_",LibID, "_PAM \n(", nsites," ", mType, " sites) positions: ", pos_str," to ", pos_str+3)
				, x=paste0("Position ", pos_str+2), y = paste0("Position ", pos_str+3), hjust = 0.5, vjust = 0.5) +
			    theme(text = element_text(family = "Courier"), strip.placement = "outside", plot.title = element_text(hjust = 0))
				
				print(myplot)
			    dev.off()
	}}}}

plot_4_PAM_pos(IdentifiedOfftargetsFile=OTfile, runDate=runDate)

