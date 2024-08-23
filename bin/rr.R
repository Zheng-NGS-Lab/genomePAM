options(width=204)
library(ggseqlogo)
library(ggplot2)
library(data.table)
library(grid)

###############################
IdentifiedOfftargetsFile=commandArgs(TRUE)[1]
bs<-c('A','C','G','T')
ss<-apply(expand.grid(bs,bs)[,2:1],1,paste,collapse='')
paste0(substr(ss,2,2),substr(ss,1,1))->ss2

data.table(two=ss,top=substr(ss,1,1),bot=substr(ss,2,2))->heng
melt(heng,id.vars='two')->heng
heng[,two:=factor(two,ss)]

px<-ggplot(heng,aes(x=two,y=variable,fill=value))+
    geom_tile(colour='white')+
    geom_text(aes(label=value))+
    scale_x_discrete(expand=c(0,0))+
    scale_y_discrete(expand=c(0,0))+
    coord_equal()

copy(heng)->z2
z2[,variable:=factor(variable,c('bot','top'))]
z2[,two:=factor(two,rev(ss),labels=rev(ss2))]

py<-ggplot(z2,aes(y=two,x=variable,fill=value))+
    geom_tile(colour='white')+
    geom_text(aes(label=value))+
    scale_x_discrete(expand=c(0,0))+
    scale_y_discrete(expand=c(0,0))+
    coord_equal()
##############################

##############################
FFF<-function(IdentifiedOfftargetsFile, runDate){
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
			
			pam_d$b12<-substr(pam_d$pam_i,1,2)
			pam_d$b34<-substr(pam_d$pam_i,3,4)
			
			pam.inter<-pam_d[,c('b12','b34','Fraction')]
			pam.inter$b34<-factor(pam.inter$b34,ss)
            pam.inter$b12<-factor(pam.inter$b12,rev(ss2))
			

            pp<-ggplot(pam.inter,aes(x=b34,y=b12,fill=Fraction))+
            geom_tile(colour='white')+
            scale_x_discrete(expand=c(0,0),position='bottom')+
            scale_y_discrete(expand=c(0,0),position='right')+
            scale_fill_gradientn(colours=rev(heat.colors(100)))+
            coord_equal()+
            theme(axis.ticks=element_blank(),axis.text=element_blank(),axis.title=element_blank())
            
            
            ggplotGrob(pp)->mat
            ggplotGrob(px)->matx
            ggplotGrob(py)->maty
            mat$grobs[[11]]<-matx$grobs[[6]]
            mat$grobs[[3]]<-maty$grobs[[6]]
            mat$heights[[7]]<-unit(2,'cm')
            mat$widths[[6]]<-unit(2,'cm')
			
		pdf(file=paste0(runDate, "_", LibID, "_PAM_", mType, "_", pos_str, "-", pos_str+3, '.pdf'))
        grid.draw(mat)
		dev.off()
	}}}
}

#IdentifiedOfftargetsFile<-'example2.tsv.txt'
#runDate<-'00'
FFF(IdentifiedOfftargetsFile,'today')
###

