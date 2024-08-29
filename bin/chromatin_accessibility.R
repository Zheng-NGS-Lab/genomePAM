

library(tidyverse)

# projDir = "~/miao/GenomePam"
# setwd(projDir)
# Set current working directory

setwd()
#df <-read.csv("data/df.csv")
# set 5m bins -------------------------------------------------------------
binsize=5e6


bin <- data.frame(chr = paste0("chr",c(1:22,"X","Y")), 
                  end = c(249250621,243199373,198022430,191154276,180915260,171115067,
                          159138663,146364022,141213431,135534747,135006516,133851895,
                          115169878,107349540,102531392,90354753,81195210,78077248,
                          59128983,63025520,48129895,51304566,155270560,59373566) ) %>%

# redirect from off-targets table
bin <- commandArgs(TRUE)[1]  
mutate(num = ceiling(end/binsize))

for (i in 1:nrow(bin)){
  tmp <- data.frame(chr = rep(bin$chr[i], bin$num[i]),
                    chrbin = 1:bin$num[i]) %>%
    mutate(start = (chrbin-1)*binsize+1,
           end = chrbin*binsize) %>%
    mutate(end = ifelse(chrbin==max(chrbin),bin$end[i],end))
  bin_df = rbind(bin_df,tmp)
  rm(tmp)
}
rm(bin,i)
bin_df$binnum = 1:nrow(bin_df)
save(bin_df, file = "data/bin_df.Rdata")


# read pam raw sequencing data --------------------------------------------
file_name = list.files("data/rawseq",full.names = T)
pam = data.frame()
for (i in 1:length(file_name)){
  tmp <- read_delim(file_name[i], delim = "\t", escape_double = FALSE,trim_ws = TRUE,
                    col_types = cols(Site_GapsAllowed.Sequence = col_character(), 
                                     Site_GapsAllowed.Length = col_number(), 
                                     Site_GapsAllowed.Score = col_number(), 
                                     Site_GapsAllowed.Substitutions = col_number(), 
                                     Site_GapsAllowed.Insertions = col_number(), 
                                     Site_GapsAllowed.Deletions = col_number(), 
                                     Site_GapsAllowed.Strand = col_character(), 
                                     Site_GapsAllowed.Start = col_double(), 
                                     Site_GapsAllowed.End = col_double()) ) %>% suppressWarnings()
  problems(tmp) %>% print
  pam  <- plyr::rbind.fill(pam, tmp)
  rm(tmp)
}
#table(str_detect(pam$`#BED_Chromosome`,"chrchr"))
#table(str_detect(pam$BED_Name,"chrchr"))
#table(str_detect(pam$BED_Site_Name,"chrchr"))
#table(str_detect(pam$BED_Site_Chromosom,"chrchr"))

pam$`#BED_Chromosome` = str_replace(pam$`#BED_Chromosome`,"chrchr","chr")
pam$BED_Name = str_replace(pam$BED_Name,"chrchr","chr")
pam$BED_Site_Name = str_replace(pam$BED_Site_Name,"chrchr","chr")
pam$BED_Site_Chromosome = str_replace(pam$BED_Site_Chromosom,"chrchr","chr")
  
pam <- pam %>% 
  mutate(cell_exp = case_when(
    str_detect(Filename,"Lib009") ~ "HepG2_1", str_detect(Filename,"Lib092") ~ "HepG2_2", str_detect(Filename,"Lib093") ~ "HepG2_3",
    str_detect(Filename,"Lib028") ~ "Huh7_1",  str_detect(Filename,"Lib029") ~ "Huh7_2", str_detect(Filename,"Lib030") ~ "Huh7_3",
    str_detect(Filename,"Lib034") ~ "HeLa_1",  str_detect(Filename,"Lib035") ~ "HeLa_2", str_detect(Filename,"Lib036") ~ "HeLa_3",
    str_detect(Filename,"Lib037") ~ "T293_1",  str_detect(Filename,"Lib038") ~ "T293_2", str_detect(Filename,"Lib039") ~ "T293_3"))
rm(file_name,i)
save(pam, file = "data/pam.Rdata")



# reads proportion by bins ------------------------------------------------
pam <- pam %>% 
  rename(chr = "#BED_Chromosome",  
         bed_start = "BED_Min.Position",
         bed_end = "BED_Max.Position", 
         reads = "total.sum",
         NumSubstitutions = "Site_SubstitutionsOnly.NumSubstitutions") %>%
  mutate(chrbin = ceiling(bed_start/binsize)) %>%
  select(chr,chrbin,bed_start,bed_end,reads,NumSubstitutions,cell_exp) %>%
  filter(reads >= 3,
         NumSubstitutions <= 7,
         chr %in% paste0("chr",c(1:22,"X","Y")))


pam_bin = data.frame()
for (i_cell in unique(pam$cell_exp)){
  
  tmp <- pam %>% filter(cell_exp == i_cell) 
  
  tmp <- merge(bin_df, tmp, by=c("chr","chrbin"), all.x = T) %>% 
    group_by(cell_exp,chr,chrbin,start,end,binnum) %>% 
    summarise(bin.reads = sum(reads, na.rm = T)) %>% ungroup() %>%
    mutate(cell.reads = sum(bin.reads, na.rm = T),
           pct = bin.reads/cell.reads,
           cell_exp = i_cell) %>%
    arrange(binnum)
  pam_bin = rbind(pam_bin,tmp)
  rm(tmp)
}



# consistency between experiments -----------------------------------------
cor_df <- pam_bin %>% 
  mutate(
    exp = case_when(
      str_detect(cell_exp,'_1') ~ 'Experiment 1',
      str_detect(cell_exp,'_2') ~ 'Experiment 2',
      str_detect(cell_exp,'_3') ~ 'Experiment 3'),
    cell = case_when(
      str_detect(cell_exp,'T293') ~ 'T293',
      str_detect(cell_exp,'HepG2') ~ 'HepG2',
      str_detect(cell_exp,'HeLa') ~ 'HeLa',
      str_detect(cell_exp,'Huh7') ~ 'Huh7')) %>%
  select(binnum,cell,exp,pct) %>%
  spread(exp,pct)

library(PerformanceAnalytics)

for (cell_line in c("HeLa",  "HepG2", "Huh7",  "T293")){
  
  tmp <-  cor_df %>% filter(cell==cell_line) %>% select('Experiment 1','Experiment 2','Experiment 3') %>% as.matrix()
  res <- Hmisc::rcorr(tmp)
  
  pdf(file = paste0("output/sup_cor_",cell_line,'.pdf'), height=5, width=5)
  
  chart.Correlation(tmp, histogram=F, pch=1)
  
  dev.off()
  
}


# PAM proportion across genome --------------------------------------------

# adjust proportion by T293 cell lines.
pamT293 <- pam_bin %>% 
  filter(str_detect(cell_exp,"T293")) %>% group_by(binnum) %>% 
  summarise(T293pct = sum(pct, na.rm = T)/3) %>% ungroup() 

pam_bin <- merge(pam_bin,pamT293,by="binnum") %>%
  mutate(pct_byt293 = ifelse(!(pct==0 | T293pct==0),log2(pct/T293pct), 0)) %>%
  filter(!str_detect(cell_exp,"T293")) 


x_break <- bin_df %>% filter(chr!="chrY") %>%  group_by(chr) %>% 
  summarize(bin = max(binnum)) %>% ungroup() %>% 
  mutate(chr = "")  

x_lab <- bin_df %>% filter(chr!="chrY") %>%  group_by(chr) %>% 
  summarize(bin = min(binnum)+((max(binnum)-min(binnum))/2) ) %>% ungroup() %>% 
  mutate(chr = str_replace(chr,"chr","")) 

x_df <- rbind(x_break,x_lab) %>% arrange(bin)

p <- pam_bin %>% filter(chr!="chrY") %>%
  ggplot(aes(x=binnum, y=pct_byt293)) + 
  facet_grid(cell_exp ~.,scales = "free_x") +
  ggh4x::stat_difference(aes(ymin = 0, ymax = pct_byt293))+
  geom_vline(xintercept = x_break$bin, linetype=3,  color = "black")+
  geom_hline(yintercept = 0, linetype=1,  color = "grey")+
  scale_x_continuous(expand = c(0.01,0),breaks = x_df$bin,
                     name = "Chromosome",labels = x_df$chr)+
  scale_y_continuous(breaks = c(-3,0,3), limits = c(-3,3))+
  scale_fill_manual(values=c("#d73027","#4575b4","black"))+
  theme_classic(base_size = 11)+
  theme(axis.text = element_text(color="black", size = 9),
        axis.ticks.x = element_blank(),
        legend.position = "none",)+
  labs(y="GUIDE-Seq read proportion, normalized by HEK293T cells, log2")

ggsave("output/Figure_5.pdf", plot=p, width=10, height=8)


writeLines(capture.output(sessionInfo()), "output/sessionInfo.txt")













