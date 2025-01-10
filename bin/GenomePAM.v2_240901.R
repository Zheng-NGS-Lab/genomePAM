# GenomePAM
# An R script to analyze and visualize the PAM sequences based on the identifiedOfftargets.txt from GUIDE-seq.
# Only the identifiedOfftargetsFile.txt file from the GUIDE-Seq output is needed.
# Computation steps:
# 1. Identify the strongest single base among the PAM bases and its position.
# 2. Extend PAM bases from Position A one base towards left and one base towards right
# 3. Calculate the difference among the extended sequences in either direction, Left and Right, respectively and record the two P values.
# 4. Compare the P values. The more significant one is selected as the enriched PAM and used as the basis for the next round of extension.
# 5. Repeat Step 2 to 4 until the ends of candidate PAM bases.
# Usage example:
# Rscript GenomePAM.R "/my/path/to/example/Lib001_identifiedOfftargets.txt" "/tools/repo/genomePAM/"
# , where the 1st argument after GenomePAM.R is the identifiedOfftargets.txt to be analyzed and the 2nd the path where the GenomePAM.R and background data (Rep-1, Rep-2, etc.) are stored [default: "/tools/repo/GenomePAM/"].

# Read command line arguments
OTfile <- commandArgs(TRUE)[1]
bgPAMpath <- commandArgs(TRUE)[2]
PAMlen < - commandArgs(TRUE)[3]
PAMpos <- commandArgs(TRUE)[4]
#PCVfile <- commandArgs(TRUE)[5]

# Set default path if bgPAMpath is not provided
if (is.na(bgPAMpath)) {
  bgPAMpath <- '/tools/repo/GenomePAM/'
}

# Load necessary libraries
library(readr)
library(dplyr)
library(gt)
library(glue)
library(purrr)
library(stringr)

# Extract library ID from the file name
LibID <- sub("_identifiedOfftargets.*", "", sub(".*/", "", OTfile))

# Read the input file
d1 <- read.table(OTfile, sep = '\t', stringsAsFactors = FALSE, colClasses = "character", comment.char = "", header = TRUE, fill = TRUE)

# Extract target sequence and spacer
target <- d1$TargetSequence[1]
spacer0 <- sub('N{2,}', '', target)

# Determine PAM direction and extract PAM sequence
if (substr(target, 1, nchar(spacer0)) == spacer0) {
    pamDir <- 3
    d1$pam <- substr(d1$Site_SubstitutionsOnly.Sequence, 1 + nchar(spacer0), nchar(target)) # 3' PAM
} else {
    pamDir <- 5
    d1$pam <- substr(d1$Site_SubstitutionsOnly.Sequence, 1, nchar(target) - nchar(spacer0)) # 5' PAM
}

# Clean up spacer and PAM candidate sequences
spacer <- sub('.*[^GCTA]', '', spacer0) # in case mixed spacers (i.e. N, W, S, Y...) are used in spacer head (in the first few bases over 20 nt)
PAMcand <- sub(spacer0, '', target)
lenPAM <- nchar(PAMcand)

# Read background / genome PAM file
bgPAMfile <- paste0(bgPAMpath, 'PAM', pamDir, '_', spacer, '.cnt')
bgPAM <- read.table(bgPAMfile, sep = '\t', stringsAsFactors = FALSE, header = TRUE)
names(bgPAM) <-  sub('Freq', 'genome', names(bgPAM))

# Filter data and classify PM/MM
d2 <- subset(d1, pam != "")
d2$PMMM <- ifelse(d2$Site_SubstitutionsOnly.NumSubstitutions == 0, 'PM',
                  ifelse(d2$Site_SubstitutionsOnly.NumSubstitutions > 0 & d2$Site_SubstitutionsOnly.NumSubstitutions < 7, 'MM', ''))

# Subset data based on match type
dataPM <- subset(d2, PMMM == 'PM')
nReads_PM = sum(as.numeric(dataPM$bi.sum.mi))
nedits <- nrow(dataPM)

if (nedits != 0) {
	pam_dr0 <- NULL
	pam_ds0 <- NULL
	
	for (pos in 1:lenPAM) {
	    for (len in 1:(lenPAM - pos + 1)) {
		pam_r <- rep(substr(dataPM$pam, pos, pos + len - 1), dataPM$bi.sum.mi) # reads
		pam_tabr <- data.frame(table(pam_r))
		pam_tabr$pos <- pos
		pam_tabr$pam2 <- paste0(pos, pam_tabr$pam_r)
		pam_tabr$tot <- sum(table(pam_r))
		pam_dr0 <- rbind(pam_dr0, pam_tabr)
		
		pam_s <- substr(dataPM$pam, pos, pos + len - 1) # edits
		pam_tabs <- data.frame(table(pam_s))
		pam_tabs$pos <- pos
		pam_tabs$pam2 <- paste0(pos, pam_tabs$pam_s)
		pam_tabs$tot <- sum(table(pam_s))
		pam_ds0 <- rbind(pam_ds0, pam_tabs)
	    }
	}
	
	pam_d0 <- merge(pam_dr0[, !names(pam_dr0) %in% c('PAMlen', 'pos', 'pam_r')], 
			pam_ds0[, !names(pam_ds0) %in% c('PAMlen', 'pos', 'pam_s')], 
			by = 'pam2', suffix = c('.r', '.s'))
	names(pam_d0) <- sub('Freq.r', 'reads', names(pam_d0))
	names(pam_d0) <- sub('Freq.s', 'edits', names(pam_d0)) # edits: edited genomic sites
	
	pam_m <- merge(bgPAM, pam_d0, by = 'pam2', all.x = TRUE)
	pam_m[is.na(pam_m)] <- 0
	pam_m$pctg <- pam_m$genome / pam_m$tot
	pam_m$pctr <- pam_m$reads / pam_m$tot.r
	pam_m$pcts <- pam_m$edits / pam_m$tot.s
	pam_m$enrr <- pam_m$pctr / pam_m$pctg # enrich fold, by reads
	pam_m$enrs <- pam_m$pcts / pam_m$pctg # enrich fold, by site
	pam_m$enr <- (pam_m$enrr + pam_m$enrs) / 2 # enrich fold, by site
	pam_m$edit2 <- pam_m$edits + round(pam_m$reads / max(pam_m$reads) * max(pam_m$edits)) # In addition to site, take into account number of reads, which is scaled to the range of size
	pam_m$enrs2 <- (pam_m$edit2 / pam_m$tot.s) / 2 / pam_m$pctg # enrich fold, by site
	pam_m$percent <- round(pam_m$edits / pam_m$genome * 100, 2)

	# Calculation of PAM cleavage value (PCV)
	pam_m$PCV0 <- pam_m$enrr
	pam_m$PCV0 <- ifelse(pam_m$enrr < 0.01 | is.na(pam_m$enrr), 0.01, pam_m$enrr) # make min as 0.01
        pam_m$PCV <- pam_m$PCV0/max(pam_m$PCV0)*100 # scale to 100
        pam_m$PCV2 <- log2(pam_m$PCV) # log2 transformation
        pam_m$PCV3 <- pam_m$PCV2 - min(pam_m$PCV2) # make lower boundary of log2 transformed PCV as 0, instead of negative
        pam_m$Relative.PCV <- pam_m$PCV3/max(pam_m$PCV3) # Relative to max (1)
	
	# Further processing and statistical analysis
	p_cutoff <- 0.0001
	enrsCutoff <- 1.5
	p1s <- NULL
	enrs1s <- NULL
	## Sub PAMlen with 4 
	pam1 <- subset(pam_m, PAMlen == 1)
	
	for (i in 1:lenPAM) {
	    p1 <- chisq.test(pam1[pam1$pos == i, c('edits', 'genome')])$p.value
	    p1s <- c(p1s, p1)
	    enrs1 <- max(pam1[pam1$pos == i, 'enrs2'])
	    enrs1s <- c(enrs1s, enrs1)
	}
	
	pos0 <- which(p1s == min(p1s))
	posSig <- which(enrs1s >= enrsCutoff)
	posMin <- min(posSig)
	posMax <- max(posSig)
	
	pam1 <- pam1[order(-(pam1$enrr > 1 & pam1$enrs > 1), -pam1$percent), ]
	pam1s <- subset(pam1, pos == pos0)
	pvalue0.r <- chisq.test(pam1s[, c('reads', 'genome')])$p.value
	pvalue0.s <- chisq.test(pam1s[, c('edits', 'genome')])$p.value
	
	# Extension iteration and Enrichment calculations
	posR <- pos0 + 1
	posL <- pos0 - 1
	sigPAM <- subset(pam1s, enrr > 1 & enrs > 1)
	topPAM <- sigPAM$pam_hg38
	sigPAM$P_value <- min(pvalue0.r, pvalue0.s)

	for (i in 1:lenPAM) {
	if (posL >= 1) {
	    pam_m$prePAM <- substr(pam_m$pam_hg38, 2, nchar(pam_m$pam_hg38))
	    data.L <- subset(pam_m, (pos == posL) & (PAMlen == i + 1) & (prePAM %in% topPAM))
	    pvL <- 1
	    data.Li <- NULL
	    
	    for (topPAMi in topPAM) {
		dataLi <- subset(data.L, prePAM == topPAMi)
		if (nrow(dataLi) >= 2) {
		    pvLi <- if (sum(dataLi$edit2) < 50) {
			fisher.test(dataLi[, c('edit2', 'genome')], workspace = 2e9)$p.value
		    } else {
			chisq.test(dataLi[, c('edit2', 'genome')])$p.value
		    }
		    data.L$P_value[data.L$prePAM == topPAMi] <- pvLi
		    if (pvLi < p_cutoff) {
			data.Li <- rbind(data.Li, data.L[data.L$prePAM == topPAMi, ])
		    }
		    if (pvLi < pvL) {
			pvL <- pvLi
		    }
		}
	    }
	}
	
	if (posR <= lenPAM) {
	    pam_m$prePAM <- substr(pam_m$pam_hg38, 1, posR - posL - 1)
	    data.R <- subset(pam_m, (pos == posL + 1) & (PAMlen == i + 1) & (prePAM %in% topPAM)) # note: use posL here for extending to the Right
	    pvR <- 1
	    data.Ri <- NULL
	    
	    for (topPAMi in topPAM) {
		dataRi <- subset(data.R, prePAM == topPAMi)
		if (nrow(dataRi) >= 2) {
		    pvRi <- if (sum(dataRi$edit2) < 50) {
			fisher.test(dataRi[, c('edit2', 'genome')], workspace = 2e9)$p.value
		    } else {
			chisq.test(dataRi[, c('edit2', 'genome')])$p.value
		    }
		    data.R$P_value[data.R$prePAM == topPAMi] <- pvRi
		    if (pvRi < p_cutoff) {
			data.Ri <- rbind(data.Ri, data.R[data.R$prePAM == topPAMi, ])
		    }
		    if (pvRi < pvR) {
			pvR <- pvRi
		    }
		}
	    }
	}
	
	cat("iteration: ", i, "; pv_L = ", pvL, "; pv_R = ", pvR, "\n")

	if (posL >= 1 & pvL <= pvR) {
	    if (pvL < p_cutoff) {
		data.Li$pctg <- data.Li$genome / sum(data.Li$genome)
		data.Li$pctr <- data.Li$reads / sum(data.Li$reads)
		data.Li$pcts <- data.Li$edits / sum(data.Li$edits)
		data.Li$enr <- (data.Li$pcts / data.Li$pctg + data.Li$pctr / data.Li$pctg) / 2
		data.Li$P_value <- pvL
		posL <- posL - 1
		sigPAM <- rbind(sigPAM, subset(data.Li[, names(sigPAM)], enr > 1))
		topPAM <- data.Li$pam_hg38[data.Li$enr > 1]
	    } else if (posL > posMin) {
		posL <- posL - 1
		sigPAM <- rbind(sigPAM, data.L[, names(sigPAM)])
		topPAM <- data.L$pam_hg38
	    }
	} else {
	    if (posR <= lenPAM) {
		if (pvR < p_cutoff) {
		    data.Ri$pctg <- data.Ri$genome / sum(data.Ri$genome)
		    data.Ri$pcts <- data.Ri$edits / sum(data.Ri$edits)
		    data.Ri$enrs <- data.Ri$pcts / data.Ri$pctg
		    data.Ri$dirs <- ifelse(data.Ri$enrs > 1, 1, 0)
		    data.Ri$P_value <- pvR
		    posR <- posR + 1
		    sigPAM <- rbind(sigPAM, subset(data.Ri[, names(sigPAM)], enr > 1))
		    topPAM <- data.Ri$pam_hg38[data.Ri$enr > 1]
		} else if (posR < posMax) {
		    posR <- posR + 1
		    sigPAM <- rbind(sigPAM, data.R[, names(sigPAM)])
		    topPAM <- data.R$pam_hg38
		} else if (posL > posMin) {
		    posL <- posL - 1
		    sigPAM <- rbind(sigPAM, data.L[, names(sigPAM)])
		    topPAM <- data.L$pam_hg38
		}
	    }
	}
    }
}

# Insignificant PAM bases denoted as '.'.
sigPAM <- subset(sigPAM, P_value <= p_cutoff)
sigPAM$PAM <- sapply(1:nrow(sigPAM), function(i) {
  paste0(
    paste0(rep('.', sigPAM$pos[i] - 1), collapse = ''),
    sigPAM$pam_hg38[i],
    paste0(rep('.', lenPAM - sigPAM$pos[i] - sigPAM$PAMlen[i] + 1), collapse = '')
  )
})

# Order and filter sigPAM
sigPAM <- sigPAM[order(sigPAM$PAMlen, sigPAM$percent), ]
sigPAM$keep <- 1

for (i in 1:nrow(sigPAM)) {
  if (sigPAM$percent[i] < max(sigPAM$percent[1:(i - 1)], na.rm = TRUE)) {
    sigPAM$keep[i] <- 0
  }
}

# Select and print the output
outPAM <- sigPAM[sigPAM$keep == 1 & sigPAM$genome >= 3, c('PAM', 'genome', 'edits', 'percent', 'P_value')]
print(outPAM)

# Write output to files
write.table(outPAM, paste0(LibID, '_GenomePAM.txt'), row.names = FALSE, quote = FALSE, sep = '\t')
dropVars <- c('PCV0', 'PCV2', 'PCV3')
write.table(pam_m[, !(names(pam_m) %in% dropVars)], paste0(LibID, '_GenomePAM_raw.txt'), row.names = FALSE, quote = FALSE, sep = '\t')

#=== plot GenomePAM table ===#

color_pam <- function(sequence) {
  sequence <- gsub("G", "<span style='color:orange; font-family: monospace;'>G</span>", sequence)
  sequence <- gsub("A", "<span style='color:green; font-family: monospace;'>A</span>", sequence)
  sequence <- gsub("T", "<span style='color:red; font-family: monospace;'>T</span>", sequence)
  sequence <- gsub("C", "<span style='color:blue; font-family: monospace;'>C</span>", sequence)
  sequence <- gsub("\\.", "<span style='color:black; font-family: monospace;'>.</span>", sequence)
  return(sequence)
}

# Data manipulation
df <- outPAM %>%
  mutate(
    PAM = sapply(PAM, color_pam),
    P_value = sprintf("%.1e", P_value)
  )

# Create the gt table
gt_table <- df %>%
  gt() %>%
  fmt_markdown(columns = PAM) %>%
  text_transform(
    locations = cells_body(columns = percent),
    fn = function(x) {
      x <- as.numeric(x)
      max_percent <- max(df$percent, na.rm = TRUE)
      formatted_numbers <- sprintf("%.1f", x)
      lapply(seq_along(x), function(i) {
        percent_bar <- round((x[i] / max_percent) * 100, 1)
        glue::glue(
          "<div style='position: relative; height: 4px;'>
            <div style='background-color: skyblue; height: 4px; width: {percent_bar}%;'></div> # nolint: line_length_linter.
            <div style='position: absolute; right: 0; bottom: 0; height: 4px; line-height: 4px; z-index: 1; color: black;'>{formatted_numbers[i]}</div> # nolint: line_length_linter.
          </div>"
        )
      })
    }
  ) %>%
  fmt_number(
    columns = P_value,
    decimals = 0
  ) %>%
  cols_align(
    align = "left",
    columns = percent
  ) %>%
  cols_label(
    PAM = ifelse(pamDir == 3, "PAM (3')", "PAM (5')"),
    genome = "Genome",
    edits = "Edits",
    percent = "Percent (%)",
    P_value = "P-value"
  ) %>%
  tab_style(
    style = cell_text(font = "monospace"),
    locations = cells_body(columns = PAM)
  ) %>%
  tab_options(
    data_row.padding = px(0.5),
    table.font.size = 12,
    heading.title.font.size = 12,
    heading.subtitle.font.size = 12
  ) %>%
  opt_table_font(font = "Courier New"
  ) %>%
  tab_header(
    title = md(LibID),
    subtitle = md(paste0("Num of perfect match reads: ", nReads_PM))
  )

gtsave(gt_table, paste0(LibID, "_GenomePAM_Tab.html"))
