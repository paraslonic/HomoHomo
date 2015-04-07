library(vioplot)
#setwd("/data4/bio/runs-manolov/ecoli_crohn/MiraAssemblies//HomoHomo/")

args = commandArgs(trailingOnly=T)
correctedFile = args[1]
infoFile = args[2]
outFile = args[3]

pdf(outFile, paper="a4r")
# C O R R E C T E D    S T A T
corrected = read.delim(correctedFile, head = F, as.is = T)
colnames(corrected) = c("contig", "pos", "evalue")
corrected_contig_pos = paste0(corrected$contig, "_", corrected$pos)
boxplot(corrected$evalue, breaks = 100, main = "blast e-value of corrected")

# I N D E L   I N F O
info = read.delim(infoFile, as.is = T)
info_contig_pos = paste0(info$contig,"_", info$pos)

nrow (info)
nrow (corrected)
old.par = par(mar = c(8,8,8,8))
pie(c(nrow (info),nrow (corrected)), labels = c("all indels", "corrected"), col = c("black", "orange"), cex = 2)
par(old.par)

# C O M P A R E 
info.crd = subset(info,  info_contig_pos %in% corrected_contig_pos)
nrow (info.crd)

# freq 
vioplot(info$altFreq, info.crd$altFreq, col = "darkolivegreen4",
        names= c("ALL INDEL","CORRECTED"))
title("ALT Variant Frequency")

# freq 
vioplot(info$DP, info.crd$DP, col = "darkolivegreen4", 
        names= c("ALL INDEL","CORRECTED"))
title("DEPTH OF COVERAGE")

# length
vioplot(nchar(info$seq), nchar(info.crd$seq), col = "darkolivegreen4",
        names= c("ALL INDEL","CORRECTED"))
title("SEQ LENGTH")

# seq length <-> frequency
par(mfrow = c(2,1))
alt = split(info$altFreq, as.factor(nchar(info$seq)))
boxplot(alt, col = "light green", xlab = "SEQ Length", ylab = "Frequency", 
        main = "ALL")

alt.crd = split(info.crd$altFreq, as.factor(nchar(info.crd$seq)))
boxplot(alt.crd, col = "light green", xlab = "SEQ Length", ylab = "Frequency", 
        main = "CORRECTED")
par(mfrow = c(1,1))

# SEQence Frequencys
seq.freq = summary(as.factor(paste0(info.crd$seq, " --> ", info.crd$alt)))
seq.freq = seq.freq[order(seq.freq, decreasing=T)]
old.par = par(mar = c(8,4,8,4))
barplot(seq.freq, las = 2, cex.names= 0.5, main = "CORRECTED", col= "black")

all.freq = summary(as.factor(paste0(info$seq, " --> ", info$alt)))
all.freq = all.freq[order(all.freq, decreasing=T)]
barplot(all.freq, las = 2, cex.names= 0.3, main = "ALL")
par(old.par)
dev.off()
