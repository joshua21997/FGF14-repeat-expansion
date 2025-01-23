args = commandArgs(trailingOnly=TRUE)
require("reshape2")
require("grid")
require(plyr)
library(reshape2)
library(grid)
library(brew)
library(cowplot)
library(dbplyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)


#Input
GAA <- read.table("GAA_summary_XXX.summary")
colnames(GAA)<- c("#line","motif","seq", "start", "end", "strand", "seqLen", "querybp", "mRatio", "m", "mm", "i", "d")


#Determination of the repeat length of the long allele
GAA <- filter(GAA, querybp > 600)
hist(GAA$querybp)
median(GAA$querybp)


#Interruption detection
filterplus <- filter(GAA, strand == "+")
filterplus <- filter(filterplus, querybp > 600)

tab <- data.frame(Position = c(1), Coverage = c(1))
test <- subset(filterplus, querybp > 0)
for(i in 1:1000) {
  VAL <- i
  test <- filter(test, querybp > VAL)
  A <- nrow(test)
  tab[i,1] <- i
  tab[i,2] <- A;
}

D <- read.table("del_+_XXX.tsv", header = FALSE, sep = "\t", col.names = paste0("V",seq_len(200)), fill = TRUE)
I <- read.table("ins_+_XXX.tsv", header = FALSE, sep = "\t", col.names = paste0("V",seq_len(200)), fill = TRUE)
A <- read.table("mm_a_+_XXX.tsv", header = FALSE, sep = "\t", col.names = paste0("V",seq_len(200)), fill = TRUE)
C <- read.table("mm_c_+_XXX.tsv", header = FALSE, sep = "\t", col.names = paste0("V",seq_len(200)), fill = TRUE)
G <- read.table("mm_g_+_XXX.tsv", header = FALSE, sep = "\t", col.names = paste0("V",seq_len(200)), fill = TRUE)
T <- read.table("mm_t_+_XXX.tsv", header = FALSE, sep = "\t", col.names = paste0("V",seq_len(200)), fill = TRUE)

names(D) <- c("seq", "start", "end", 1:197)
names(I) <- c("seq", "start", "end", 1:197)
names(A) <- c("seq", "start", "end", 1:197)
names(C) <- c("seq", "start", "end", 1:197)
names(G) <- c("seq", "start", "end", 1:197)
names(T) <- c("seq", "start", "end", 1:197)

DEL <- semi_join(D, filterplus, by = "seq")
INS <- semi_join(I, filterplus, by = "seq")
A <- semi_join(A, filterplus, by = "seq")
C <- semi_join(C, filterplus, by = "seq")
G <- semi_join(G, filterplus, by = "seq")
T <- semi_join(T, filterplus, by = "seq")

del <- D[,4:ncol(D)]
ins <- I[,4:ncol(I)]
a <- A[,4:ncol(A)]
c <- C[,4:ncol(C)]
g <- G[,4:ncol(G)]
t <- T[,4:ncol(T)]

bins <- seq(0,1000,by=1)
x <- seq(0, 999,by=1)

delall <- data.frame(V4 = c(t(del)))
delall <- delall[!is.na(delall)]
insall <- data.frame(V4 = c(t(ins)))
insall <- insall[!is.na(insall)]
Aall <- data.frame(V4 = c(t(a)))
Aall <- Aall[!is.na(Aall)]
Call <- data.frame(V4 = c(t(c)))
Call <- Call[!is.na(Call)]
Gall <- data.frame(V4 = c(t(g)))
Gall <- Gall[!is.na(Gall)]
Tall <- data.frame(V4 = c(t(t)))
Tall <- Tall[!is.na(Tall)]

scores_del <- cut(delall, bins)
freq_del <- transform(table(scores_del))
freq_del$xa <- x 
scores_ins <- cut(insall, bins)
freq_ins <- transform(table(scores_ins))
scores_a <- cut(Aall, bins)
freq_a <- transform(table(scores_a))
scores_c <- cut(Call, bins)
freq_c <- transform(table(scores_c))
scores_g <- cut(Gall, bins)
freq_g <- transform(table(scores_g))
scores_t <- cut(Tall, bins)
freq_t <- transform(table(scores_t))

freq_del$INS <- freq_ins$Freq
freq_del$A <- freq_a$Freq
freq_del$C <- freq_c$Freq
freq_del$G <- freq_g$Freq
freq_del$T <- freq_t$Freq
freq_del <- rename(freq_del,c("DEL" = "Freq"))

freq_del$COV <- tab$Coverage
freq_all <- melt(freq_del[,-1], id="xa")

ggplot(data=freq_all, aes(x=xa, y=value, group=variable, colour=variable)) +
  geom_line(linewidth=0.7) + theme_minimal(base_size = 25) + 
  theme(axis.text.x = element_text(angle = 90)) +
  theme(title = element_text(size = 20), axis.text = element_text(size = 20), axis.text.x = element_text(angle=90), legend.text=element_text(size = 20)) +
  ggtitle("") + 
  xlab("Sequence Position") + ylab("Number of Reads") 




#minus
filterminus <- filter(GAA, strand == "-")
filterminus <- filter(filterminus, querybp > 600)

tab <- data.frame(Position = c(1), Coverage = c(1))
test <- subset(filterminus, querybp > 0)
for(i in 1:1000) {
  VAL <- i
  test <- filter(test, querybp > VAL)
  A <- nrow(test)
  tab[i,1] <- i
  tab[i,2] <- A;
}

D <- read.table("del_-_XXX.tsv", header = FALSE, sep = "\t", col.names = paste0("V",seq_len(200)), fill = TRUE)
I <- read.table("ins_-_XXX.tsv", header = FALSE, sep = "\t", col.names = paste0("V",seq_len(200)), fill = TRUE)
A <- read.table("mm_a_-_XXX.tsv", header = FALSE, sep = "\t", col.names = paste0("V",seq_len(200)), fill = TRUE)
C <- read.table("mm_c_-_XXX.tsv", header = FALSE, sep = "\t", col.names = paste0("V",seq_len(200)), fill = TRUE)
G <- read.table("mm_g_-_XXX.tsv", header = FALSE, sep = "\t", col.names = paste0("V",seq_len(200)), fill = TRUE)
T <- read.table("mm_t_-_XXX.tsv", header = FALSE, sep = "\t", col.names = paste0("V",seq_len(200)), fill = TRUE)

names(D) <- c("seq", "start", "end", 1:197)
names(I) <- c("seq", "start", "end", 1:197)
names(A) <- c("seq", "start", "end", 1:197)
names(C) <- c("seq", "start", "end", 1:197)
names(G) <- c("seq", "start", "end", 1:197)
names(T) <- c("seq", "start", "end", 1:197)

DEL <- semi_join(D, filterminus, by = "seq")
INS <- semi_join(I, filterminus, by = "seq")
A <- semi_join(A, filterminus, by = "seq")
C <- semi_join(C, filterminus, by = "seq")
G <- semi_join(G, filterminus, by = "seq")
T <- semi_join(T, filterminus, by = "seq")

del <- DEL[,4:ncol(DEL)]
ins <- INS[,4:ncol(INS)]

a <- A[,4:ncol(A)]
c <- C[,4:ncol(C)]
g <- G[,4:ncol(G)]
t <- T[,4:ncol(T)]

bins <- seq(0,1000,by=1)
x <- seq(0, 999, by=1)

delall <- data.frame(V4 = c(t(del)))
delall <- delall[!is.na(delall)]
insall <- data.frame(V4 = c(t(ins)))
insall <- insall[!is.na(insall)]

Aall <- data.frame(V4 = c(t(a)))
Aall <- Aall[!is.na(Aall)]
Call <- data.frame(V4 = c(t(c)))
Call <- Call[!is.na(Call)]
Gall <- data.frame(V4 = c(t(g)))
Gall <- Gall[!is.na(Gall)]
Tall <- data.frame(V4 = c(t(t)))
Tall <- Tall[!is.na(Tall)]

scores_del <- cut(delall, bins)
freq_del <- transform(table(scores_del))
freq_del$xa <- x 
scores_ins <- cut(insall, bins)
freq_ins <- transform(table(scores_ins))
scores_a <- cut(Aall, bins)
freq_a <- transform(table(scores_a))
scores_c <- cut(Call, bins)
freq_c <- transform(table(scores_c))
scores_g <- cut(Gall, bins)
freq_g <- transform(table(scores_g))
scores_t <- cut(Tall, bins)
freq_t <- transform(table(scores_t))

freq_del$INS <- freq_ins$Freq
freq_del$A <- freq_a$Freq
freq_del$C <- freq_c$Freq
freq_del$G <- freq_g$Freq
freq_del$T <- freq_t$Freq

freq_del <- rename(freq_del,c("DEL" = "Freq"))
freq_del$COV <- tab$Coverage
freq_all <- melt(freq_del[,-1], id="xa")

ggplot(data=freq_all, aes(x=xa, y=value, group=variable, colour=variable)) +
  geom_line(linewidth=0.7) + theme_minimal(base_size = 25) + 
  theme(axis.text.x = element_text(angle = 90)) +
  theme(title = element_text(size = 20), axis.text = element_text(size = 20), axis.text.x = element_text(angle=90), legend.text=element_text(size = 20)) +
  ggtitle("") + 
  xlab("Sequence Position") + ylab("Number of Reads") 


