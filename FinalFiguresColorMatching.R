#set working directory
setwd("~/04_Champion_Lab/02_N-terminal_Acetylation/PDIM/Msphere/FinalFigures")

#load packages
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(ggrepel)
library(stringr)

#Accessions of proteins of interest to label
labels <- c("MMAR_5450", "MMAR_5457", "MMAR_5448", "MMAR_5449", "MMAR_4166", "MMAR_2894", "MMAR_5440", "MMAR_5439", "MMAR_5453", "MMAR_5455")
lowers <- c("PPE68", "EspA", "EspE", "EspF")

#####functions

#imports one biorep in formated dataframe, with optional gene to normalize to (default is RpoB)
#filters out contaminants
importBioRep <- function(file, sheet, biorep, type, normGene = "MMAR_0995") {
  df <- read_xlsx(file,
                  sheet = sheet)
  #filter Contaminants
  df <- df[!grepl("CONTAM", df$Accession),]
  
  #convert data to long format
  df2 <- df %>% pivot_longer(which(grepl("Area", names(df))),
                             values_to = "Area",
                             names_to = "Inj")
  
  #parse out replicates, and sample types, accessions, genes
  df2$Inj <- str_remove(df2$Inj, " Area")
  df2$techrep <- as.numeric(str_extract(df2$Inj, "[0-9]$"))
  df2$strain <- str_remove_all(str_extract(df2$Inj, "\\_.*\\_"), "\\_")
  df2$biorep <- biorep
  df2$type <- type
  df2$Accession2 <- str_remove(str_extract(df2$Accession, "^[^\\|]*\\|"), "\\|")
  df2$gene <- str_remove(str_remove(str_extract(df2$Accession, "^[^\\|]*\\|[^\\|]*\\|"),
                                    "^[^\\|]*\\|"), "\\|")
  
  #pull out normalization gene
  NormFactors <- df2[df2$Accession2 == normGene,]
  df2$AreaNorm <- NA
  
  #normalize areas
  for(i in 1:nrow(df2)) {
    df2$AreaNorm[i] <- df2$Area[i] / NormFactors$Area[NormFactors$techrep == df2$techrep[i] &
                                                        NormFactors$strain == df2$strain[i]]
  }
  return(df2)
}


#creates volcano plot dataframe
plotVolcano <- function(type,
                        comparison,
                        sigValue = 0.05,
                        L2FCcutoff = 1, return = "plot") {
  subdf <- comb[comb$type == type &
                  (comb$strain == comparison | comb$strain == "WT"),]
  
  subdf <- subdf %>% group_by(Accession2, gene) %>%
    mutate(CountWT = length(Area[strain == "WT"]),
           CountTest = length(Area[strain != "WT"]),
           maxCount = max(CountWT, CountTest),
           minCount = min(CountWT, CountTest))
  
  subdf <- subdf[subdf$maxCount >= 3,]
  compromisedDF <- subdf[subdf$minCount < 2,]
  subdf <- subdf[subdf$minCount >= 2,]
  
  
  volcano9 <- subdf %>% group_by(Accession2, gene) %>%
    summarise(totalCount = n(),
              CountWT = length(Area[strain == "WT"]),
              CountTest = length(Area[strain != "WT"]),
              meanWT = mean(Area[strain == "WT"], na.rm = T),
              meanTest = mean(Area[strain != "WT"], na.rm = T),
              foldChange = meanTest / meanWT,
              L2FC = log(foldChange, base = 2),
              pval = t.test(Area[strain == "WT"],
                            Area[strain != "WT"])
              [["p.value"]],
              pvalNeg10 = -(log(pval, base = 10)),
              compromised = F)
  
  volcano9Compromised <- compromisedDF %>% group_by(Accession2, gene) %>%
    summarise(totalCount = n(),
              CountWT = length(Area[strain == "WT"]),
              CountTest = length(Area[strain != "WT"]),
              meanWT = mean(Area[strain == "WT"], na.rm = T),
              meanTest = mean(Area[strain != "WT"], na.rm = T),
              foldChange = meanTest / meanWT,
              L2FC = log(foldChange, base = 2),
              pval = NA,
              pvalNeg10 = -(log(pval, base = 10)),
              compromised = T)
  
  volcano9Compromised$foldChange[volcano9Compromised$CountTest == 0] <- -Inf
  volcano9Compromised$foldChange[volcano9Compromised$CountWT == 0] <- Inf
  volcano9Compromised$L2FC[!is.finite(volcano9Compromised$foldChange)] <- volcano9Compromised$foldChange[!is.finite(volcano9Compromised$foldChange)]
  
  
  
  volcano9 <- volcano9[order(volcano9$pval),]
  volcano9$pvalRank <- 1:nrow(volcano9)
  volcano9$critValue <- (volcano9$pvalRank / nrow(volcano9[!is.na(volcano9$pval),])) * sigValue
  volcano9$pvalLessCritValue <- volcano9$pval <= volcano9$critValue
  
  sigDF <- volcano9[volcano9$pvalLessCritValue &
                      !is.na(volcano9$pvalLessCritValue),]
  
  if (nrow(sigDF) == 0) {
    sigValue <- min(volcano9$pval, na.rm = T) * 0.9
  } else {
    sigValue <- max(sigDF$pval, na.rm = T)
  }
  
  volcano9$Up <- (volcano9$L2FC > 1 & volcano9$pval <= sigValue)
  volcano9$Down <- (volcano9$L2FC < -1 & volcano9$pval <= sigValue)
  
  volcano9$plotCategory[volcano9$L2FC > L2FCcutoff & volcano9$pval <= sigValue] <- "UP"
  volcano9$plotCategory[volcano9$L2FC < -(L2FCcutoff) & volcano9$pval <= sigValue] <- "DOWN"
  
  if (nrow(volcano9[!is.na(volcano9$plotCategory),]) < 1) {
    volcano9$plotCategory <- "None"
  }
  
  
  volcano9$type <- type
  volcano9$comparison <- comparison
  
  volcano9Compromised$type <- type
  volcano9Compromised$comparison <- comparison
  
  xlim <- max(max(abs(volcano9$L2FC), na.rm = T), max(abs(volcano9Compromised$L2FC[is.finite(volcano9Compromised$L2FC)]), na.rm = T))
  ylim <- max(volcano9$pvalNeg10)
  plot <- ggplot() +
    geom_point(data = volcano9, aes(x = L2FC, y = pvalNeg10,
                                    color = plotCategory)) +
    geom_hline(yintercept = -(log(sigValue, base = 10)),
               linetype = 2) +
    geom_vline(xintercept = c(L2FCcutoff, - L2FCcutoff),
               linetype = 2) +
    geom_point(data = volcano9Compromised[is.finite(volcano9Compromised$L2FC),],
               aes(x = L2FC),
               y = 0) +
    geom_point(data = volcano9Compromised[is.infinite(volcano9Compromised$L2FC) & volcano9Compromised$L2FC > 0,],
               y = runif(nrow(volcano9Compromised[is.infinite(volcano9Compromised$L2FC) & volcano9Compromised$L2FC > 0,]), 0, ylim),
               x = xlim, color = "blue") +
    geom_point(data = volcano9Compromised[is.infinite(volcano9Compromised$L2FC) & volcano9Compromised$L2FC < 0,],
               y = runif(nrow(volcano9Compromised[is.infinite(volcano9Compromised$L2FC) & volcano9Compromised$L2FC < 0,]), 0, ylim),
               x = -xlim, color = "red") +
    coord_cartesian(xlim = c(-xlim, xlim)) +
    labs(x = "Fold Change (log2)",
         y = "Significance (-10log(p))",
         title = paste0(comparison, " vs WT BH corrected ", type)) +
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "none") +
    scale_color_manual(values = list("DOWN" = "#f8766d", "UP" = "#00BFC4", "None" = "grey45"))
  
  
  if (return != "plot") {
    return(list(dfvolcano = volcano9, dfcomp = volcano9Compromised, sigValue = sigValue))
  } else {
    return(plot)
  }
}


#create XY plot dataframe
returnScatterDF <- function(type, yaxis) {
  filtered <- comb[comb$type == type &
                     (comb$strain == "WT" | comb$strain == yaxis),]
  
  
  filtered2 <- filtered %>% group_by(Accession2, gene) %>%
    summarise(totalCount = n(),
              CountX = length(Area[strain == "WT"]),
              CountY = length(Area[strain == yaxis]),
              meanX = mean(Area[strain == "WT"], na.rm = T),
              meanY = mean(Area[strain == yaxis], na.rm = T))
  
  filtered2$type = type
  filtered2$comparison = yaxis
  
  return(filtered2)
}







#import bioreps
P1 <- importBioRep("PDIM_RAWs Processed.xlsx",
                   sheet = "P2 Raw",
                   biorep = 1,
                   type = "Cell Associated")

P2 <- importBioRep("PDIM_RAWs Processed.xlsx",
                   sheet = "P3 Raw",
                   biorep = 2,
                   type = "Cell Associated")
P3 <- importBioRep("PDIM_RAWs Processed.xlsx",
                   sheet = "P4 Raw",
                   biorep = 3,
                   type = "Cell Associated")


S1 <- importBioRep("PDIM_RAWs Processed.xlsx",
                   sheet = "S2 Raw",
                   biorep = 1,
                   type = "Secreted")

S2 <- importBioRep("PDIM_RAWs Processed.xlsx",
                   sheet = "S3 Raw",
                   biorep = 2,
                   type = "Secreted")
S3 <- importBioRep("PDIM_RAWs Processed.xlsx",
                   sheet = "S4 Raw",
                   biorep = 3,
                   type = "Secreted")





master <- bind_rows(P1,P2,P3,S1,S2,S3)



#strains to pull for figures 1
interest <- c("WT", "dmas", "dmasComp", "dCb")
comb <- master[master$strain %in% interest,]
#remove rows with no area
comb <- comb[comb$Area > 0,]



#calculate significance and fold change for each protein pairwise with wildtype
CAdmasComp <- plotVolcano("Cell Associated", "dmasComp", return = "df")
CAdCb <- plotVolcano("Cell Associated", "dCb", return = "df")
CAdmas <- plotVolcano("Cell Associated", "dmas", return = "df")
SdmasComp <- plotVolcano("Secreted", "dmasComp", return = "df")
SdCb <- plotVolcano("Secreted", "dCb", return = "df")
Sdmas <- plotVolcano("Secreted", "dmas", return = "df")


#create dataframes
volcanoComb <- bind_rows(CAdmasComp[[1]],
                         CAdCb[[1]],
                         CAdmas[[1]],
                         SdmasComp[[1]],
                         SdCb[[1]],
                         Sdmas[[1]])

volcanoCombComp <- bind_rows(CAdmasComp[[2]],
                             CAdCb[[2]],
                             CAdmas[[2]],
                             SdmasComp[[2]],
                             SdCb[[2]],
                             Sdmas[[2]])

sigs <- data.frame(sigValue = c(CAdmasComp[[3]],
                                CAdCb[[3]],
                                CAdmas[[3]],
                                SdmasComp[[3]],
                                SdCb[[3]],
                                Sdmas[[3]]),
                   type = c(rep("Cell Associated",3),
                            rep("Secreted",3)),
                   comparison = c("dmasComp",
                                  "dCb",
                                  "dmas",
                                  "dmasComp",
                                  "dCb",
                                  "dmas"))

#calculate plot limits
xlim <- max(max(abs(volcanoComb$L2FC), na.rm = T), max(abs(volcanoCombComp$L2FC[is.finite(volcanoCombComp$L2FC)]), na.rm = T))
ylim <- max(volcanoComb$pvalNeg10)

#fold change cutoff and dot size for plot
L2FCcutoff <- 1
dotSize <- 1.2


#calculate xy coords for compromised (infinite) values
volcanoCombComp$L2FC[is.infinite(volcanoCombComp$L2FC) & volcanoCombComp$L2FC > 0] <- xlim
volcanoCombComp$L2FC[is.infinite(volcanoCombComp$L2FC) & volcanoCombComp$L2FC < 0] <- -xlim
volcanoCombComp$pvalNeg10[is.infinite(volcanoCombComp$foldChange)] <- 
  runif(nrow(volcanoCombComp[is.infinite(volcanoCombComp$foldChange),]), 0, ylim)
volcanoCombComp$pvalNeg10[is.na(volcanoCombComp$pvalNeg10)] <- 0


#pull out proteins to label
Label <- bind_rows(volcanoComb, volcanoCombComp)
Label <- Label[Label$Accession2 %in% labels,]

#convert labels to prettier format
volcanoComb[volcanoComb == "dCb"] <- "∆eccCb1"
volcanoCombComp[volcanoCombComp == "dCb"] <- "∆eccCb1"
sigs[sigs == "dCb"] <- "∆eccCb1"
Label[Label == "dCb"] <- "∆eccCb1"

volcanoComb[volcanoComb == "dmas"] <- "∆mas"
volcanoCombComp[volcanoCombComp == "dmas"] <- "∆mas"
sigs[sigs == "dmas"] <- "∆mas"
Label[Label == "dmas"] <- "∆mas"

volcanoComb[volcanoComb == "dmasComp"] <- "∆mas/comp"
volcanoCombComp[volcanoCombComp == "dmasComp"] <- "∆mas/comp"
sigs[sigs == "dmasComp"] <- "∆mas/comp"
Label[Label == "dmasComp"] <- "∆mas/comp"

Label[Label == "MMAR_5448"] <- "PPE68"
Label[Label == "MMAR_2894"] <- "2894"
Label[Label == "MMAR_5440"] <- "espF"
Label[Label == "MMAR_5439"] <- "espE"
Label[Label == "MMAR_5453"] <- "espJ"
Label[Label == "MMAR_5455"] <- "espK"

Label$gene <- str_replace_all(Label$gene, "esp", "Esp")
Label$gene <- str_replace_all(Label$gene, "esx", "Esx")


#create plot
volcano1 <- ggplot() +
  geom_point(data = volcanoComb, aes(x = L2FC, y = pvalNeg10,
                                     color = plotCategory), size = dotSize) +
  geom_hline(data = sigs, aes(yintercept = -(log(sigValue, base = 10))),
             linetype = 2) +
  geom_vline(xintercept = c(L2FCcutoff, - L2FCcutoff),
             linetype = 2) +
  geom_hline(yintercept = 0) +
  geom_point(data = volcanoCombComp[is.finite(volcanoCombComp$foldChange),],
             aes(x = L2FC, y = pvalNeg10),
             y = 0, size = dotSize) +
  geom_point(data = volcanoCombComp[is.infinite(volcanoCombComp$foldChange) & volcanoCombComp$foldChange > 0,],
             aes(x = L2FC, y = pvalNeg10), color = "blue", size = dotSize) +
  geom_point(data = volcanoCombComp[is.infinite(volcanoCombComp$foldChange) & volcanoCombComp$foldChange < 0,],
             aes(x = L2FC, y = pvalNeg10), color = "red", size = dotSize) +
  coord_cartesian(xlim = c(-xlim, xlim)) +
  labs(x = "Fold Change (log2)",
       y = "Significance (-10log(p))") +
  theme_bw(base_size = 25)+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "white"),
        strip.text.y = element_text(face = "italic")) +
  scale_color_manual(values = list("DOWN" = "#f8766d", "UP" = "#00BFC4", "None" = "grey45")) +
  facet_grid(comparison~type) +
  geom_label_repel(data = Label, aes(x = L2FC, y = pvalNeg10, label = gene), nudge_x = 0,
                   nudge_y =3, size = 4,fill = alpha(c("white"),0.8))
volcano1
# 
# ggsave("plots/Volcano1-noNorm.png", height = 13, width = 11)
# write.csv(bind_rows(volcanoComb, volcanoCombComp),
#           "csvs/Volcano1-noNorm.csv", quote = F, row.names = F)
# 

#create dataframes for scatter
SdCb <- returnScatterDF("Secreted", "dCb")
Sdmas <-returnScatterDF("Secreted", "dmas")
SdmasComp <- returnScatterDF("Secreted", "dmasComp")

CAdCb <- returnScatterDF("Cell Associated", "dCb")
CAdmas <-returnScatterDF("Cell Associated", "dmas")
CAdmasComp <- returnScatterDF("Cell Associated", "dmasComp")

#combined into one dataframe
combscatter <- bind_rows(SdCb, Sdmas, SdmasComp, CAdCb, CAdmas, CAdmasComp)

#calculate log areas
combscatter$log10X <- log(combscatter$meanX, base = 10)
combscatter$log10Y <- log(combscatter$meanY, base = 10)

#put missing values in at zero to show on plot
combscatter$log10X[is.na(combscatter$meanX)] <- 0
combscatter$log10Y[is.na(combscatter$meanY)] <- 0

#pull out proteins of interest
Label <- combscatter[combscatter$Accession2 %in% labels,]


#convert to prettier labels
combscatter[combscatter == "dCb"] <- "∆eccCb1"
Label[Label == "dCb"] <- "∆eccCb1"

combscatter[combscatter == "dmas"] <- "∆mas"
Label[Label == "dmas"] <- "∆mas"

combscatter[combscatter == "dmasComp"] <- "∆mas/comp"
Label[Label == "dmasComp"] <- "∆mas/comp"

Label[Label == "MMAR_5448"] <- "PPE68"
Label[Label == "MMAR_2894"] <- "2894"
Label[Label == "MMAR_5440"] <- "espF"
Label[Label == "MMAR_5439"] <- "espE"
Label[Label == "MMAR_5453"] <- "espJ"
Label[Label == "MMAR_5455"] <- "espK"

Label$gene <- str_replace_all(Label$gene, "esp", "Esp")
Label$gene <- str_replace_all(Label$gene, "esx", "Esx")


group1 <- c("EsxA", "EsxB", "2894", "PPE68")
group2 <- c("EspJ", "EspB", "EspK")
group3 <- c("EspE", "EspF")
Other <- c("EspA")

Label$group[Label$gene %in% group1] <- "Group 1"
Label$group[Label$gene %in% group2] <- "Group 2"
Label$group[Label$gene %in% group3] <- "Group 3"
Label$group[Label$gene %in% Other] <- "other"


GroupColors <- c("#f00a00", "#399092", "#6c33a1", "grey30")
names(GroupColors) <- c("Group 1", "Group 2", "Group 3", "other")


XY1 <- ggplot() +
  geom_point(data = combscatter,
             aes(x = log10X, y = log10Y), size = 0.8, color = "grey60", shape = 21, fill = "grey60") +
  theme_bw(base_size = 25) +
  labs(x = "Wild Type (log 10 transformed)", y = "Experimental (log 10 transformed)") +
  facet_grid(comparison~type) +
  theme(strip.background = element_rect(fill = "white"),
        strip.text.y = element_text(face = "italic"),
        legend.position = "none") +
  geom_label_repel(data = Label[Label$gene %in% lowers &
                                  Label$type == "Secreted",], aes(x = log10X, y = log10Y, label = gene,
                                                                  color = group), nudge_x = 0,
                   nudge_y =-2, size = 4, fill = alpha(c("white"), 0.7),
                   segment.color = "grey30", seed = 1234, label.size = 1) +
  geom_label_repel(data = Label[!(Label$gene %in% lowers) &
                                  Label$type == "Secreted",], aes(x = log10X, y = log10Y, label = gene,
                                                                  color = group), nudge_x = -3,
                   nudge_y =1, size = 4,fill = alpha(c("white"),0.7),
                   segment.color = "grey30", seed = 1234, label.size = 1) +
  geom_label_repel(data = Label[Label$gene %in% lowers &
                                  Label$type == "Cell Associated",], aes(x = log10X, y = log10Y, label = gene,
                                                                         color = group), nudge_x = 0,
                   nudge_y =-2, size = 4, fill = alpha(c("white"), 0.7),
                   segment.color = "grey30", seed = 1234, label.size = 1) +
  geom_label_repel(data = Label[!(Label$gene %in% lowers)&
                                  Label$type == "Cell Associated",], aes(x = log10X, y = log10Y, label = gene,
                                                                         color = group), nudge_x = -2,
                   nudge_y =1, size = 4,fill = alpha(c("white"),0.7),
                   segment.color = "grey30", seed = 1234, label.size = 1) +
  geom_label_repel(data = Label[Label$gene %in% lowers &
                                  Label$type == "Secreted",], aes(x = log10X, y = log10Y, label = gene), nudge_x = 0,
                   nudge_y =-2, size = 4, fill = alpha(c("white"), 0.7),
                   segment.color = "grey30", color = "black", seed = 1234, label.size = NA) +
  geom_label_repel(data = Label[!(Label$gene %in% lowers) &
                                  Label$type == "Secreted",], aes(x = log10X, y = log10Y, label = gene), nudge_x = -3,
                   nudge_y =1, size = 4,fill = alpha(c("white"),0.7),
                   segment.color = "grey30", color = "black", seed = 1234, label.size = NA) +
  geom_label_repel(data = Label[Label$gene %in% lowers &
                                  Label$type == "Cell Associated",], aes(x = log10X, y = log10Y, label = gene), nudge_x = 0,
                   nudge_y =-2, size = 4, fill = alpha(c("white"), 0.7),
                   segment.color = "grey30", color = "black", seed = 1234, label.size = NA) +
  geom_label_repel(data = Label[!(Label$gene %in% lowers)&
                                  Label$type == "Cell Associated",], aes(x = log10X, y = log10Y, label = gene), nudge_x = -2,
                   nudge_y =1, size = 4,fill = alpha(c("white"),0.7),
                   segment.color = "grey30", color = "black", seed = 1234, label.size = NA) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "blue") +
  geom_point(data = Label,
             aes(x = log10X, y = log10Y,
                 color = group), size = 3, fill = "black", shape = 21, stroke = 2) +
  scale_color_manual(values = GroupColors)
XY1

ggsave("plots/Fancy/XY1-noNormMatchingcolors2.png", height = 11, width = 11)

# 
# ggsave("plots/XY1-noNorm.png", height = 11, width = 11)
# write.csv(combscatter,
#           "csvs/XY1-noNorm.csv", quote = F, row.names = F)






#SECOND FIGURES


#strains to pull for figures 2
interest2 <- c("WT", "ddrr", "dppsC", "dmas")
comb <- master[master$strain %in% interest2,]
#remove rows with no area
comb <- comb[comb$Area > 0,]



#calculate significance and fold change for each protein pairwise with wildtype
CAddrr <- plotVolcano("Cell Associated", "ddrr", return = "df")
CAdppsC <- plotVolcano("Cell Associated", "dppsC", return = "df")
CAdmas <- plotVolcano("Cell Associated", "dmas", return = "df")
Sddrr <- plotVolcano("Secreted", "ddrr", return = "df")
SdppsC <- plotVolcano("Secreted", "dppsC", return = "df")
Sdmas <- plotVolcano("Secreted", "dmas", return = "df")


#create dataframes
volcanoComb <- bind_rows(CAddrr[[1]],
                         CAdppsC[[1]],
                         CAdmas[[1]],
                         Sddrr[[1]],
                         SdppsC[[1]],
                         Sdmas[[1]])

volcanoCombComp <- bind_rows(CAddrr[[2]],
                             CAdppsC[[2]],
                             CAdmas[[2]],
                             Sddrr[[2]],
                             SdppsC[[2]],
                             Sdmas[[2]])

sigs <- data.frame(sigValue = c(CAddrr[[3]],
                                CAdppsC[[3]],
                                CAdmas[[3]],
                                Sddrr[[3]],
                                SdppsC[[3]],
                                Sdmas[[3]]),
                   type = c(rep("Cell Associated",3),
                            rep("Secreted",3)),
                   comparison = c("ddrr",
                                  "dppsC",
                                  "dmas",
                                  "ddrr",
                                  "dppsC",
                                  "dmas"))

#calculate plot limits
xlim <- max(max(abs(volcanoComb$L2FC), na.rm = T), max(abs(volcanoCombComp$L2FC[is.finite(volcanoCombComp$L2FC)]), na.rm = T))
ylim <- max(volcanoComb$pvalNeg10)

#fold change cutoff and dot size for plot
L2FCcutoff <- 1
dotSize <- 1.2


#calculate xy coords for compromised (infinite) values
volcanoCombComp$L2FC[is.infinite(volcanoCombComp$L2FC) & volcanoCombComp$L2FC > 0] <- xlim
volcanoCombComp$L2FC[is.infinite(volcanoCombComp$L2FC) & volcanoCombComp$L2FC < 0] <- -xlim
volcanoCombComp$pvalNeg10[is.infinite(volcanoCombComp$foldChange)] <- 
  runif(nrow(volcanoCombComp[is.infinite(volcanoCombComp$foldChange),]), 0, ylim)
volcanoCombComp$pvalNeg10[is.na(volcanoCombComp$pvalNeg10)] <- 0


#pull out proteins to label
Label <- bind_rows(volcanoComb, volcanoCombComp)
Label <- Label[Label$Accession2 %in% labels,]

#convert labels to prettier format
volcanoComb[volcanoComb == "ddrr"] <- "∆drr"
volcanoCombComp[volcanoCombComp == "ddrr"] <- "∆drr"
sigs[sigs == "ddrr"] <- "∆drr"
Label[Label == "ddrr"] <- "∆drr"

volcanoComb[volcanoComb == "dmas"] <- "∆mas"
volcanoCombComp[volcanoCombComp == "dmas"] <- "∆mas"
sigs[sigs == "dmas"] <- "∆mas"
Label[Label == "dmas"] <- "∆mas"

volcanoComb[volcanoComb == "dppsC"] <- "∆ppsC"
volcanoCombComp[volcanoCombComp == "dppsC"] <- "∆ppsC"
sigs[sigs == "dppsC"] <- "∆ppsC"
Label[Label == "dppsC"] <- "∆ppsC"


volcanoComb$comparison <- factor(volcanoComb$comparison, levels = c("∆mas", "∆drr", "∆ppsC"))
volcanoCombComp$comparison <- factor(volcanoCombComp$comparison, levels = c("∆mas", "∆drr", "∆ppsC"))
Label$comparison <- factor(Label$comparison, levels = c("∆mas", "∆drr", "∆ppsC"))
sigs$comparison <- factor(sigs$comparison, levels = c("∆mas", "∆drr", "∆ppsC"))


Label[Label == "MMAR_5448"] <- "PPE68"
Label[Label == "MMAR_2894"] <- "2894"
Label[Label == "MMAR_5440"] <- "espF"
Label[Label == "MMAR_5439"] <- "espE"
Label[Label == "MMAR_5453"] <- "espJ"
Label[Label == "MMAR_5455"] <- "espK"


Label$gene <- str_replace_all(Label$gene, "esp", "Esp")
Label$gene <- str_replace_all(Label$gene, "esx", "Esx")

#create plot
volcano2 <- ggplot() +
  geom_point(data = volcanoComb, aes(x = L2FC, y = pvalNeg10,
                                     color = plotCategory), size = dotSize) +
  geom_hline(data = sigs, aes(yintercept = -(log(sigValue, base = 10))),
             linetype = 2) +
  geom_vline(xintercept = c(L2FCcutoff, - L2FCcutoff),
             linetype = 2) +
  geom_hline(yintercept = 0) +
  geom_point(data = volcanoCombComp[is.finite(volcanoCombComp$foldChange),],
             aes(x = L2FC, y = pvalNeg10),
             y = 0, size = dotSize) +
  geom_point(data = volcanoCombComp[is.infinite(volcanoCombComp$foldChange) & volcanoCombComp$foldChange > 0,],
             aes(x = L2FC, y = pvalNeg10), color = "blue", size = dotSize) +
  geom_point(data = volcanoCombComp[is.infinite(volcanoCombComp$foldChange) & volcanoCombComp$foldChange < 0,],
             aes(x = L2FC, y = pvalNeg10), color = "red", size = dotSize) +
  coord_cartesian(xlim = c(-xlim, xlim)) +
  labs(x = "Fold Change (log2)",
       y = "Significance (-10log(p))") +
  theme_bw(base_size = 25)+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "white"),
        strip.text.y = element_text(face = "italic")) +
  scale_color_manual(values = list("DOWN" = "#f8766d", "UP" = "#00BFC4", "None" = "grey45")) +
  facet_grid(comparison~type) +
  geom_label_repel(data = Label, aes(x = L2FC, y = pvalNeg10, label = gene), nudge_x = 0,
                   nudge_y =3, size = 4,fill = alpha(c("white"),0.8))
volcano2

# ggsave("plots/Volcano2-noNorm.png", height = 13, width = 11)
# write.csv(bind_rows(volcanoComb, volcanoCombComp),
#           "csvs/Volcano2-noNorm.csv", quote = F, row.names = F)

#create dataframes for scatter
Sddrr <- returnScatterDF("Secreted", "ddrr")
Sdmas <-returnScatterDF("Secreted", "dmas")
SdppsC <- returnScatterDF("Secreted", "dppsC")

CAddrr <- returnScatterDF("Cell Associated", "ddrr")
CAdmas <-returnScatterDF("Cell Associated", "dmas")
CAdppsC <- returnScatterDF("Cell Associated", "dppsC")

#combined into one dataframe
combscatter <- bind_rows(Sddrr, Sdmas, SdppsC, CAddrr, CAdmas, CAdppsC)

#calculate log areas
combscatter$log10X <- log(combscatter$meanX, base = 10)
combscatter$log10Y <- log(combscatter$meanY, base = 10)

#put missing values in at zero to show on plot
combscatter$log10X[is.na(combscatter$meanX)] <- 0
combscatter$log10Y[is.na(combscatter$meanY)] <- 0

#pull out proteins of interest
Label <- combscatter[combscatter$Accession2 %in% labels,]


#convert to prettier labels
combscatter[combscatter == "ddrr"] <- "∆drr"
Label[Label == "ddrr"] <- "∆drr"

combscatter[combscatter == "dmas"] <- "∆mas"
Label[Label == "dmas"] <- "∆mas"

combscatter[combscatter == "dppsC"] <- "∆ppsC"
Label[Label == "dppsC"] <- "∆ppsC"

combscatter$comparison <- factor(combscatter$comparison, levels = c("∆mas", "∆drr", "∆ppsC"))
Label$comparison <- factor(Label$comparison, levels = c("∆mas", "∆drr", "∆ppsC"))


Label[Label == "MMAR_5448"] <- "PPE68"
Label[Label == "MMAR_2894"] <- "2894"
Label[Label == "MMAR_5440"] <- "espF"
Label[Label == "MMAR_5439"] <- "espE"
Label[Label == "MMAR_5453"] <- "espJ"
Label[Label == "MMAR_5455"] <- "espK"


Label$gene <- str_replace_all(Label$gene, "esp", "Esp")
Label$gene <- str_replace_all(Label$gene, "esx", "Esx")


group1 <- c("EsxA", "EsxB", "2894", "PPE68")
group2 <- c("EspJ", "EspB", "EspK")
group3 <- c("EspE", "EspF")
Other <- c("EspA")

Label$group[Label$gene %in% group1] <- "Group 1"
Label$group[Label$gene %in% group2] <- "Group 2"
Label$group[Label$gene %in% group3] <- "Group 3"
Label$group[Label$gene %in% Other] <- "other"


GroupColors <- c("#f00a00", "#399092", "#6c33a1", "grey30")
names(GroupColors) <- c("Group 1", "Group 2", "Group 3", "other")

#create plot
XY2 <- ggplot() +
  geom_point(data = combscatter,
             aes(x = log10X, y = log10Y), size = 0.8, color = "grey60", shape = 21, fill = "grey60") +
  theme_bw(base_size = 25) +
  labs(x = "Wild Type (log 10 transformed)", y = "Experimental (log 10 transformed)") +
  facet_grid(comparison~type) +
  theme(strip.background = element_rect(fill = "white"),
        strip.text.y = element_text(face = "italic"),
        legend.position = "none") +
  geom_label_repel(data = Label[Label$gene %in% lowers &
                                  Label$type == "Secreted",], aes(x = log10X, y = log10Y, label = gene,
                                                                  color = group), nudge_x = 0,
                   nudge_y =-2, size = 4, fill = alpha(c("white"), 0.7),
                   segment.color = "grey30", seed = 1234, label.size = 1) +
  geom_label_repel(data = Label[!(Label$gene %in% lowers) &
                                  Label$type == "Secreted",], aes(x = log10X, y = log10Y, label = gene,
                                                                  color = group), nudge_x = -3,
                   nudge_y =1, size = 4,fill = alpha(c("white"),0.7),
                   segment.color = "grey30", seed = 1234, label.size = 1) +
  geom_label_repel(data = Label[Label$gene %in% lowers &
                                  Label$type == "Cell Associated",], aes(x = log10X, y = log10Y, label = gene,
                                                                         color = group), nudge_x = 0,
                   nudge_y =-2, size = 4, fill = alpha(c("white"), 0.7),
                   segment.color = "grey30", seed = 1234, label.size = 1) +
  geom_label_repel(data = Label[!(Label$gene %in% lowers)&
                                  Label$type == "Cell Associated",], aes(x = log10X, y = log10Y, label = gene,
                                                                         color = group), nudge_x = -2,
                   nudge_y =1, size = 4,fill = alpha(c("white"),0.7),
                   segment.color = "grey30", seed = 1234, label.size = 1) +
  geom_label_repel(data = Label[Label$gene %in% lowers &
                                  Label$type == "Secreted",], aes(x = log10X, y = log10Y, label = gene), nudge_x = 0,
                   nudge_y =-2, size = 4, fill = alpha(c("white"), 0.7),
                   segment.color = "grey30", color = "black", seed = 1234, label.size = NA) +
  geom_label_repel(data = Label[!(Label$gene %in% lowers) &
                                  Label$type == "Secreted",], aes(x = log10X, y = log10Y, label = gene), nudge_x = -3,
                   nudge_y =1, size = 4,fill = alpha(c("white"),0.7),
                   segment.color = "grey30", color = "black", seed = 1234, label.size = NA) +
  geom_label_repel(data = Label[Label$gene %in% lowers &
                                  Label$type == "Cell Associated",], aes(x = log10X, y = log10Y, label = gene), nudge_x = 0,
                   nudge_y =-2, size = 4, fill = alpha(c("white"), 0.7),
                   segment.color = "grey30", color = "black", seed = 1234, label.size = NA) +
  geom_label_repel(data = Label[!(Label$gene %in% lowers)&
                                  Label$type == "Cell Associated",], aes(x = log10X, y = log10Y, label = gene), nudge_x = -2,
                   nudge_y =1, size = 4,fill = alpha(c("white"),0.7),
                   segment.color = "grey30", color = "black", seed = 1234, label.size = NA) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "blue") +
  geom_point(data = Label,
             aes(x = log10X, y = log10Y,
                 color = group), size = 3, fill = "black", shape = 21, stroke = 2) +
  scale_color_manual(values = GroupColors)
XY2

ggsave("plots/Fancy/XY2-noNormMatchingcolors2.png", height = 11, width = 11)

