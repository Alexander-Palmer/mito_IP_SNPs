###Determining SNP frequency of mouse brain tissues###
###For Ina Kirmes and Bartlett Lab, QBI###

##########
#Packages#
##########

library(ggplot2)
library(wesanderson)
library(xlsx)

############
#Data input#
############

features <- read.csv("mtDNA_features.csv")

S1825 <- read.csv("mutect2_1825_SNVs.csv")[,c(1,3:5,7,10,15)]
S1829 <- read.csv("mutect2_1829_SNVs.csv")[,c(1,3:5,7,10,15)]
S3444 <- read.csv("mutect2_3444_SNVs.csv")[,c(1,3:5,7,10,15)]
S3446 <- read.csv("mutect2_3446_SNVs.csv")[,c(1,3:5,7,10,15)]

S3450 <- read.csv("mutect2_3450_SNVs.csv")[,c(1,3:5,7,10,15)]
S3451 <- read.csv("mutect2_3451_SNVs.csv")[,c(1,3:5,7,10,15)]
S3453 <- read.csv("mutect2_3453_SNVs.csv")[,c(1,3:5,7,10,15)]
S3454 <- read.csv("mutect2_3454_SNVs.csv")[,c(1,3:5,7,10,15)]

S3458 <- read.csv("mutect2_3458_SNVs.csv")[,c(1,3:5,7,10,15)]
S3460 <- read.csv("mutect2_3460_SNVs.csv")[,c(1,3:5,7,10,15)]
S3461 <- read.csv("mutect2_3461_SNVs.csv")[,c(1,3:5,7,10,14)]
S3476 <- read.csv("mutect2_3476_SNVs.csv")[,c(1,3:5,7,10,14)]

################
#Data wrangling#
################

#Miscellaneous
#Annotate SNP by gene
data <- list(S1825, S1829, S3444, S3446, S3450, S3451, S3453, S3454, S3458,
             S3460, S3461, S3476)

data_out <- lapply(data, function(df) {
  df$gene <- ifelse(df[,1] >= 1 & df[,1] <= 68, print("tRNA-Phe"), 
                    (ifelse(df[,1] >= 70 & df[,1] <= 1024, print("rRNA-12S"),
                    ifelse(df[,1] >= 1025 & df[,1] <= 1093, print("tRNA-Val"),
                    ifelse(df[,1] >= 1094 & df[,1] <= 2675, print("rRNA-16S"), 
                    
                    ifelse(df[,1] >= 2676 & df[,1] <= 2750, print("tRNA-Leu"), 
                    ifelse(df[,1] >= 2751 & df[,1] <= 3707, print("ND1"), 
                    ifelse(df[,1] >= 3706 & df[,1] <= 3774, print("tRNA-Ile"), 
                    ifelse(df[,1] >= 3772 & df[,1] <= 3842, print("tRNA-Gln"),
                    ifelse(df[,1] >= 3845 & df[,1] <= 3913, print("tRNA-Met"), 
                    ifelse(df[,1] >= 3914 & df[,1] <= 4951, print("ND2"), 
                    ifelse(df[,1] >= 4950 & df[,1] <= 5016, print("tRNA-Trp"), 
                    ifelse(df[,1] >= 5018 & df[,1] <= 5086, print("tRNA-Ala"),
                                                                                            
                    ifelse(df[,1] >= 5089 & df[,1] <= 5159, print("tRNA-Asn"), 
                    ifelse(df[,1] >= 5160 & df[,1] <= 5191, print("L-strand ORI"), 
                    ifelse(df[,1] >= 5192 & df[,1] <= 5257, print("tRNA-Cys"), 
                    ifelse(df[,1] >= 5260 & df[,1] <= 5326, print("tRNA-Tyr"),
                    ifelse(df[,1] >= 5328 & df[,1] <= 6872, print("COXI"), 
                    ifelse(df[,1] >= 6870 & df[,1] <= 6938, print("tRNA-Ser"), 
                    ifelse(df[,1] >= 6942 & df[,1] <= 7011, print("tRNA-Asp"), 
                    ifelse(df[,1] >= 7013 & df[,1] <= 7696, print("COXII"),       
                            
                    ifelse(df[,1] >= 7700 & df[,1] <= 7764, print("tRNA-Lys"), 
                    ifelse(df[,1] >= 7766 & df[,1] <= 7969, print("ATP8"), 
                    ifelse(df[,1] >= 7927 & df[,1] <= 8607, print("ATP6"), 
                    ifelse(df[,1] >= 8607 & df[,1] <= 9390, print("COXIII"),
                    ifelse(df[,1] >= 9391 & df[,1] <= 9458, print("tRNA-Gly"), 
                    ifelse(df[,1] >= 9459 & df[,1] <= 9806, print("ND3"), 
                    ifelse(df[,1] >= 9808 & df[,1] <= 9875, print("tRNA-Arg"), 
                    ifelse(df[,1] >= 9821 & df[,1] <= 9830, print("Polymorphic region"),        
                            
                    ifelse(df[,1] >= 9877 & df[,1] <= 10173, print("ND4L"), 
                    ifelse(df[,1] >= 10167 & df[,1] <= 11544, print("ND4"), 
                    ifelse(df[,1] >= 11545 & df[,1] <= 11612, print("tRNA-His"), 
                    ifelse(df[,1] >= 11613 & df[,1] <= 11671, print("tRNA-Ser"),
                    ifelse(df[,1] >= 11672 & df[,1] <= 11742, print("tRNA-Leu"), 
                    ifelse(df[,1] >= 11742 & df[,1] <= 13565, print("ND5"), 
                    ifelse(df[,1] >= 13552 & df[,1] <= 14070, print("ND6"), 
                    ifelse(df[,1] >= 14071 & df[,1] <= 14139, print("tRNA-Glu"),        
                      
                    ifelse(df[,1] >= 14145 & df[,1] <= 15228, print("CYTB"), 
                    ifelse(df[,1] >= 15289 & df[,1] <= 15355, print("tRNA-Thr"), 
                    ifelse(df[,1] >= 15356 & df[,1] <= 15422, print("tRNA-Pro"),
                    ifelse(df[,1] >= 15423 & df[,1] <= 16299, print("D-loop"),
                    ifelse(df[,1] >= 15451 & df[,1] <= 15509, print("ETAS1"), 
                    ifelse(df[,1] >= 15515 & df[,1] <= 15558, print("ETAS2"), 
                    ifelse(df[,1] >= 16035 & df[,1] <= 16058, print("CSB1"), 
                    ifelse(df[,1] >= 16089 & df[,1] <= 16104, print("CSB2"),       
                    ifelse(df[,1] >= 16114 & df[,1] <= 16131, print("CSB3"), print("intergenic")       
                    ))))))))))))))))))))))))))))))))))))))))))))))
                    
  return(df)
})

S1825 <- as.data.frame(data_out[1])
S1829 <- as.data.frame(data_out[2])
S3444 <- as.data.frame(data_out[3])
S3446 <- as.data.frame(data_out[4])

S3450 <- as.data.frame(data_out[5])
S3451 <- as.data.frame(data_out[6])
S3453 <- as.data.frame(data_out[7])
S3454 <- as.data.frame(data_out[8])

S3458 <- as.data.frame(data_out[1])
S3460 <- as.data.frame(data_out[10])
S3461 <- as.data.frame(data_out[11])
S3476 <- as.data.frame(data_out[12])

#Introduce TLOD cutoff of >= 5
#In literature, TLOD cutoff >= 2, 
S1825 <- subset(S1825, TLOD >= 5)
S1829 <- subset(S1829, TLOD >= 5)
S3444 <- subset(S3444, TLOD >= 5)
S3446 <- subset(S3446, TLOD >= 5)

S3450 <- subset(S3450, TLOD >= 5)
S3451 <- subset(S3451, TLOD >= 5)
S3453 <- subset(S3453, TLOD >= 5)
S3454 <- subset(S3454, TLOD >= 5)

S3458 <- subset(S3458, TLOD >= 5)
S3460 <- subset(S3460, TLOD >= 5)
S3461 <- subset(S3461, TLOD >= 5)
S3476 <- subset(S3476, TLOD >= 5)

#Total SNP count
S1825_count <- nrow(S1825)
S1829_count <- nrow(S1829)
S3444_count <- nrow(S3444)
S3446_count <- nrow(S3446)

S3450_count <- nrow(S3450)
S3451_count <- nrow(S3451)
S3453_count <- nrow(S3453)
S3454_count <- nrow(S3454)

S3458_count <- nrow(S3458)
S3460_count <- nrow(S3460)
S3461_count <- nrow(S3461)
S3476_count <- nrow(S3476)

sample_count <- data.frame(
  "Sample"=c("S1825", "S1829", "S3444", "S3446", "S3450", "S3451", "S3453", 
             "S3454", "S3458", "S3460", "S3461", "S3476"),
  "Count"=c(S1825_count, S1829_count, S3444_count, S3446_count, S3450_count,
          S3451_count, S3453_count, S3454_count, S3458_count, S3460_count, 
          S3461_count, S3476_count))

ggplot(sample_count, aes(x=Sample, y=Count, fill=Sample)) + 
  geom_col(show.legend=FALSE) +
  scale_fill_manual(values = wes_palette(12, name="Darjeeling1", "continuous"))+
  xlab("Library") +
  ylab("Count (All SNPs)")

#Mapped bases number mtDNA (from 'qualimap_out.html')
#Used, as total mapped reads may map across different chromosomes
S1825_lib_size <- 172126609
S1829_lib_size <- 24114528
S3444_lib_size <- 68736322
S3446_lib_size <- 15683186
  
S3450_lib_size <- 64856438
S3451_lib_size <- 68749374
S3453_lib_size <- 112180561
S3454_lib_size <- 76371138
  
S3458_lib_size <- 452490814
S3460_lib_size <- 80416364
S3461_lib_size <- 84692966
S3476_lib_size <- 34089596

library_count <- data.frame(
  "Sample"=c("S1825", "S1829", "S3444", "S3446", "S3450", "S3451", "S3453", 
             "S3454", "S3458", "S3460", "S3461", "S3476"),
  "Count"=c(S1825_lib_size, S1829_lib_size, S3444_lib_size, S3446_lib_size, 
            S3450_lib_size, S3451_lib_size, S3453_lib_size, S3454_lib_size, 
            S3458_lib_size, S3460_lib_size, S3461_lib_size, S3476_lib_size))

ggplot(library_count, aes(x=Sample, y=Count, fill=Sample)) + 
  geom_col(show.legend=FALSE) +
  scale_fill_manual(values = wes_palette(12, name="Darjeeling1", "continuous"))+
  xlab("Library") +
  ylab("Mapped bases")

#Normalize total SNP count by library size for true SNP count/sample
S1825_norm <- S1825_count / S1825_lib_size
S1829_norm <- S1829_count / S1829_lib_size
S3444_norm <- S3444_count / S3444_lib_size
S3446_norm <- S3446_count / S3446_lib_size

S3450_norm <- S3450_count / S3450_lib_size
S3451_norm <- S3451_count / S3451_lib_size
S3453_norm <- S3453_count / S3453_lib_size
S3454_norm <- S3454_count / S3454_lib_size

S3458_norm <- S3458_count / S3458_lib_size
S3460_norm <- S3460_count / S3460_lib_size
S3461_norm <- S3461_count / S3461_lib_size
S3476_norm <- S3476_count / S3476_lib_size

norm_count <- data.frame(
  "Sample"=c("S1825", "S1829", "S3444", "S3446", "S3450", "S3451", "S3453", 
             "S3454", "S3458", "S3460", "S3461", "S3476"),
  "Count"=c(S1825_norm, S1829_norm, S3444_norm, S3446_norm, S3450_norm, 
            S3451_norm, S3453_norm, S3454_norm, S3458_norm, S3460_norm, 
            S3461_norm, S3476_norm))

ggplot(norm_count, aes(x=Sample, y=Count, fill=Sample)) + 
  geom_col(show.legend=FALSE) +
  scale_fill_manual(values = wes_palette(12, name="Darjeeling1", "continuous"))+
  xlab("Library") +
  ylab("SNPs/mapped bases") 

################################################################################
################################################################################
#Compare individual SNP percentages

#Create sum tables
SNP_counts <- lapply(data, function(df) {
  
  SNP_list <- as.data.frame(table(df[,3]))
  
  A_T <- sum(ifelse(SNP_list[,1]== "A/T", SNP_list[,2], 0))
  A_C <- sum(ifelse(SNP_list[,1]== "A/C", SNP_list[,2], 0))
  A_G <- sum(ifelse(SNP_list[,1]== "A/G", SNP_list[,2], 0))
  
  T_A <- sum(ifelse(SNP_list[,1]== "T/A", SNP_list[,2], 0))
  T_C <- sum(ifelse(SNP_list[,1]== "T/C", SNP_list[,2], 0))
  T_G <- sum(ifelse(SNP_list[,1]== "T/G", SNP_list[,2], 0))
  
  C_A <- sum(ifelse(SNP_list[,1]== "C/A", SNP_list[,2], 0))
  C_T <- sum(ifelse(SNP_list[,1]== "C/T", SNP_list[,2], 0))
  C_G <- sum(ifelse(SNP_list[,1]== "C/G", SNP_list[,2], 0))
  
  G_A <- sum(ifelse(SNP_list[,1]== "G/A", SNP_list[,2], 0))
  G_T <- sum(ifelse(SNP_list[,1]== "G/T", SNP_list[,2], 0))
  G_C <- sum(ifelse(SNP_list[,1]== "G/C", SNP_list[,2], 0))
  
  table <- as.data.frame(rbind(A_T, A_C, A_G, T_A, T_C, T_G, 
                               C_A, C_T, C_G, G_A, G_T, G_C))
  
  table$Allele <- rownames(table)
  table <- table[,c(2,1)]
  colnames(table) <- c("Allele", "Count")
  rownames(table) <- c(1:12)
  
  return(table)
})

#Create individual tables  
S1825_SNP <- as.data.frame(SNP_counts[1])
S1829_SNP <- as.data.frame(SNP_counts[2])
S3444_SNP <- as.data.frame(SNP_counts[3])
S3446_SNP <- as.data.frame(SNP_counts[4])

S3450_SNP <- as.data.frame(SNP_counts[5])
S3451_SNP <- as.data.frame(SNP_counts[6])
S3453_SNP <- as.data.frame(SNP_counts[7])
S3454_SNP <- as.data.frame(SNP_counts[8])

S3458_SNP <- as.data.frame(SNP_counts[9])
S3460_SNP <- as.data.frame(SNP_counts[10])
S3461_SNP <- as.data.frame(SNP_counts[11])
S3476_SNP <- as.data.frame(SNP_counts[12])

#Create combined table
SNP_list <- as.data.frame(cbind(S1825_SNP$Count, S1829_SNP$Count, 
                                S3444_SNP$Count, S3446_SNP$Count, 
                                S3450_SNP$Count, S3451_SNP$Count,
                                S3453_SNP$Count, S3454_SNP$Count,
                                S3458_SNP$Count, S3460_SNP$Count, 
                                S3461_SNP$Count, S3476_SNP$Count))

rownames(SNP_list) <- (c("A/T", "A/C", "A/G", "T/A", "T/C", "T/G",
                         "C/A", "C/T", "C/G", "G/A", "G/T", "G/C"))
colnames(SNP_list) <- (c("S1825", "S1829", "S3444", "S3446", "S3450", "S3451", 
                         "S3453", "S3454", "S3458", "S3460", "S3461", "S3476"))



SNP_list_combined <- as.data.frame(rowSums(SNP_list))
SNP_list_combined$Allele <- rownames(SNP_list_combined)
SNP_list_combined <- SNP_list_combined[,c(2,1)]
colnames(SNP_list_combined) <- c("Allele", "Count")

SNP_list$Allele <- rownames(SNP_list)
SNP_list <- SNP_list[,c(13,1:12)]

#Plot combined table
ggplot(SNP_list_combined, aes(x=Allele, y=Count, fill=Allele)) + 
  geom_col(show.legend=FALSE) +
  scale_fill_manual(values = wes_palette(12, name="Darjeeling1", "continuous"))+
  xlab("Allele") +
  ylab("Number of SNPs") 

#Plot individual tables 

ggplot(S3476_SNP, aes(x=Allele, y=Count, fill=Allele)) + 
  geom_col(show.legend=FALSE) +
  scale_fill_manual(values = wes_palette(12, name="Darjeeling1", "continuous"))+
  xlab("Allele") +
  ylab("Number of SNPs") 

#.....
  
################################################################################
################################################################################ 
#Analyse SNP position

#Create sum tables
data <- list(S1825, S1829, S3444, S3446, S3450, S3451, S3453, S3454, S3458,
             S3460, S3461, S3476)

gene_counts <- lapply(data, function(df) {
  
  gene_list <- as.data.frame(table(df[,8]))
  
  "tRNA_Phe" <- sum(ifelse(gene_list[,1]== "tRNA-Phe", gene_list[,2], 0))
  "rRNA_12S" <- sum(ifelse(gene_list[,1]== "rRNA-12S", gene_list[,2], 0))
  "tRNA_Val" <- sum(ifelse(gene_list[,1]== "tRNA-Val", gene_list[,2], 0))
  "rRNA_16S" <- sum(ifelse(gene_list[,1]== "rRNA-16S", gene_list[,2], 0))
  "tRNA_Leu" <- sum(ifelse(gene_list[,1]== "tRNA-Leu", gene_list[,2], 0))
  "ND1" <- sum(ifelse(gene_list[,1]== "ND1", gene_list[,2], 0))
  "tRNA_Ile" <- sum(ifelse(gene_list[,1]== "tRNA-Ile", gene_list[,2], 0))
  "tRNA_Gln" <- sum(ifelse(gene_list[,1]== "tRNA-Gln", gene_list[,2], 0))
  
  "tRNA_Met" <- sum(ifelse(gene_list[,1]== "tRNA-Met", gene_list[,2], 0))
  "ND2" <- sum(ifelse(gene_list[,1]== "ND2", gene_list[,2], 0))
  "tRNA_Trp" <- sum(ifelse(gene_list[,1]== "tRNA-Trp", gene_list[,2], 0))
  "tRNA_Ala" <- sum(ifelse(gene_list[,1]== "tRNA-Ala", gene_list[,2], 0))
  "tRNA_Asn" <- sum(ifelse(gene_list[,1]== "tRNA-Asn", gene_list[,2], 0))
  "L_strand_ORI" <- sum(ifelse(gene_list[,1]== "L-strand ORI", gene_list[,2], 0))
  "tRNA_Cys" <- sum(ifelse(gene_list[,1]== "tRNA-Cys", gene_list[,2], 0))
  "tRNA_Tyr" <- sum(ifelse(gene_list[,1]== "tRNA-Tyr", gene_list[,2], 0))
  
  "COXI" <- sum(ifelse(gene_list[,1]== "COXI", gene_list[,2], 0))
  "tRNA_Ser" <- sum(ifelse(gene_list[,1]== "tRNA-Ser", gene_list[,2], 0))
  "tRNA_Asp" <- sum(ifelse(gene_list[,1]== "tRNA-Asp", gene_list[,2], 0))
  "COXII" <- sum(ifelse(gene_list[,1]== "COXII", gene_list[,2], 0))
  "tRNA_Lys" <- sum(ifelse(gene_list[,1]== "tRNA-Lys", gene_list[,2], 0))
  "ATP8" <- sum(ifelse(gene_list[,1]== "ATP8", gene_list[,2], 0))
  "ATP6" <- sum(ifelse(gene_list[,1]== "ATP6", gene_list[,2], 0))
  "COXIII" <- sum(ifelse(gene_list[,1]== "COXIII", gene_list[,2], 0))
  
  "tRNA_Gly" <- sum(ifelse(gene_list[,1]== "tRNA-Gly", gene_list[,2], 0))
  "ND3" <- sum(ifelse(gene_list[,1]== "ND3", gene_list[,2], 0))
  "tRNA_Arg" <- sum(ifelse(gene_list[,1]== "tRNA-Arg", gene_list[,2], 0))
  "Polymorphic_region" <- sum(ifelse(gene_list[,1]== "Polymorphic region", gene_list[,2], 0))
  "ND4L" <- sum(ifelse(gene_list[,1]== "ND4L", gene_list[,2], 0))
  "ND4" <- sum(ifelse(gene_list[,1]== "ND4", gene_list[,2], 0))
  "tRNA_His" <- sum(ifelse(gene_list[,1]== "tRNA-His", gene_list[,2], 0))
  "tRNA_Ser" <- sum(ifelse(gene_list[,1]== "tRNA-Ser", gene_list[,2], 0))
  
  "tRNA_Leu" <- sum(ifelse(gene_list[,1]== "tRNA-Leu", gene_list[,2], 0))
  "ND5" <- sum(ifelse(gene_list[,1]== "ND5", gene_list[,2], 0))
  "ND6" <- sum(ifelse(gene_list[,1]== "ND6", gene_list[,2], 0))
  "tRNA_Glu" <- sum(ifelse(gene_list[,1]== "tRNA-Glu", gene_list[,2], 0))
  "CYTB" <- sum(ifelse(gene_list[,1]== "CYTB", gene_list[,2], 0))
  "tRNA_Thr" <- sum(ifelse(gene_list[,1]== "tRNA-Thr", gene_list[,2], 0))
  "tRNA_Pro" <- sum(ifelse(gene_list[,1]== "tRNA-Pro", gene_list[,2], 0))
  "D_Loop" <- sum(ifelse(gene_list[,1]== "D-loop", gene_list[,2], 0))

  
  "ETAS1" <- sum(ifelse(gene_list[,1]== "ETAS1", gene_list[,2], 0))
  
  "ETAS2" <- sum(ifelse(gene_list[,1]== "ETAS2", gene_list[,2], 0))
  "CSB1" <- sum(ifelse(gene_list[,1]== "CSB1", gene_list[,2], 0))
  "CSB2" <- sum(ifelse(gene_list[,1]== "CSB2", gene_list[,2], 0))
  "CSB3" <- sum(ifelse(gene_list[,1]== "CSB3", gene_list[,2], 0))
  "intergenic" <- sum(ifelse(gene_list[,1]== "intergenic", gene_list[,2], 0))
  
  table <- as.data.frame(rbind(tRNA_Phe, rRNA_12S, tRNA_Val, rRNA_16S, tRNA_Leu, 
                               ND1, tRNA_Ile, tRNA_Gln, tRNA_Met, ND2, tRNA_Trp,
                               tRNA_Ala, tRNA_Asn, L_strand_ORI, tRNA_Cys,
                               tRNA_Tyr, COXI, tRNA_Ser, tRNA_Asp, COXII,
                               tRNA_Lys, ATP8, ATP6, COXIII, tRNA_Gly, ND3,
                               tRNA_Arg, Polymorphic_region, ND4L, ND4,
                               tRNA_His, tRNA_Ser, tRNA_Leu, ND5, ND6, tRNA_Glu,
                               CYTB, tRNA_Thr, tRNA_Pro, D_Loop, ETAS1, ETAS2, 
                               CSB1, CSB2, CSB3, intergenic))
  
  table$Gene <- rownames(table)
  table <- table[,c(2,1)]
  colnames(table) <- c("Gene", "Count")
  rownames(table) <- c(1:46)
  
  return(table)
})

#Create individual tables  
S1825_gene <- as.data.frame(gene_counts[1])
S1829_gene <- as.data.frame(gene_counts[2])
S3444_gene <- as.data.frame(gene_counts[3])
S3446_gene <- as.data.frame(gene_counts[4])

S3450_gene <- as.data.frame(gene_counts[5])
S3451_gene <- as.data.frame(gene_counts[6])
S3453_gene <- as.data.frame(gene_counts[7])
S3454_gene <- as.data.frame(gene_counts[8])

S3458_gene <- as.data.frame(gene_counts[9])
S3460_gene <- as.data.frame(gene_counts[10])
S3461_gene <- as.data.frame(gene_counts[11])
S3476_gene <- as.data.frame(gene_counts[12])

#Plot individual tables 
ggplot(S3476_gene, aes(x=Gene, y=Count, fill=Gene)) + 
  geom_col(show.legend=FALSE) +
  scale_fill_manual(values = wes_palette(46, name="Darjeeling1", "continuous"))+
  xlab("Gene") +
  ylab("Number of SNPs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#.....

#Create combined table
gene_list <- as.data.frame(cbind(S1825_gene$Count, S1829_gene$Count, 
                                 S3444_gene$Count, S3446_gene$Count, 
                                 S3450_gene$Count, S3451_gene$Count,
                                 S3453_gene$Count, S3454_gene$Count,
                                 S3458_gene$Count, S3460_gene$Count, 
                                 S3461_gene$Count, S3476_gene$Count))

rownames(gene_list) <- S1825_gene$Gene
colnames(gene_list) <- (c("S1825", "S1829", "S3444", "S3446", "S3450", "S3451", 
                          "S3453", "S3454", "S3458", "S3460", "S3461", "S3476"))

gene_list_combined <- as.data.frame(rowSums(gene_list))
gene_list_combined$Gene <- rownames(gene_list_combined)
gene_list_combined <- gene_list_combined[,c(2,1)]
colnames(gene_list_combined) <- c("Gene", "Count")

gene_list$Gene <- rownames(gene_list)
gene_list <- gene_list[,c(13,1:12)]

#Plot combined table
ggplot(gene_list_combined, aes(x=Gene, y=Count, fill=Gene)) + 
  geom_col(show.legend=FALSE) +
  scale_fill_manual(values = wes_palette(46, name="Darjeeling1", "continuous"))+
  xlab("Gene") +
  ylab("Number of SNPs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

################################################################################

#Normalise counts to gene length (generates FPKM)
"tRNA_Phe_length" <- 68
"rRNA_12S_length" <- 955
"tRNA_Val_length" <- 69
"rRNA_16S_length" <- 1582
"tRNA_Leu_length" <- 75
"ND1_length" <- 957
"tRNA_Ile_length" <- 69
"tRNA_Gln_length" <- 71

"tRNA_Met_length" <- 69
"ND2_length" <- 1038
"tRNA_Trp_length" <- 67 
"tRNA_Ala_length" <- 69
"tRNA_Asn_length" <- 71
"L_strand_ORI_length" <- 32
"tRNA_Cys_length" <- 66
"tRNA_Tyr_length" <- 67

"COXI_length" <- 1545
"tRNA_Ser_length" <- 69
"tRNA_Asp_length" <- 70
"COXII_length" <- 684
"tRNA_Lys_length" <- 65
"ATP8_length" <- 204
"ATP6_length" <- 681
"COXIII_length" <- 784

"tRNA_Gly_length" <- 68
"ND3_length" <- 348
"tRNA_Arg_length" <- 68
"Polymorphic_region_length" <- 10
"ND4L_length" <- 297
"ND4_length" <- 1378
"tRNA_His_length" <- 68
"tRNA_Ser_length" <- 59

"tRNA_Leu_length" <- 71
"ND5_length" <- 1824
"ND6_length" <- 519
"tRNA_Glu_length" <- 69
"CYTB_length" <- 1144
"tRNA_Thr_length" <- 67
"tRNA_Pro_length" <- 67
"D_loop_length" <- 877

"ETAS1_length" <- 59
"ETAS2_length" <- 44
"CSB1_length" <- 24
"CSB2_length" <- 16
"CSB3_length" <- 18
"intergenic_length" <- 1

#Could have just put this in a separate graph... but oh well, too late now 
#Normalise individual libraries
gene_list_norm <- gene_list  

gene_list_norm[1,] <- gene_list[1,] / tRNA_Phe_length
gene_list_norm[2,] <- gene_list[2,] / rRNA_12S_length
gene_list_norm[3,] <- gene_list[3,] / tRNA_Val_length
gene_list_norm[4,] <- gene_list[4,] / rRNA_16S_length
gene_list_norm[5,] <- gene_list[5,] / tRNA_Leu_length
gene_list_norm[6,] <- gene_list[6,] / ND1_length
gene_list_norm[7,] <- gene_list[7,] / tRNA_Ile_length
gene_list_norm[8,] <- gene_list[8,] / tRNA_Gln_length

gene_list_norm[9,] <- gene_list[9,] / tRNA_Met_length
gene_list_norm[10,] <- gene_list[10,] / ND2_length
gene_list_norm[11,] <- gene_list[11,] / tRNA_Trp_length
gene_list_norm[12,] <- gene_list[12,] / tRNA_Ala_length
gene_list_norm[13,] <- gene_list[13,] / tRNA_Asn_length
gene_list_norm[14,] <- gene_list[14,] / L_strand_ORI_length
gene_list_norm[15,] <- gene_list[15,] / tRNA_Cys_length
gene_list_norm[16,] <- gene_list[16,] / tRNA_Tyr_length

gene_list_norm[17,] <- gene_list[17,] / COXI_length
gene_list_norm[18,] <- gene_list[18,] / tRNA_Ser_length
gene_list_norm[19,] <- gene_list[19,] / tRNA_Asp_length
gene_list_norm[20,] <- gene_list[20,] / COXII_length
gene_list_norm[21,] <- gene_list[21,] / tRNA_Lys_length
gene_list_norm[22,] <- gene_list[22,] / ATP8_length
gene_list_norm[23,] <- gene_list[23,] / ATP6_length
gene_list_norm[24,] <- gene_list[24,] / COXIII_length

gene_list_norm[25,] <- gene_list[25,] / tRNA_Gly_length
gene_list_norm[26,] <- gene_list[26,] / ND3_length
gene_list_norm[27,] <- gene_list[27,] / tRNA_Arg_length
gene_list_norm[28,] <- gene_list[28,] / Polymorphic_region_length
gene_list_norm[29,] <- gene_list[29,] / ND4L_length
gene_list_norm[30,] <- gene_list[30,] / ND4_length
gene_list_norm[31,] <- gene_list[31,] / tRNA_His_length
gene_list_norm[32,] <- gene_list[32,] / tRNA_Ser_length

gene_list_norm[33,] <- gene_list[33,] / tRNA_Leu_length
gene_list_norm[34,] <- gene_list[34,] / ND5_length
gene_list_norm[35,] <- gene_list[35,] / ND6_length
gene_list_norm[36,] <- gene_list[36,] / tRNA_Glu_length
gene_list_norm[37,] <- gene_list[37,] / CYTB_length
gene_list_norm[38,] <- gene_list[38,] / tRNA_Thr_length
gene_list_norm[39,] <- gene_list[39,] / tRNA_Pro_length
gene_list_norm[40,] <- gene_list[40,] / D_loop_length

gene_list_norm[41,] <- gene_list[41,] / ETAS1_length
gene_list_norm[42,] <- gene_list[42,] / ETAS2_length
gene_list_norm[43,] <- gene_list[43,] / CSB1_length
gene_list_norm[44,] <- gene_list[44,] / CSB2_length
gene_list_norm[45,] <- gene_list[45,] / CSB3_length
gene_list_norm[46,] <- gene_list[46,] / intergenic_length

gene_list_norm$Gene <- rownames(gene_list_norm)
gene_list_norm <- gene_list_norm[,c(13,1:12)]

S1825_gene_norm <- as.data.frame(gene_list_norm[,c(1,2)])
S1829_gene_norm <- as.data.frame(gene_list_norm[,c(1,3)])
S3444_gene_norm <- as.data.frame(gene_list_norm[,c(1,4)])
S3446_gene_norm <- as.data.frame(gene_list_norm[,c(1,5)])

S3450_gene_norm <- as.data.frame(gene_list_norm[,c(1,6)])
S3451_gene_norm <- as.data.frame(gene_list_norm[,c(1,7)])
S3453_gene_norm <- as.data.frame(gene_list_norm[,c(1,8)])
S3454_gene_norm <- as.data.frame(gene_list_norm[,c(1,9)])

S3458_gene_norm <- as.data.frame(gene_list_norm[,c(1,10)])
S3460_gene_norm <- as.data.frame(gene_list_norm[,c(1,11)])
S3461_gene_norm <- as.data.frame(gene_list_norm[,c(1,12)])
S3476_gene_norm <- as.data.frame(gene_list_norm[,c(1,13)])

colnames(S1825_gene_norm) <- c("Gene", "Count")
colnames(S1829_gene_norm) <- c("Gene", "Count")
colnames(S3444_gene_norm) <- c("Gene", "Count")
colnames(S3446_gene_norm) <- c("Gene", "Count")

colnames(S3450_gene_norm) <- c("Gene", "Count")
colnames(S3451_gene_norm) <- c("Gene", "Count")
colnames(S3453_gene_norm) <- c("Gene", "Count")
colnames(S3454_gene_norm) <- c("Gene", "Count")

colnames(S3458_gene_norm) <- c("Gene", "Count")
colnames(S3460_gene_norm) <- c("Gene", "Count")
colnames(S3461_gene_norm) <- c("Gene", "Count")
colnames(S3476_gene_norm) <- c("Gene", "Count")

S1825_gene_norm$Gene <- factor(S1825_gene_norm$Gene, 
                               levels = S1825_gene_norm$Gene)
S1829_gene_norm$Gene <- factor(S1829_gene_norm$Gene, 
                               levels = S1829_gene_norm$Gene)
S3444_gene_norm$Gene <- factor(S3444_gene_norm$Gene, 
                               levels = S3444_gene_norm$Gene)
S3446_gene_norm$Gene <- factor(S3446_gene_norm$Gene, 
                               levels = S3446_gene_norm$Gene)

S3450_gene_norm$Gene <- factor(S3450_gene_norm$Gene, 
                               levels = S3450_gene_norm$Gene)
S3451_gene_norm$Gene <- factor(S3451_gene_norm$Gene, 
                               levels = S3451_gene_norm$Gene)
S3453_gene_norm$Gene <- factor(S3453_gene_norm$Gene, 
                               levels = S3453_gene_norm$Gene)
S3454_gene_norm$Gene <- factor(S3454_gene_norm$Gene, 
                               levels = S3454_gene_norm$Gene)

S3458_gene_norm$Gene <- factor(S3458_gene_norm$Gene, 
                               levels = S3458_gene_norm$Gene)
S3460_gene_norm$Gene <- factor(S3460_gene_norm$Gene, 
                               levels = S3460_gene_norm$Gene)
S3461_gene_norm$Gene <- factor(S3461_gene_norm$Gene, 
                               levels = S3461_gene_norm$Gene)
S3476_gene_norm$Gene <- factor(S3476_gene_norm$Gene, 
                               levels = S3476_gene_norm$Gene)

#Plot individual normalised table
ggplot(S3476_gene_norm, aes(x=Gene, y=Count, fill=Gene)) + 
  geom_col(show.legend=FALSE) +
  scale_fill_manual(values = wes_palette(46, name="Darjeeling1", "continuous"))+
  xlab("Gene") +
  ylab("Number of SNPs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#.....

#Create combined normalised table
gene_list_norm_combined <- as.data.frame(rowSums(gene_list_norm))
gene_list_norm_combined$Gene <- rownames(gene_list_norm_combined)
gene_list_norm_combined <- gene_list_norm_combined[,c(2,1)]
colnames(gene_list_norm_combined) <- c("Gene", "Count")

gene_list_norm_combined$Gene <- factor(gene_list_norm_combined$Gene, 
                                       levels = gene_list_norm_combined$Gene)

#Plot combined normalised table
ggplot(gene_list_norm_combined, aes(x=Gene, y=Count, fill=Gene)) + 
  geom_col(show.legend=FALSE) +
  scale_fill_manual(values = wes_palette(46, name="Darjeeling1", "continuous"))+
  xlab("Gene") +
  ylab("Number of SNPs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

################################################################################
################################################################################ 
#Data output for further analysis
write.xlsx(sample_count, file="Mouse brain mtDNA mRNA-seq.xlsx", sheetName = "SNP count", append=T, row.names=F)
write.xlsx(library_count, file="Mouse brain mtDNA mRNA-seq.xlsx", sheetName = "Library size", append=T, row.names=F)
write.xlsx(norm_count, file="Mouse brain mtDNA mRNA-seq.xlsx", sheetName = "SNP count (normalised)", append=T, row.names=F)

write.xlsx(SNP_list, file="Mouse brain mtDNA mRNA-seq.xlsx", sheetName = "SNP alleles", append=T, row.names=F)

write.xlsx(features, file="Mouse brain mtDNA mRNA-seq.xlsx", sheetName = "mtDNA features", append=T, row.names=F)
write.xlsx(gene_list, file="Mouse brain mtDNA mRNA-seq.xlsx", sheetName = "Genes hits", append=T, row.names=F)
write.xlsx(gene_list_norm, file="Mouse brain mtDNA mRNA-seq.xlsx", sheetName = "Gene hits (normalised)", append=T, row.names=F)


