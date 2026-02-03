#### Load packages ####
library(readxl)
library(compositions)
library(tidyverse)
library(ccmm)
library(regmed)
library(gtools)

#### Process immune markers ####

immune <- read_xlsx("~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Data/Merged/Merged Immune Markers.xlsx", sheet = "metadata_immune markers")
immune <- immune[, names(immune) %in% c("Subject_ID", "Group_ID", "Serum_IL-2", "Serum_Fractalkine")]

immune <- immune[-which(immune$Subject_ID%in%c(43,51)),]

set.seed(823543)
immune_working <- immune[-c(sample(which(immune$Group_ID==1),1),
    sample(which(immune$Group_ID==7),1)),]


immune_1to2 <- immune_working[immune_working$Group_ID%in%c(1,2),]
immune_1to3 <- immune_working[immune_working$Group_ID%in%c(1,3),]
immune_1to4 <- immune_working[immune_working$Group_ID%in%c(1,4),]
immune_1to5 <- immune_working[immune_working$Group_ID%in%c(1,5),]
immune_1to10 <- immune_working[immune_working$Group_ID%in%c(1,10),]

immune_5to6 <- immune_working[immune_working$Group_ID%in%c(5,6),]
immune_5to7 <- immune_working[immune_working$Group_ID%in%c(5,7),]
immune_5to8 <- immune_working[immune_working$Group_ID%in%c(5,8),]

#### Handle Misc stuffs ####
### Set search grid for penalty parameters
lambda.grid <- seq(from = 0.4, to = 0.01, by = -0.05)

### Source
source("~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Code files/mvregmed.dat.check.R")

### Save times

times_sl <- times_gl <- matrix(ncol=3, nrow=16)
colnames(times_sl)<- colnames(times_gl) <- c("user", "system", "elapsed")
rownames(times_gl) <- c("gl_1to2_IL2", "gl_1to2_frac", "gl_1to3_IL2", "gl_1to3_frac", 
"gl_1to4_IL2", "gl_1to4_frac", "gl_1to5_IL2", "gl_1to5_frac", 
"gl_1to10_IL2", "gl_1to10_frac", "gl_5to6_IL2", "gl_5to6_frac", 
"gl_5to7_IL2", "gl_5to7_frac", "gl_5to8_IL2", "gl_5to8_frac")

rownames(times_sl) <- c("sl_1to2_IL2", "sl_1to2_frac", "sl_1to3_IL2", "sl_1to3_frac", 
"sl_1to4_IL2", "sl_1to4_frac", "sl_1to5_IL2", "sl_1to5_frac", 
"sl_1to10_IL2", "sl_1to10_frac", "sl_5to6_IL2", "sl_5to6_frac", 
"sl_5to7_IL2", "sl_5to7_frac", "sl_5to8_IL2", "sl_5to8_frac")


#### Process Genus-level microbiomes ####
micro_path <- "~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Data/Merged/"
microbiome <- read_xlsx(paste0(micro_path,"Merged Genus OTUs.xlsx", sheet="Oral vs Gut (No new 8)")

only_OTUs <- microbiome[,-c(1:5)]

# microbiome filtering
## transpose data for convenience
otu_table <- only_OTUs %>% t()

## remove samples with fewer than 100 reads
otu_table <- otu_table[,which(colSums(otu_table)>=100)]

## remove OTUs with fewer than 10 reads
otu_table <- otu_table[which(rowSums(otu_table)>=10),]

## remove OTUs present in less than 1% of samples
otu_table <- otu_table[which(rowSums(otu_table>0)>=ncol(otu_table)*0.01),]

## transpose back to original dimensions
OTUs_filtered <- t(otu_table)

## 5 taxa were filtered out; those included and dropped are stored below
included_taxa <- names(OTUs_filtered)
removed_taxa <- names(only_OTUs)[which(!(names(only_OTUs) %in% rownames(otu_table)))]

# add back ID columns
microbiome <- cbind(microbiome[,1:5], OTUs_filtered)

microbiome_oral <- microbiome[microbiome$Oral==1,]
microbiome_gut <- microbiome[microbiome$Oral==0,]

# Remove observations corresponding to the outlier immune markers (IDs 43 and 51)
microbiome_oral <- microbiome_oral[-which(microbiome_oral$Subject_ID%in%c(43,51)),]
microbiome_gut <- microbiome_gut[-which(microbiome_gut$Subject_ID%in%c(43,51)),]


# randomly drop 1 obs each from groups 1 and 7 in the oral microbiome to match dimensions of gut microbiome
set.seed(823543)
microbiome_oral = microbiome_oral[-c(sample(which(microbiome_oral$Group==1), 1), 
             sample(which(microbiome_oral$Group==7), 1)),]


# Convert rows to proportions
microbiome_oral[,-(1:5)] <- sweep(microbiome_oral[,-(1:5)], 1, rowSums(microbiome_oral[,-(1:5)]), FUN="/")
microbiome_gut[,-(1:5)] <- sweep(microbiome_gut[,-(1:5)], 1, rowSums(microbiome_gut[,-(1:5)]), FUN="/")

# Convert proportions to clr-scale
## These are not for use directly. Rather, we will want to apply clr to the closure of each
## subset. For example:
## medi1to2_oral <- microbiome_oral[]
microbiome_oral_clr <- cbind(microbiome_oral[,1:5], t(apply(microbiome_oral[,-(1:5)], 1, clr)))
microbiome_gut_clr <- cbind(microbiome_gut[,1:5], t(apply(microbiome_gut[,-(1:5)], 1, clr)))



#### Genus-level 1to2 ####
#### Oral
micro1to2_oral <- microbiome_oral[which(microbiome_oral$Group %in% c(1,2)),]
micro1to2_oral_clo <- apply(micro1to2_oral[-(1:5)], 1, clo)
micro1to2_oral_clr <- apply(micro1to2_oral_clo, 1, clr)
colnames(micro1to2_oral_clr) <- paste0(colnames(micro1to2_oral_clr), "Oral")

#### Gut
micro1to2_gut <- microbiome_oral[which(microbiome_gut$Group %in% c(1,2)),]
micro1to2_gut_clo <- apply(micro1to2_gut[-(1:5)], 1, clo)
micro1to2_gut_clr <- apply(micro1to2_gut_clo, 1, clr)
colnames(micro1to2_gut_clr) <- paste0(colnames(micro1to2_gut_clr), "Gut")

#### Combined
medi1to2_clr <- cbind(micro1to2_oral_clr[,-which(colSums(micro1to2_oral_clr)==0)], 
                      micro1to2_gut_clr[,-which(colSums(micro1to2_gut_clr)==0)])

#### Treatment and Outcome
x_1to2 <- c(rep(1,6), rep(0,8)) # code grp 1 = 1, grp 2 = 0
y1to2_IL2 <- immune_1to2$`Serum_IL-2` ; y1to2_frac <- immune_1to2$Serum_Fractalkine


#### IL-2
dat1to2_IL2_clr <- mvregmed.dat.check(x = x_1to2, y = y1to2_IL2, mediator = medi1to2_clr)
x1to2_IL2_clean <- dat1to2_IL2_clr$x; medi1to2_IL2_clean <- dat1to2_IL2_clr$mediator; y1to2_IL2_clean <- dat1to2_IL2_clr$y

s2 <- proc.time() 
gl_grid_sem1to2_IL2 <- regmed.grid(x = x1to2_IL2_clean, mediator = medi1to2_IL2_clean, y = y1to2_IL2_clean, 
           lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(gl_grid_sem1to2_IL2, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Genus level/gl_grid_sem1to2_IL2.Rdata")
times_gl[1,] <- e2[1:3]

#### Fractalkine
dat1to2_frac_clr <- mvregmed.dat.check(x = x_1to2, y = y1to2_frac, mediator = medi1to2_clr)
x1to2_frac_clean <- dat1to2_frac_clr$x; medi1to2_frac_clean <- dat1to2_frac_clr$mediator; y1to2_frac_clean <- dat1to2_frac_clr$y

s2 <- proc.time() 
gl_grid_sem1to2_frac <- regmed.grid(x = x1to2_frac_clean, mediator = medi1to2_frac_clean, y = y1to2_frac_clean, 
            lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(gl_grid_sem1to2_frac, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Genus level/gl_grid_sem1to2_frac.Rdata")
times_gl[2,] <- e2[1:3]


#### Genus-level 1to3 ####
#### Oral
micro1to3_oral <- microbiome_oral[which(microbiome_oral$Group %in% c(1,3)),]
micro1to3_oral_clo <- apply(micro1to3_oral[-(1:5)], 1, clo)
micro1to3_oral_clr <- apply(micro1to3_oral_clo, 1, clr)
colnames(micro1to3_oral_clr) <- paste0(colnames(micro1to3_oral_clr), "Oral")

#### Gut
micro1to3_gut <- microbiome_oral[which(microbiome_gut$Group %in% c(1,3)),]
micro1to3_gut_clo <- apply(micro1to3_gut[-(1:5)], 1, clo)
micro1to3_gut_clr <- apply(micro1to3_gut_clo, 1, clr)
colnames(micro1to3_gut_clr) <- paste0(colnames(micro1to3_gut_clr), "Gut")

#### Combined
medi1to3_clr <- cbind(micro1to3_oral_clr[,-which(colSums(micro1to3_oral_clr)==0)], 
                      micro1to3_gut_clr[,-which(colSums(micro1to3_gut_clr)==0)])

#### Treatment and Outcome
x_1to3 <- c(rep(1,6), rep(0,8)) # code grp 1 = 1, grp 2 = 0
y1to3_IL2 <- immune_1to3$`Serum_IL-2` ; y1to3_frac <- immune_1to3$Serum_Fractalkine


#### IL-2
dat1to3_IL2_clr <- mvregmed.dat.check(x = x_1to3, y = y1to3_IL2, mediator = medi1to3_clr)
x1to3_IL2_clean <- dat1to3_IL2_clr$x; medi1to3_IL2_clean <- dat1to3_IL2_clr$mediator; y1to3_IL2_clean <- dat1to3_IL2_clr$y

s2 <- proc.time() 
gl_grid_sem1to3_IL2 <- regmed.grid(x = x1to3_IL2_clean, mediator = medi1to3_IL2_clean, y = y1to3_IL2_clean, 
           lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(gl_grid_sem1to3_IL2, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Genus level/gl_grid_sem1to3_IL2.Rdata")
times_gl[3,] <- e2[1:3]

#### Fractalkine
dat1to3_frac_clr <- mvregmed.dat.check(x = x_1to3, y = y1to3_frac, mediator = medi1to3_clr)
x1to3_frac_clean <- dat1to3_frac_clr$x; medi1to3_frac_clean <- dat1to3_frac_clr$mediator; y1to3_frac_clean <- dat1to3_frac_clr$y

s2 <- proc.time() 
gl_grid_sem1to3_frac <- regmed.grid(x = x1to3_frac_clean, mediator = medi1to3_frac_clean, y = y1to3_frac_clean, 
            lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(gl_grid_sem1to3_frac, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Genus level/gl_grid_sem1to3_frac.Rdata")
times_gl[4,] <- e2[1:3]

#### Genus-level 1to4 ####
#### Oral
micro1to4_oral <- microbiome_oral[which(microbiome_oral$Group %in% c(1,4)),]
micro1to4_oral_clo <- apply(micro1to4_oral[-(1:5)], 1, clo)
micro1to4_oral_clr <- apply(micro1to4_oral_clo, 1, clr)
colnames(micro1to4_oral_clr) <- paste0(colnames(micro1to4_oral_clr), "Oral")

#### Gut
micro1to4_gut <- microbiome_oral[which(microbiome_gut$Group %in% c(1,4)),]
micro1to4_gut_clo <- apply(micro1to4_gut[-(1:5)], 1, clo)
micro1to4_gut_clr <- apply(micro1to4_gut_clo, 1, clr)
colnames(micro1to4_gut_clr) <- paste0(colnames(micro1to4_gut_clr), "Gut")

#### Combined
medi1to4_clr <- cbind(micro1to4_oral_clr[,-which(colSums(micro1to4_oral_clr)==0)], 
                      micro1to4_gut_clr[,-which(colSums(micro1to4_gut_clr)==0)])

#### Treatment and Outcome
x_1to4 <- c(rep(1,6), rep(0,8)) # code grp 1 = 1, grp 2 = 0
y1to4_IL2 <- immune_1to4$`Serum_IL-2` ; y1to4_frac <- immune_1to4$Serum_Fractalkine


#### IL-2
dat1to4_IL2_clr <- mvregmed.dat.check(x = x_1to4, y = y1to4_IL2, mediator = medi1to4_clr)
x1to4_IL2_clean <- dat1to4_IL2_clr$x; medi1to4_IL2_clean <- dat1to4_IL2_clr$mediator; y1to4_IL2_clean <- dat1to4_IL2_clr$y

s2 <- proc.time() 
gl_grid_sem1to4_IL2 <- regmed.grid(x = x1to4_IL2_clean, mediator = medi1to4_IL2_clean, y = y1to4_IL2_clean, 
           lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(gl_grid_sem1to4_IL2, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Genus level/gl_grid_sem1to4_IL2.Rdata")
times_gl[5,] <- e2[1:3]

#### Fractalkine
dat1to4_frac_clr <- mvregmed.dat.check(x = x_1to4, y = y1to4_frac, mediator = medi1to4_clr)
x1to4_frac_clean <- dat1to4_frac_clr$x; medi1to4_frac_clean <- dat1to4_frac_clr$mediator; y1to4_frac_clean <- dat1to4_frac_clr$y

s2 <- proc.time() 
gl_grid_sem1to4_frac <- regmed.grid(x = x1to4_frac_clean, mediator = medi1to4_frac_clean, y = y1to4_frac_clean, 
            lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(gl_grid_sem1to4_frac, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Genus level/gl_grid_sem1to4_frac.Rdata")
times_gl[6,] <- e2[1:3]

#### Genus-level 1to5 ####
#### Oral
micro1to5_oral <- microbiome_oral[which(microbiome_oral$Group %in% c(1,5)),]
micro1to5_oral_clo <- apply(micro1to5_oral[-(1:5)], 1, clo)
micro1to5_oral_clr <- apply(micro1to5_oral_clo, 1, clr)
colnames(micro1to5_oral_clr) <- paste0(colnames(micro1to5_oral_clr), "Oral")

#### Gut
micro1to5_gut <- microbiome_oral[which(microbiome_gut$Group %in% c(1,5)),]
micro1to5_gut_clo <- apply(micro1to5_gut[-(1:5)], 1, clo)
micro1to5_gut_clr <- apply(micro1to5_gut_clo, 1, clr)
colnames(micro1to5_gut_clr) <- paste0(colnames(micro1to5_gut_clr), "Gut")

#### Combined
medi1to5_clr <- cbind(micro1to5_oral_clr[,-which(colSums(micro1to5_oral_clr)==0)], 
                      micro1to5_gut_clr[,-which(colSums(micro1to5_gut_clr)==0)])

#### Treatment and Outcome
x_1to5 <- c(rep(1,6), rep(0,8)) # code grp 1 = 1, grp 2 = 0
y1to5_IL2 <- immune_1to5$`Serum_IL-2` ; y1to5_frac <- immune_1to5$Serum_Fractalkine


#### IL-2
dat1to5_IL2_clr <- mvregmed.dat.check(x = x_1to5, y = y1to5_IL2, mediator = medi1to5_clr)
x1to5_IL2_clean <- dat1to5_IL2_clr$x; medi1to5_IL2_clean <- dat1to5_IL2_clr$mediator; y1to5_IL2_clean <- dat1to5_IL2_clr$y

s2 <- proc.time() 
gl_grid_sem1to5_IL2 <- regmed.grid(x = x1to5_IL2_clean, mediator = medi1to5_IL2_clean, y = y1to5_IL2_clean, 
           lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(gl_grid_sem1to5_IL2, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Genus level/gl_grid_sem1to5_IL2.Rdata")
times_gl[7,] <- e2[1:3]

#### Fractalkine
dat1to5_frac_clr <- mvregmed.dat.check(x = x_1to5, y = y1to5_frac, mediator = medi1to5_clr)
x1to5_frac_clean <- dat1to5_frac_clr$x; medi1to5_frac_clean <- dat1to5_frac_clr$mediator; y1to5_frac_clean <- dat1to5_frac_clr$y

s2 <- proc.time() 
gl_grid_sem1to5_frac <- regmed.grid(x = x1to5_frac_clean, mediator = medi1to5_frac_clean, y = y1to5_frac_clean, 
            lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(gl_grid_sem1to5_frac, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Genus level/gl_grid_sem1to5_frac.Rdata")
times_gl[8,] <- e2[1:3]

#### Genus-level 1to10 ####
#### Oral
micro1to10_oral <- microbiome_oral[which(microbiome_oral$Group %in% c(1,2)),]
micro1to10_oral_clo <- apply(micro1to10_oral[-(1:5)], 1, clo)
micro1to10_oral_clr <- apply(micro1to10_oral_clo, 1, clr)
colnames(micro1to10_oral_clr) <- paste0(colnames(micro1to10_oral_clr), "Oral")

#### Gut
micro1to10_gut <- microbiome_oral[which(microbiome_gut$Group %in% c(1,2)),]
micro1to10_gut_clo <- apply(micro1to10_gut[-(1:5)], 1, clo)
micro1to10_gut_clr <- apply(micro1to10_gut_clo, 1, clr)
colnames(micro1to10_gut_clr) <- paste0(colnames(micro1to10_gut_clr), "Gut")

#### Combined
medi1to10_clr <- cbind(micro1to10_oral_clr[,-which(colSums(micro1to10_oral_clr)==0)], 
                       micro1to10_gut_clr[,-which(colSums(micro1to10_gut_clr)==0)])

#### Treatment and Outcome
x_1to10 <- c(rep(1,6), rep(0,8)) # code grp 1 = 1, grp 2 = 0
y1to10_IL2 <- immune_1to10$`Serum_IL-2` ; y1to10_frac <- immune_1to10$Serum_Fractalkine


#### IL-2
dat1to10_IL2_clr <- mvregmed.dat.check(x = x_1to10, y = y1to10_IL2, mediator = medi1to10_clr)
x1to10_IL2_clean <- dat1to10_IL2_clr$x; medi1to10_IL2_clean <- dat1to10_IL2_clr$mediator; y1to10_IL2_clean <- dat1to10_IL2_clr$y

s2 <- proc.time() 
gl_grid_sem1to10_IL2 <- regmed.grid(x = x1to10_IL2_clean, mediator = medi1to10_IL2_clean, y = y1to10_IL2_clean, 
            lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(gl_grid_sem1to10_IL2, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Genus level/gl_grid_sem1to10_IL2.Rdata")
times_gl[9,] <- e2[1:3]

#### Fractalkine
dat1to10_frac_clr <- mvregmed.dat.check(x = x_1to10, y = y1to10_frac, mediator = medi1to10_clr)
x1to10_frac_clean <- dat1to10_frac_clr$x; medi1to10_frac_clean <- dat1to10_frac_clr$mediator; y1to10_frac_clean <- dat1to10_frac_clr$y

s2 <- proc.time() 
gl_grid_sem1to10_frac <- regmed.grid(x = x1to10_frac_clean, mediator = medi1to10_frac_clean, y = y1to10_frac_clean, 
             lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(gl_grid_sem1to10_frac, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Genus level/gl_grid_sem1to10_frac.Rdata")
times_gl[10,] <- e2[1:3]

#### Genus-level 5to6 ####
#### Oral
micro5to6_oral <- microbiome_oral[which(microbiome_oral$Group %in% c(5,6)),]
micro5to6_oral_clo <- apply(micro5to6_oral[-(1:5)], 1, clo)
micro5to6_oral_clr <- apply(micro5to6_oral_clo, 1, clr)
colnames(micro5to6_oral_clr) <- paste0(colnames(micro5to6_oral_clr), "Oral")

#### Gut
micro5to6_gut <- microbiome_oral[which(microbiome_gut$Group %in% c(5,6)),]
micro5to6_gut_clo <- apply(micro5to6_gut[-(1:5)], 1, clo)
micro5to6_gut_clr <- apply(micro5to6_gut_clo, 1, clr)
colnames(micro5to6_gut_clr) <- paste0(colnames(micro5to6_gut_clr), "Gut")

#### Combined
medi5to6_clr <- cbind(micro5to6_oral_clr[,-which(colSums(micro5to6_oral_clr)==0)], 
                      micro5to6_gut_clr[,-which(colSums(micro5to6_gut_clr)==0)])

#### Treatment and Outcome
x_5to6 <- c(rep(1,8), rep(0,7)) # code grp 1 = 1, grp 2 = 0
y5to6_IL2 <- immune_5to6$`Serum_IL-2` ; y5to6_frac <- immune_5to6$Serum_Fractalkine


#### IL-2
dat5to6_IL2_clr <- mvregmed.dat.check(x = x_5to6, y = y5to6_IL2, mediator = medi5to6_clr)
x5to6_IL2_clean <- dat5to6_IL2_clr$x; medi5to6_IL2_clean <- dat5to6_IL2_clr$mediator; y5to6_IL2_clean <- dat5to6_IL2_clr$y

s2 <- proc.time() 
gl_grid_sem5to6_IL2 <- regmed.grid(x = x5to6_IL2_clean, mediator = medi5to6_IL2_clean, y = y5to6_IL2_clean, 
           lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(gl_grid_sem5to6_IL2, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Genus level/gl_grid_sem5to6_IL2.Rdata")
times_gl[11,] <- e2[1:3]

#### Fractalkine
dat5to6_frac_clr <- mvregmed.dat.check(x = x_5to6, y = y5to6_frac, mediator = medi5to6_clr)
x5to6_frac_clean <- dat5to6_frac_clr$x; medi5to6_frac_clean <- dat5to6_frac_clr$mediator; y5to6_frac_clean <- dat5to6_frac_clr$y

s2 <- proc.time() 
gl_grid_sem5to6_frac <- regmed.grid(x = x5to6_frac_clean, mediator = medi5to6_frac_clean, y = y5to6_frac_clean, 
            lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(gl_grid_sem5to6_frac, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Genus level/gl_grid_sem5to6_frac.Rdata")
times_gl[12,] <- e2[1:3]

#### Genus-level 5to7 ####
#### Oral
micro5to7_oral <- microbiome_oral[which(microbiome_oral$Group %in% c(5,7)),]
micro5to7_oral_clo <- apply(micro5to7_oral[-(1:5)], 1, clo)
micro5to7_oral_clr <- apply(micro5to7_oral_clo, 1, clr)
colnames(micro5to7_oral_clr) <- paste0(colnames(micro5to7_oral_clr), "Oral")

#### Gut
micro5to7_gut <- microbiome_oral[which(microbiome_gut$Group %in% c(5,7)),]
micro5to7_gut_clo <- apply(micro5to7_gut[-(1:5)], 1, clo)
micro5to7_gut_clr <- apply(micro5to7_gut_clo, 1, clr)
colnames(micro5to7_gut_clr) <- paste0(colnames(micro5to7_gut_clr), "Gut")

#### Combined
medi5to7_clr <- cbind(micro5to7_oral_clr[,-which(colSums(micro5to7_oral_clr)==0)], 
                      micro5to7_gut_clr[,-which(colSums(micro5to7_gut_clr)==0)])

#### Treatment and Outcome
x_5to7 <- c(rep(1,8), rep(0,6)) # code grp 1 = 1, grp 2 = 0
y5to7_IL2 <- immune_5to7$`Serum_IL-2` ; y5to7_frac <- immune_5to7$Serum_Fractalkine


#### IL-2
dat5to7_IL2_clr <- mvregmed.dat.check(x = x_5to7, y = y5to7_IL2, mediator = medi5to7_clr)
x5to7_IL2_clean <- dat5to7_IL2_clr$x; medi5to7_IL2_clean <- dat5to7_IL2_clr$mediator; y5to7_IL2_clean <- dat5to7_IL2_clr$y

s2 <- proc.time() 
gl_grid_sem5to7_IL2 <- regmed.grid(x = x5to7_IL2_clean, mediator = medi5to7_IL2_clean, y = y5to7_IL2_clean, 
           lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(gl_grid_sem5to7_IL2, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Genus level/gl_grid_sem5to7_IL2.Rdata")
times_gl[13,] <- e2[1:3]

#### Fractalkine
dat5to7_frac_clr <- mvregmed.dat.check(x = x_5to7, y = y5to7_frac, mediator = medi5to7_clr)
x5to7_frac_clean <- dat5to7_frac_clr$x; medi5to7_frac_clean <- dat5to7_frac_clr$mediator; y5to7_frac_clean <- dat5to7_frac_clr$y

s2 <- proc.time() 
gl_grid_sem5to7_frac <- regmed.grid(x = x5to7_frac_clean, mediator = medi5to7_frac_clean, y = y5to7_frac_clean, 
            lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(gl_grid_sem5to7_frac, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Genus level/gl_grid_sem5to7_frac.Rdata")
times_gl[14,] <- e2[1:3]

#### Genus-level 5to8 ####
#### Oral
micro5to8_oral <- microbiome_oral[which(microbiome_oral$Group %in% c(5,8)),]
micro5to8_oral_clo <- apply(micro5to8_oral[-(1:5)], 1, clo)
micro5to8_oral_clr <- apply(micro5to8_oral_clo, 1, clr)
colnames(micro5to8_oral_clr) <- paste0(colnames(micro5to8_oral_clr), "Oral")

#### Gut
micro5to8_gut <- microbiome_oral[which(microbiome_gut$Group %in% c(5,8)),]
micro5to8_gut_clo <- apply(micro5to8_gut[-(1:5)], 1, clo)
micro5to8_gut_clr <- apply(micro5to8_gut_clo, 1, clr)
colnames(micro5to8_gut_clr) <- paste0(colnames(micro5to8_gut_clr), "Gut")

#### Combined
medi5to8_clr <- cbind(micro5to8_oral_clr[,-which(colSums(micro5to8_oral_clr)==0)], 
                      micro5to8_gut_clr[,-which(colSums(micro5to8_gut_clr)==0)])

#### Treatment and Outcome
x_5to8 <- c(rep(1,8), rep(0,8)) # code grp 1 = 1, grp 2 = 0
y5to8_IL2 <- immune_5to8$`Serum_IL-2` ; y5to8_frac <- immune_5to8$Serum_Fractalkine


#### IL-2
dat5to8_IL2_clr <- mvregmed.dat.check(x = x_5to8, y = y5to8_IL2, mediator = medi5to8_clr)
x5to8_IL2_clean <- dat5to8_IL2_clr$x; medi5to8_IL2_clean <- dat5to8_IL2_clr$mediator; y5to8_IL2_clean <- dat5to8_IL2_clr$y

s2 <- proc.time() 
gl_grid_sem5to8_IL2 <- regmed.grid(x = x5to8_IL2_clean, mediator = medi5to8_IL2_clean, y = y5to8_IL2_clean, 
           lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(gl_grid_sem5to8_IL2, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Genus level/gl_grid_sem5to8_IL2.Rdata")
times_gl[15,] <- e2[1:3]

#### Fractalkine
dat5to8_frac_clr <- mvregmed.dat.check(x = x_5to8, y = y5to8_frac, mediator = medi5to8_clr)
x5to8_frac_clean <- dat5to8_frac_clr$x; medi5to8_frac_clean <- dat5to8_frac_clr$mediator; y5to8_frac_clean <- dat5to8_frac_clr$y

s2 <- proc.time() 
gl_grid_sem5to8_frac <- regmed.grid(x = x5to8_frac_clean, mediator = medi5to8_frac_clean, y = y5to8_frac_clean, 
            lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(gl_grid_sem5to8_frac, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Genus level/gl_grid_sem5to8_frac.Rdata")
times_gl[16,] <- e2[1:3]


#### Process Species-level microbiomes ####
microbiome <- read_xlsx("~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Data/Merged/Merged Species OTUs.xlsx", 
sheet = "Oral vs Gut (no new 8)")



only_OTUs <- microbiome[,-c(1:5)]

# microbiome filtering
## transpose data for convenience
otu_table <- only_OTUs %>% t()

## remove samples with fewer than 100 reads
otu_table <- otu_table[,which(colSums(otu_table)>=100)]

## remove OTUs with fewer than 10 reads
otu_table <- otu_table[which(rowSums(otu_table)>=10),]

## remove OTUs present in less than 1% of samples
otu_table <- otu_table[which(rowSums(otu_table>0)>=ncol(otu_table)*0.01),]

## transpose back to original dimensions
OTUs_filtered <- t(otu_table)

## 13 taxa were filtered out; those included and dropped are stored below
# included_taxa <- names(OTUs_filtered)
removed_taxa <- names(only_OTUs)[which(!(names(only_OTUs) %in% rownames(otu_table)))]

# add back ID columns
microbiome <- cbind(microbiome[,1:5], OTUs_filtered)

microbiome_oral <- microbiome[microbiome$Oral==1,]
microbiome_gut <- microbiome[microbiome$Oral==0,]

# Remove observations corresponding to the outlier immune markers (IDs 43 and 51)
microbiome_oral <- microbiome_oral[-which(microbiome_oral$Subject_ID%in%c(43,51)),]
microbiome_gut <- microbiome_gut[-which(microbiome_gut$Subject_ID%in%c(43,51)),]


# randomly drop 1 obs each from groups 1 and 7 in the oral microbiome to match dimensions of gut microbiome
set.seed(823543)
microbiome_oral = microbiome_oral[-c(sample(which(microbiome_oral$Group==1), 1), 
             sample(which(microbiome_oral$Group==7), 1)),]



# Convert rows to proportions
microbiome_oral[,-(1:5)] <- sweep(microbiome_oral[,-(1:5)], 1, rowSums(microbiome_oral[,-(1:5)]), FUN="/")
microbiome_gut[,-(1:5)] <- sweep(microbiome_gut[,-(1:5)], 1, rowSums(microbiome_gut[,-(1:5)]), FUN="/")

# Convert proportions to clr-scale
## These are not for use directly. Rather, we will want to apply clr to the closure of each
## subset. For example:
## medi1to2_oral <- microbiome_oral[]
microbiome_oral_clr <- cbind(microbiome_oral[,1:5], t(apply(microbiome_oral[,-(1:5)], 1, clr)))
microbiome_gut_clr <- cbind(microbiome_gut[,1:5], t(apply(microbiome_gut[,-(1:5)], 1, clr)))

#### Species-level 1to2 ####
#### Oral
micro1to2_oral <- microbiome_oral[which(microbiome_oral$Group %in% c(1,2)),]
micro1to2_oral_clo <- apply(micro1to2_oral[-(1:5)], 1, clo)
micro1to2_oral_clr <- apply(micro1to2_oral_clo, 1, clr)
colnames(micro1to2_oral_clr) <- paste0(colnames(micro1to2_oral_clr), "Oral")

#### Gut
micro1to2_gut <- microbiome_oral[which(microbiome_gut$Group %in% c(1,2)),]
micro1to2_gut_clo <- apply(micro1to2_gut[-(1:5)], 1, clo)
micro1to2_gut_clr <- apply(micro1to2_gut_clo, 1, clr)
colnames(micro1to2_gut_clr) <- paste0(colnames(micro1to2_gut_clr), "Gut")

#### Combined
medi1to2_clr <- cbind(micro1to2_oral_clr[,-which(colSums(micro1to2_oral_clr)==0)], 
                      micro1to2_gut_clr[,-which(colSums(micro1to2_gut_clr)==0)])

#### Treatment and Outcome
x_1to2 <- c(rep(1,6), rep(0,8)) # code grp 1 = 1, grp 2 = 0
y1to2_IL2 <- immune_1to2$`Serum_IL-2` ; y1to2_frac <- immune_1to2$Serum_Fractalkine


#### IL-2
dat1to2_IL2_clr <- mvregmed.dat.check(x = x_1to2, y = y1to2_IL2, mediator = medi1to2_clr)
x1to2_IL2_clean <- dat1to2_IL2_clr$x; medi1to2_IL2_clean <- dat1to2_IL2_clr$mediator; y1to2_IL2_clean <- dat1to2_IL2_clr$y

s2 <- proc.time() 
sl_grid_sem1to2_IL2 <- regmed.grid(x = x1to2_IL2_clean, mediator = medi1to2_IL2_clean, y = y1to2_IL2_clean, 
           lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(sl_grid_sem1to2_IL2, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Species level/sl_grid_sem1to2_IL2.Rdata")
times_sl[1,] <- e2[1:3]

#### Fractalkine
dat1to2_frac_clr <- mvregmed.dat.check(x = x_1to2, y = y1to2_frac, mediator = medi1to2_clr)
x1to2_frac_clean <- dat1to2_frac_clr$x; medi1to2_frac_clean <- dat1to2_frac_clr$mediator; y1to2_frac_clean <- dat1to2_frac_clr$y

s2 <- proc.time() 
sl_grid_sem1to2_frac <- regmed.grid(x = x1to2_frac_clean, mediator = medi1to2_frac_clean, y = y1to2_frac_clean, 
            lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(sl_grid_sem1to2_frac, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Species level/sl_grid_sem1to2_frac.Rdata")
times_sl[2,] <- e2[1:3]


#### Species-level 1to3 ####
#### Oral
micro1to3_oral <- microbiome_oral[which(microbiome_oral$Group %in% c(1,3)),]
micro1to3_oral_clo <- apply(micro1to3_oral[-(1:5)], 1, clo)
micro1to3_oral_clr <- apply(micro1to3_oral_clo, 1, clr)
colnames(micro1to3_oral_clr) <- paste0(colnames(micro1to3_oral_clr), "Oral")

#### Gut
micro1to3_gut <- microbiome_oral[which(microbiome_gut$Group %in% c(1,3)),]
micro1to3_gut_clo <- apply(micro1to3_gut[-(1:5)], 1, clo)
micro1to3_gut_clr <- apply(micro1to3_gut_clo, 1, clr)
colnames(micro1to3_gut_clr) <- paste0(colnames(micro1to3_gut_clr), "Gut")

#### Combined
medi1to3_clr <- cbind(micro1to3_oral_clr[,-which(colSums(micro1to3_oral_clr)==0)], 
                      micro1to3_gut_clr[,-which(colSums(micro1to3_gut_clr)==0)])

#### Treatment and Outcome
x_1to3 <- c(rep(1,6), rep(0,8)) # code grp 1 = 1, grp 2 = 0
y1to3_IL2 <- immune_1to3$`Serum_IL-2` ; y1to3_frac <- immune_1to3$Serum_Fractalkine


#### IL-2
dat1to3_IL2_clr <- mvregmed.dat.check(x = x_1to3, y = y1to3_IL2, mediator = medi1to3_clr)
x1to3_IL2_clean <- dat1to3_IL2_clr$x; medi1to3_IL2_clean <- dat1to3_IL2_clr$mediator; y1to3_IL2_clean <- dat1to3_IL2_clr$y

s2 <- proc.time() 
sl_grid_sem1to3_IL2 <- regmed.grid(x = x1to3_IL2_clean, mediator = medi1to3_IL2_clean, y = y1to3_IL2_clean, 
           lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(sl_grid_sem1to3_IL2, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Species level/sl_grid_sem1to3_IL2.Rdata")
times_sl[3,] <- e2[1:3]

#### Fractalkine
dat1to3_frac_clr <- mvregmed.dat.check(x = x_1to3, y = y1to3_frac, mediator = medi1to3_clr)
x1to3_frac_clean <- dat1to3_frac_clr$x; medi1to3_frac_clean <- dat1to3_frac_clr$mediator; y1to3_frac_clean <- dat1to3_frac_clr$y

s2 <- proc.time() 
sl_grid_sem1to3_frac <- regmed.grid(x = x1to3_frac_clean, mediator = medi1to3_frac_clean, y = y1to3_frac_clean, 
            lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(sl_grid_sem1to3_frac, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Species level/sl_grid_sem1to3_frac.Rdata")
times_sl[4,] <- e2[1:3]

#### Species-level 1to4 ####
#### Oral
micro1to4_oral <- microbiome_oral[which(microbiome_oral$Group %in% c(1,4)),]
micro1to4_oral_clo <- apply(micro1to4_oral[-(1:5)], 1, clo)
micro1to4_oral_clr <- apply(micro1to4_oral_clo, 1, clr)
colnames(micro1to4_oral_clr) <- paste0(colnames(micro1to4_oral_clr), "Oral")

#### Gut
micro1to4_gut <- microbiome_oral[which(microbiome_gut$Group %in% c(1,4)),]
micro1to4_gut_clo <- apply(micro1to4_gut[-(1:5)], 1, clo)
micro1to4_gut_clr <- apply(micro1to4_gut_clo, 1, clr)
colnames(micro1to4_gut_clr) <- paste0(colnames(micro1to4_gut_clr), "Gut")

#### Combined
medi1to4_clr <- cbind(micro1to4_oral_clr[,-which(colSums(micro1to4_oral_clr)==0)], 
                      micro1to4_gut_clr[,-which(colSums(micro1to4_gut_clr)==0)])

#### Treatment and Outcome
x_1to4 <- c(rep(1,6), rep(0,8)) # code grp 1 = 1, grp 2 = 0
y1to4_IL2 <- immune_1to4$`Serum_IL-2` ; y1to4_frac <- immune_1to4$Serum_Fractalkine


#### IL-2
dat1to4_IL2_clr <- mvregmed.dat.check(x = x_1to4, y = y1to4_IL2, mediator = medi1to4_clr)
x1to4_IL2_clean <- dat1to4_IL2_clr$x; medi1to4_IL2_clean <- dat1to4_IL2_clr$mediator; y1to4_IL2_clean <- dat1to4_IL2_clr$y

s2 <- proc.time() 
sl_grid_sem1to4_IL2 <- regmed.grid(x = x1to4_IL2_clean, mediator = medi1to4_IL2_clean, y = y1to4_IL2_clean, 
           lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(sl_grid_sem1to4_IL2, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Species level/sl_grid_sem1to4_IL2.Rdata")
times_sl[5,] <- e2[1:3]

#### Fractalkine
dat1to4_frac_clr <- mvregmed.dat.check(x = x_1to4, y = y1to4_frac, mediator = medi1to4_clr)
x1to4_frac_clean <- dat1to4_frac_clr$x; medi1to4_frac_clean <- dat1to4_frac_clr$mediator; y1to4_frac_clean <- dat1to4_frac_clr$y

s2 <- proc.time() 
sl_grid_sem1to4_frac <- regmed.grid(x = x1to4_frac_clean, mediator = medi1to4_frac_clean, y = y1to4_frac_clean, 
            lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(sl_grid_sem1to4_frac, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Species level/sl_grid_sem1to4_frac.Rdata")
times_sl[6,] <- e2[1:3]

#### Species-level 1to5 ####
#### Oral
micro1to5_oral <- microbiome_oral[which(microbiome_oral$Group %in% c(1,5)),]
micro1to5_oral_clo <- apply(micro1to5_oral[-(1:5)], 1, clo)
micro1to5_oral_clr <- apply(micro1to5_oral_clo, 1, clr)
colnames(micro1to5_oral_clr) <- paste0(colnames(micro1to5_oral_clr), "Oral")

#### Gut
micro1to5_gut <- microbiome_oral[which(microbiome_gut$Group %in% c(1,5)),]
micro1to5_gut_clo <- apply(micro1to5_gut[-(1:5)], 1, clo)
micro1to5_gut_clr <- apply(micro1to5_gut_clo, 1, clr)
colnames(micro1to5_gut_clr) <- paste0(colnames(micro1to5_gut_clr), "Gut")

#### Combined
medi1to5_clr <- cbind(micro1to5_oral_clr[,-which(colSums(micro1to5_oral_clr)==0)], 
                      micro1to5_gut_clr[,-which(colSums(micro1to5_gut_clr)==0)])

#### Treatment and Outcome
x_1to5 <- c(rep(1,6), rep(0,8)) # code grp 1 = 1, grp 2 = 0
y1to5_IL2 <- immune_1to5$`Serum_IL-2` ; y1to5_frac <- immune_1to5$Serum_Fractalkine


#### IL-2
dat1to5_IL2_clr <- mvregmed.dat.check(x = x_1to5, y = y1to5_IL2, mediator = medi1to5_clr)
x1to5_IL2_clean <- dat1to5_IL2_clr$x; medi1to5_IL2_clean <- dat1to5_IL2_clr$mediator; y1to5_IL2_clean <- dat1to5_IL2_clr$y

s2 <- proc.time() 
sl_grid_sem1to5_IL2 <- regmed.grid(x = x1to5_IL2_clean, mediator = medi1to5_IL2_clean, y = y1to5_IL2_clean, 
           lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(sl_grid_sem1to5_IL2, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Species level/sl_grid_sem1to5_IL2.Rdata")
times_sl[7,] <- e2[1:3]

#### Fractalkine
dat1to5_frac_clr <- mvregmed.dat.check(x = x_1to5, y = y1to5_frac, mediator = medi1to5_clr)
x1to5_frac_clean <- dat1to5_frac_clr$x; medi1to5_frac_clean <- dat1to5_frac_clr$mediator; y1to5_frac_clean <- dat1to5_frac_clr$y

s2 <- proc.time() 
sl_grid_sem1to5_frac <- regmed.grid(x = x1to5_frac_clean, mediator = medi1to5_frac_clean, y = y1to5_frac_clean, 
            lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(sl_grid_sem1to5_frac, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Species level/sl_grid_sem1to5_frac.Rdata")
times_sl[8,] <- e2[1:3]

#### Species-level 1to10 ####
#### Oral
micro1to10_oral <- microbiome_oral[which(microbiome_oral$Group %in% c(1,2)),]
micro1to10_oral_clo <- apply(micro1to10_oral[-(1:5)], 1, clo)
micro1to10_oral_clr <- apply(micro1to10_oral_clo, 1, clr)
colnames(micro1to10_oral_clr) <- paste0(colnames(micro1to10_oral_clr), "Oral")

#### Gut
micro1to10_gut <- microbiome_oral[which(microbiome_gut$Group %in% c(1,2)),]
micro1to10_gut_clo <- apply(micro1to10_gut[-(1:5)], 1, clo)
micro1to10_gut_clr <- apply(micro1to10_gut_clo, 1, clr)
colnames(micro1to10_gut_clr) <- paste0(colnames(micro1to10_gut_clr), "Gut")

#### Combined
medi1to10_clr <- cbind(micro1to10_oral_clr[,-which(colSums(micro1to10_oral_clr)==0)], 
                       micro1to10_gut_clr[,-which(colSums(micro1to10_gut_clr)==0)])

#### Treatment and Outcome
x_1to10 <- c(rep(1,6), rep(0,8)) # code grp 1 = 1, grp 2 = 0
y1to10_IL2 <- immune_1to10$`Serum_IL-2` ; y1to10_frac <- immune_1to10$Serum_Fractalkine


#### IL-2
dat1to10_IL2_clr <- mvregmed.dat.check(x = x_1to10, y = y1to10_IL2, mediator = medi1to10_clr)
x1to10_IL2_clean <- dat1to10_IL2_clr$x; medi1to10_IL2_clean <- dat1to10_IL2_clr$mediator; y1to10_IL2_clean <- dat1to10_IL2_clr$y

s2 <- proc.time() 
sl_grid_sem1to10_IL2 <- regmed.grid(x = x1to10_IL2_clean, mediator = medi1to10_IL2_clean, y = y1to10_IL2_clean, 
            lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(sl_grid_sem1to10_IL2, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Species level/sl_grid_sem1to10_IL2.Rdata")
times_sl[9,] <- e2[1:3]

#### Fractalkine
dat1to10_frac_clr <- mvregmed.dat.check(x = x_1to10, y = y1to10_frac, mediator = medi1to10_clr)
x1to10_frac_clean <- dat1to10_frac_clr$x; medi1to10_frac_clean <- dat1to10_frac_clr$mediator; y1to10_frac_clean <- dat1to10_frac_clr$y

s2 <- proc.time() 
sl_grid_sem1to10_frac <- regmed.grid(x = x1to10_frac_clean, mediator = medi1to10_frac_clean, y = y1to10_frac_clean, 
             lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(sl_grid_sem1to10_frac, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Species level/sl_grid_sem1to10_frac.Rdata")
times_sl[10,] <- e2[1:3]

#### Species-level 5to6 ####
#### Oral
micro5to6_oral <- microbiome_oral[which(microbiome_oral$Group %in% c(5,6)),]
micro5to6_oral_clo <- apply(micro5to6_oral[-(1:5)], 1, clo)
micro5to6_oral_clr <- apply(micro5to6_oral_clo, 1, clr)
colnames(micro5to6_oral_clr) <- paste0(colnames(micro5to6_oral_clr), "Oral")

#### Gut
micro5to6_gut <- microbiome_oral[which(microbiome_gut$Group %in% c(5,6)),]
micro5to6_gut_clo <- apply(micro5to6_gut[-(1:5)], 1, clo)
micro5to6_gut_clr <- apply(micro5to6_gut_clo, 1, clr)
colnames(micro5to6_gut_clr) <- paste0(colnames(micro5to6_gut_clr), "Gut")

#### Combined
medi5to6_clr <- cbind(micro5to6_oral_clr[,-which(colSums(micro5to6_oral_clr)==0)], 
                      micro5to6_gut_clr[,-which(colSums(micro5to6_gut_clr)==0)])

#### Treatment and Outcome
x_5to6 <- c(rep(1,8), rep(0,7)) # code grp 1 = 1, grp 2 = 0
y5to6_IL2 <- immune_5to6$`Serum_IL-2` ; y5to6_frac <- immune_5to6$Serum_Fractalkine


#### IL-2
dat5to6_IL2_clr <- mvregmed.dat.check(x = x_5to6, y = y5to6_IL2, mediator = medi5to6_clr)
x5to6_IL2_clean <- dat5to6_IL2_clr$x; medi5to6_IL2_clean <- dat5to6_IL2_clr$mediator; y5to6_IL2_clean <- dat5to6_IL2_clr$y

s2 <- proc.time() 
sl_grid_sem5to6_IL2 <- regmed.grid(x = x5to6_IL2_clean, mediator = medi5to6_IL2_clean, y = y5to6_IL2_clean, 
           lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(sl_grid_sem5to6_IL2, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Species level/sl_grid_sem5to6_IL2.Rdata")
times_sl[11,] <- e2[1:3]

#### Fractalkine
dat5to6_frac_clr <- mvregmed.dat.check(x = x_5to6, y = y5to6_frac, mediator = medi5to6_clr)
x5to6_frac_clean <- dat5to6_frac_clr$x; medi5to6_frac_clean <- dat5to6_frac_clr$mediator; y5to6_frac_clean <- dat5to6_frac_clr$y

s2 <- proc.time() 
sl_grid_sem5to6_frac <- regmed.grid(x = x5to6_frac_clean, mediator = medi5to6_frac_clean, y = y5to6_frac_clean, 
            lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(sl_grid_sem5to6_frac, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Species level/sl_grid_sem5to6_frac.Rdata")
times_sl[12,] <- e2[1:3]

#### Species-level 5to7 ####
#### Oral
micro5to7_oral <- microbiome_oral[which(microbiome_oral$Group %in% c(5,7)),]
micro5to7_oral_clo <- apply(micro5to7_oral[-(1:5)], 1, clo)
micro5to7_oral_clr <- apply(micro5to7_oral_clo, 1, clr)
colnames(micro5to7_oral_clr) <- paste0(colnames(micro5to7_oral_clr), "Oral")

#### Gut
micro5to7_gut <- microbiome_oral[which(microbiome_gut$Group %in% c(5,7)),]
micro5to7_gut_clo <- apply(micro5to7_gut[-(1:5)], 1, clo)
micro5to7_gut_clr <- apply(micro5to7_gut_clo, 1, clr)
colnames(micro5to7_gut_clr) <- paste0(colnames(micro5to7_gut_clr), "Gut")

#### Combined
medi5to7_clr <- cbind(micro5to7_oral_clr[,-which(colSums(micro5to7_oral_clr)==0)], 
                      micro5to7_gut_clr[,-which(colSums(micro5to7_gut_clr)==0)])

#### Treatment and Outcome
x_5to7 <- c(rep(1,8), rep(0,6)) # code grp 1 = 1, grp 2 = 0
y5to7_IL2 <- immune_5to7$`Serum_IL-2` ; y5to7_frac <- immune_5to7$Serum_Fractalkine


#### IL-2
dat5to7_IL2_clr <- mvregmed.dat.check(x = x_5to7, y = y5to7_IL2, mediator = medi5to7_clr)
x5to7_IL2_clean <- dat5to7_IL2_clr$x; medi5to7_IL2_clean <- dat5to7_IL2_clr$mediator; y5to7_IL2_clean <- dat5to7_IL2_clr$y

s2 <- proc.time() 
sl_grid_sem5to7_IL2 <- regmed.grid(x = x5to7_IL2_clean, mediator = medi5to7_IL2_clean, y = y5to7_IL2_clean, 
           lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(sl_grid_sem5to7_IL2, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Species level/sl_grid_sem5to7_IL2.Rdata")
times_sl[13,] <- e2[1:3]

#### Fractalkine
dat5to7_frac_clr <- mvregmed.dat.check(x = x_5to7, y = y5to7_frac, mediator = medi5to7_clr)
x5to7_frac_clean <- dat5to7_frac_clr$x; medi5to7_frac_clean <- dat5to7_frac_clr$mediator; y5to7_frac_clean <- dat5to7_frac_clr$y

s2 <- proc.time() 
sl_grid_sem5to7_frac <- regmed.grid(x = x5to7_frac_clean, mediator = medi5to7_frac_clean, y = y5to7_frac_clean, 
            lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(sl_grid_sem5to7_frac, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Species level/sl_grid_sem5to7_frac.Rdata")
times_sl[14,] <- e2[1:3]

#### Species-level 5to8 ####
#### Oral
micro5to8_oral <- microbiome_oral[which(microbiome_oral$Group %in% c(5,8)),]
micro5to8_oral_clo <- apply(micro5to8_oral[-(1:5)], 1, clo)
micro5to8_oral_clr <- apply(micro5to8_oral_clo, 1, clr)
colnames(micro5to8_oral_clr) <- paste0(colnames(micro5to8_oral_clr), "Oral")

#### Gut
micro5to8_gut <- microbiome_oral[which(microbiome_gut$Group %in% c(5,8)),]
micro5to8_gut_clo <- apply(micro5to8_gut[-(1:5)], 1, clo)
micro5to8_gut_clr <- apply(micro5to8_gut_clo, 1, clr)
colnames(micro5to8_gut_clr) <- paste0(colnames(micro5to8_gut_clr), "Gut")

#### Combined
medi5to8_clr <- cbind(micro5to8_oral_clr[,-which(colSums(micro5to8_oral_clr)==0)], 
                      micro5to8_gut_clr[,-which(colSums(micro5to8_gut_clr)==0)])

#### Treatment and Outcome
x_5to8 <- c(rep(1,8), rep(0,8)) # code grp 1 = 1, grp 2 = 0
y5to8_IL2 <- immune_5to8$`Serum_IL-2` ; y5to8_frac <- immune_5to8$Serum_Fractalkine


#### IL-2
dat5to8_IL2_clr <- mvregmed.dat.check(x = x_5to8, y = y5to8_IL2, mediator = medi5to8_clr)
x5to8_IL2_clean <- dat5to8_IL2_clr$x; medi5to8_IL2_clean <- dat5to8_IL2_clr$mediator; y5to8_IL2_clean <- dat5to8_IL2_clr$y

s2 <- proc.time() 
sl_grid_sem5to8_IL2 <- regmed.grid(x = x5to8_IL2_clean, mediator = medi5to8_IL2_clean, y = y5to8_IL2_clean, 
           lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(sl_grid_sem5to8_IL2, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Species level/sl_grid_sem5to8_IL2.Rdata")
times_sl[15,] <- e2[1:3]

#### Fractalkine
dat5to8_frac_clr <- mvregmed.dat.check(x = x_5to8, y = y5to8_frac, mediator = medi5to8_clr)
x5to8_frac_clean <- dat5to8_frac_clr$x; medi5to8_frac_clean <- dat5to8_frac_clr$mediator; y5to8_frac_clean <- dat5to8_frac_clr$y

s2 <- proc.time() 
sl_grid_sem5to8_frac <- regmed.grid(x = x5to8_frac_clean, mediator = medi5to8_frac_clean, y = y5to8_frac_clean, 
            lambda.grid, frac.lasso=0)
e2 <- proc.time() - s2

save(sl_grid_sem5to8_frac, file="~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Species level/sl_grid_sem5to8_frac.Rdata")
times_sl[16,] <- e2[1:3]


#### Save times ####
times_all <- rbind(times_gl, times_sl)
save(times_all, file = "~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Results/v6 results/Path times.Rdata")