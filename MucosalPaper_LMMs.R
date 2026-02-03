#### Load packages ####
library(lme4)
library(lmerTest)
library(readxl)
library(dplyr)

#### Read data ####
# read in the data
dentin <- read_xlsx(paste0("~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Data/",
                           "Surfaces.xlsx"), sheet = "Dentin")
enamel <- read_xlsx(paste0("~/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Data/",
                           "Surfaces.xlsx"), sheet = "Enamel")

# drop columns "Subject_ID_original" and "Study_Round"
dentin <- dentin[,-c(1,3)]; enamel <- enamel[,-c(1,3)]

# Subset groups for each comparison
dentin1345 <- dentin[which(dentin$Group_ID%in%c(1,3,4,5)),] %>% 
  mutate(Subject_ID = factor(Subject_ID),
         Group_ID = factor(Group_ID),
         Side = factor(Side))

dentin25678 <- dentin[which(dentin$Group_ID%in%c(2,5,6,7,8)),] %>% 
  mutate(Subject_ID = factor(Subject_ID),
         Group_ID = factor(Group_ID),
         Side = factor(Side))

enamel1345 <- enamel[which(enamel$Group_ID%in%c(1,3,4,5)),] %>% 
  mutate(Subject_ID = factor(Subject_ID),
         Group_ID = factor(Group_ID),
         Side = factor(Side))

enamel25678 <- enamel[which(enamel$Group_ID%in%c(2,5,6,7,8)),] %>% 
  mutate(Subject_ID = factor(Subject_ID),
         Group_ID = factor(Group_ID),
         Side = factor(Side))


#### Dentin models ####
dent1345_mod = lmer(Mandible_Dentin ~ Group_ID + (1|Side) + (1|Subject_ID), data = dentin1345,
                    control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-6)))
anova(dent1345_mod, ddf="Kenward-Roger")
ls_means(dent1345_mod, pairwise = T, ddf = "Kenward-Roger")

dent25678_mod = lmer(Mandible_Dentin ~ Group_ID + (1|Side) + (1|Subject_ID), data = dentin25678,
                     control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-6)))
anova(dent25678_mod, ddf="Kenward-Roger")
ls_means(dent25678_mod, pairwise = T, ddf = "Kenward-Roger") 

#### Enamel models ####
enam1345_mod = lmer(Mandible_Enamel ~ Group_ID + (1|Side) + (1|Subject_ID), data = enamel1345,
                    control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-6)))
anova(enam1345_mod, ddf="Kenward-Roger")
ls_means(enam1345_mod, pairwise = T, ddf = "Kenward-Roger") 

enamel25678_mod = lmer(Mandible_Enamel ~ Group_ID + (1|Side) + (1|Subject_ID), data = enamel25678,
                       control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-6)))

anova(enamel25678_mod, ddf="Kenward-Roger") %>% round(4) %>% knitr::kable(format = "latex")
ls_means(enamel25678_mod, pairwise = T, ddf = "Kenward-Roger") %>% round(4) %>% knitr::kable(format = "latex")


#### Mandible dentin model ####
dent1345_mod = lmer(Mandible_Dentin ~ Group_ID + (1|Side) + (1|Subject_ID), 
                    data = dentin1345,
                    control = lmerControl(check.conv.singular = .makeCC(action = "ignore", 
                                                                        tol = 1e-6)))
anova(dent1345_mod, ddf="Kenward-Roger")
ls_means(dent1345_mod, pairwise = T, ddf = "Kenward-Roger")
