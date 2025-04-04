# Script to analyse distributions of cytokine data

# Written by Love Ahnström

# Clear environment
rm(list=ls())

# 1. Read data

# Require packages
library(foreign)
library(tidyverse)
library(fitdistrplus)
library(haven)
library(moments)
library(readr)
library(readxl)
library(dplyr)

# Read data from files
Abhimanyu_df <- read.csv2("Abhimanyu_derivation.csv")
Sothern_df <- read.csv2("Sothern1995IndividualData.csv")
Imaeda_df <- read.csv2("Imaeda derivation.csv", na.strings = c("", "NA"), check.names = F, encoding = "UTF-8")
Imaeda_v_df <- read.csv2("Imaeda validation.csv", na.strings = c("", "NA"), check.names = F, encoding = "UTF-8")
Wand_df <- read.csv2("Wand (2016) data.csv")
MIDUS_df <- load("29282-0001-Data.rda")
IL6MIDUS <- data.frame(da29282.0001$B4BMSDIL6)
IL6MIDUS2 <- data.frame(da29282.0001$B4BIL6)
MIDJA_df <- load("34969-0001-Data.rda")
IL6MIDJA <- data.frame(da34969.0001$J2BIL6)
Lekander <- read.dta("psd_srh_il6.dta")
Lekander_df <- data.frame("psd_srh_il6.dta")
Lasselin <- read_sav("Lasselin.sav")
MIDJA2_df <- read_tsv("MIDJA2.tsv")
MIDUS3 <- read_sav("MIDUS3.sav")
load("MIDUS-REFRESHER.rda")
# There is also $RA4BMSDIL6 for another method of measurement. RA4BIL6 is ELISA.
IL6_MIDUS_REF <- data.frame(da36901.0001$RA4BIL6)
Yang <- read_xlsx("Yang.xlsx")
IL6_Yang <- Yang$`IL6(pg/mL)`

SEBAS1 <- read.csv("SEBAS-1.csv")
IL6_SEBAS1 <- rep(SEBAS1$val, times = SEBAS1$freq)
IL6_SEBAS1 <- data.frame(SEBAS1 = IL6_SEBAS1)
IL6_SEBAS1 <- IL6_SEBAS1 %>% filter(SEBAS1 != 0)

SEBAS2 <- read.csv("SEBAS-2.csv")
IL6_SEBAS2 <- rep(SEBAS2$val, times = SEBAS2$freq)
IL6_SEBAS2 <- data.frame(SEBAS2 = IL6_SEBAS2)
IL6_SEBAS2 <- IL6_SEBAS2 %>% filter(SEBAS2 != 0)

IL6_SEBAS_COMBINED <- bind_rows(IL6_SEBAS1, IL6_SEBAS2)
IL6_SEBAS_COMBINED <- data.frame(stack(IL6_SEBAS_COMBINED))
IL6_SEBAS_COMBINED <- data.frame(SEBAS_COMB = IL6_SEBAS_COMBINED$values)
IL6_SEBAS_COMBINED <- IL6_SEBAS_COMBINED %>% filter(SEBAS_COMB != 0)

# Abhimanyu - Extract IL6 and omit NAs
IL6Abhimanyu <- data.frame(Abhimanyu_df$IL.6)
IL6Abhimanyu[IL6Abhimanyu == 0] <- NA
IL6Abhimanyu <- na.omit(IL6Abhimanyu)

# Wand - extract IL-6
IL6Wand <- data.frame(Wand_df$IL6.d1.4h, Wand_df$IL6.d1.16h, Wand_df$IL6.d2.4h, Wand_df$IL6.d2.16h, Wand_df$IL6.d3.4h, Wand_df$IL6.d3.16h, Wand_df$IL6.d4.4h, Wand_df$IL6.d4.16h)
IL6Wand_stacked <- data.frame(stack(IL6Wand))[1]

# Imaeda - extract IL-6
IL6Imaeda <- data.frame(Imaeda_df$`Interleukin-6, pg/mL`)
IL6Imaeda_v <- data.frame(Imaeda_v_df$`Interleukin-6, pg/mL`)

# Sothern - re-arrange data into one column and extract IL6
Sothern_df <- data.frame(stack(Sothern_df))
IL6Sothern <- data.frame(Sothern_df$values)

# Lekander - extract IL6 and omit NAs
IL6Lekander <- data.frame(Lekander$il6m)
IL6Lekander<- na.omit(IL6Lekander)

# Lasselin - extract IL-6
IL6Lasselin <- data.frame(Lasselin$T0_IL6pgmL.LPS)

# MIDJA2 - extract IL-6
IL6MIDJA2 <- data.frame(MIDJA2_df$K2BIL6)
IL6MIDJA2 <- na.omit(IL6MIDJA2)

# MIDUS3 - extract IL-6
IL6MIDUS3 <- data.frame(as.numeric(MIDUS3$RA4BIL6))
IL6MIDUS3 <- na.omit(IL6MIDUS3)

# MIDUS REFRESHER - extract IL-6
IL6_MIDUS_REF <- data.frame(as.numeric(IL6_MIDUS_REF$da36901.0001.RA4BIL6))
IL6_MIDUS_REF <- na.omit(IL6_MIDUS_REF)

#Yang
IL6_Yang <- data.frame(as.numeric(IL6_Yang))
IL6_Yang <- na.omit(IL6_Yang)


# Stack data in a 2 column df 
IL6_only <- bind_rows(IL6Abhimanyu, IL6Wand_stacked, IL6Imaeda, IL6Imaeda_v, IL6Sothern, IL6Lekander, IL6MIDUS, IL6MIDUS2, IL6MIDJA, IL6Lasselin, IL6MIDJA2, IL6MIDUS3, IL6_MIDUS_REF, IL6_Yang, IL6_SEBAS1, IL6_SEBAS2, IL6_SEBAS_COMBINED)
colnames(IL6_only) <- c("Abhimanyu", "Wand", "Imaeda", "Imaeda_v", "Sothern", "Lekander", "IL6MIDUS", "IL6MIDUS2", "IL6MIDJA", "IL6Lasselin", "IL6MIDJA2", "IL6MIDUS3", "IL6_MIDUS_REF", "IL6_Yang", "IL6_SEBAS1", "IL6_SEBAS2", "IL6_SEBAS_COMBINED")
IL6_only_stacked <- data.frame(stack(IL6_only)) %>% filter(!is.na(values))

# Write results to file
write.csv(IL6_only_stacked, file="IL6-only-stacked.csv", row.names = F)

