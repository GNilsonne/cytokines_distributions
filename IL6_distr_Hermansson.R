# Script to analyse distributions of cytokine data

# Written by Jimmy Hermansson
# Appended by Love Ahnström

# Clear environment
rm(list=ls())

# 1. Read data

# Require packages
require(foreign)
require(tidyverse)
require(fitdistrplus)

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

# Stack data in a 2 column df 
IL6_only <- bind_rows(IL6Abhimanyu, IL6Wand_stacked, IL6Imaeda, IL6Imaeda_v, IL6Sothern, IL6Lekander, IL6MIDUS, IL6MIDUS2, IL6MIDJA)
colnames(IL6_only) <- c("Abhimanyu", "Wand", "Imaeda", "Imaeda_v", "Sothern", "Lekander", "IL6MIDUS", "IL6MIDUS2", "IL6MIDJA")
IL6_only_stacked <- data.frame(stack(IL6_only)) %>% filter(!is.na(values))

# To simplify the syntax when writing files
setwd("../output")

# Group by dataset, for each group run the model and export csv file, once for dnorm and once for dexp
IL6_only_stacked %>% 
  group_by(ind) %>%
  group_walk(~{
    
    # Fitdistr
    test_dnorm <- fitdist(.x$values, distr = dnorm)
    model_summary_dnorm <- summary(test_dnorm)
    
    test_dexp <- fitdist(.x$values/1000, distr = dexp)
    model_summary_dexp <- summary(test_dexp)
    
    # Shapiro-wilk normality test
    for_shapiro <- .x$values
    for_shapiro_log <- log(.x$values)
    shapiro_test <- shapiro.test(for_shapiro)
    shapiro_test_log <- shapiro.test(for_shapiro_log)
    
    
    
    #group_walk silences output to console, going around this by writing csv files instead
    df <- data.frame(
      dataset=.y$ind[1], 
      n=model_summary_dnorm$n, 
      w_raw = shapiro_test$statistic, 
      p_raw = shapiro_test$p.value, 
      w_log = shapiro_test_log$statistic, 
      p_log = shapiro_test_log$p.value, 
      loglike_norm=model_summary_dnorm$loglik, 
      aic_norm=model_summary_dnorm$aic, 
      bic_norm=model_summary_dnorm$bic, 
      loglike_exp=model_summary_dexp$loglik, 
      aic_exp=model_summary_dexp$aic, 
      bic_exp=model_summary_dexp$bic)
    write.csv(df, sprintf("%s-fitdistr.csv", .y$ind[1]), row.names = F)
    
    
    # plot(test_dnorm)
    # plot(test_dexp)
    
  })


#Read and merge the CSV files again
outputs <- list.files(path='.') %>% 
  lapply(read_csv) %>% 
  bind_rows 
write.csv(outputs, "../all_outputs.csv", row.names = F)

# 2. Analyse distributions and create tables


# Wand analysis and tables

# Count nr of rows for later
Wandn <- nrow(IL6Wand)
# Make values numerical
IL6Wand <- as.numeric(IL6Wand$Wand_df.IL6.d1.4h, IL6Wand$Wand_df.IL6.d1.16, IL6Wand$Wand_df.IL6.d2.4h, IL6Wand$Wand_df.IL6.d2.16h, IL6Wand$Wand_df.IL6.d3.4h, IL6Wand$Wand_df.IL6.d3.16h, IL6Wand$Wand_df.IL6.d4.4h, IL6Wand$Wand_df.IL6.d4.16h)
# Assign data to matrix and convert to log
IL6Wand = data.matrix(IL6Wand)
logIL6Wand <- log(IL6Wand)
# Plot histograms
histIL6wand <- hist(IL6Wand, main = "Raw data", xlab = "IL-6 (pg/?l)")
histlogIL6wand <- hist(logIL6Wand, main = "Log-transformed", xlab = "IL-6 (pg/?l)")
# Normal QQ-plot
qqIL6Wand <- qqnorm(IL6Wand)
# Log QQ-plot
qqlogIL6Wand <- qqnorm(logIL6Wand)
# Write plots
pdf("histIL6wand.pdf", width = 3, height = 3) 
plot(histIL6wand, main = "Raw data", xlab = "IL-6 (pg/?l)") 
dev.off()
pdf("histlogIL6wand.pdf", width = 3, height = 3) 
plot(histlogIL6wand, main = "Log-transformed", xlab = "IL-6 (pg/?l)") 
dev.off()
pdf("qqIL6wand.pdf", width = 3, height = 3) 
plot(qqIL6Wand, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles") 
dev.off() 
pdf("qqlogIL6wand.pdf", width = 3, height = 3) 
plot(qqlogIL6Wand, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles") 
dev.off()
# Shapiro-wilk normality test
WandSW <- shapiro.test(IL6Wand)
logWandSW <- shapiro.test(logIL6Wand)
# Create dataframes for W statistic and p value
WandSWw <- as.data.frame(WandSW$statistic)
WandSWp <- as.data.frame(WandSW$p.value)
# Merge W and p statistic into table
Wandtable <- merge(WandSWp, WandSWw)
# Create dataframes for W statistic and p value - log
logWandSWw <- as.data.frame(logWandSW$statistic)
logWandSWp <- as.data.frame(logWandSW$p.value)
# Merge w and p statistic into table - log
Wandtablelog <- merge(logWandSWp, logWandSWw)
# Merge norm and log tables and add nr of obs.
Wandtable <- merge(Wandtable, Wandtablelog)
Wandtable <- cbind(Wandn, Wandtable)
# Table names
names(Wandtable) <- c("n", "p (norm)", "W (norm)", "p (log)", "W (log)")
Wandtable$Study <- "Wand"
Wandtable <- Wandtable[, c("Study", "n", "W (norm)", "p (norm)", "W (log)", "p (log)")]


# Imaeda analysis and tables

Imaedan <- nrow(IL6Imaeda)

IL6Imaeda <- as.numeric(IL6Imaeda$Imaeda_df.Interleukin.6..pg.mL)
logIL6Imaeda <- log(IL6Imaeda)

histIL6Imaeda <- hist(IL6Imaeda, main = "Raw data", xlab = "IL-6 (pg/ml)")
histlogIL6Imaeda <- hist(logIL6Imaeda, main = "Log-transformed", xlab = "IL-6 (pg/ml)")

qqIL6Imaeda <- qqnorm(IL6Imaeda)
qqlogIL6Imaeda <- qqnorm(logIL6Imaeda)

pdf("histIL6Imaeda.pdf", width = 3, height = 3) 
plot(histIL6Imaeda, main = "Raw data", xlab = "IL-6 (pg/ml)") 
dev.off()
pdf("histlogIL6Imaeda.pdf", width = 3, height = 3) 
plot(histlogIL6Imaeda, main = "Log-transformed", xlab = "IL-6 (pg/ml)") 
dev.off()
pdf("qqIL6Imaeda.pdf", width = 3, height = 3) 
plot(qqIL6Imaeda, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
dev.off() 
pdf("qqlogIL6Imaeda.pdf", width = 3, height = 3) 
plot(qqlogIL6Imaeda, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
dev.off()

ImaedaSW <- shapiro.test(IL6Imaeda)
logImaedaSW <- shapiro.test(logIL6Imaeda)

ImaedaSWw <- as.data.frame(ImaedaSW$statistic)
ImaedaSWp <- as.data.frame(ImaedaSW$p.value)

Imaedatable <- merge(ImaedaSWp, ImaedaSWw)

logImaedaSWw <- as.data.frame(logImaedaSW$statistic)
logImaedaSWp <- as.data.frame(logImaedaSW$p.value)

Imaedatablelog <- merge(logImaedaSWp, logImaedaSWw)

Imaedatable <- merge(Imaedatable, Imaedatablelog)
Imaedatable <- cbind(Imaedan, Imaedatable)
names(Imaedatable) <- c("n", "p (norm)", "W (norm)", "p (log)", "W (log)")
Imaedatable$Study <- "Imaeda"
Imaedatable <- Imaedatable[, c("Study", "n", "W (norm)", "p (norm)", "W (log)", "p (log)")]

Imaeda_vn <- nrow(IL6Imaeda_v)

IL6Imaeda_v <- as.numeric(IL6Imaeda_v$Imaeda_v_df.Interleukin.6..pg.mL)
logIL6Imaeda_v <- log(IL6Imaeda_v)

histIL6Imaeda_v <- hist(IL6Imaeda_v, main = "Raw data", xlab = "IL-6 (pg/ml)")
histlogIL6Imaeda_v <- hist(logIL6Imaeda_v, main = "Log-transformed", xlab = "IL-6 (pg/ml)")

qqIL6Imaeda_v <- qqnorm(IL6Imaeda_v)
qqlogIL6Imaeda_v <- qqnorm(logIL6Imaeda_v)

Imaeda_vSW <- shapiro.test(IL6Imaeda_v)
logImaeda_vSW <- shapiro.test(logIL6Imaeda_v)

pdf("histIL6Imaeda_v.pdf", width = 3, height = 3) 
plot(histIL6Imaeda_v, main = "Raw data", xlab = "IL-6 (pg/ml)") 
dev.off()
pdf("histlogIL6Imaeda_v.pdf", width = 3, height = 3) 
plot(histlogIL6Imaeda_v, main = "Log-transformed", xlab = "IL-6 (pg/ml)") 
dev.off()
pdf("qqIL6Imaeda_v.pdf", width = 3, height = 3) 
plot(qqIL6Imaeda_v, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
dev.off() 
pdf("qqlogIL6Imaeda_v.pdf", width = 3, height = 3) 
plot(qqlogIL6Imaeda_v, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
dev.off()

Imaeda_vSWw <- as.data.frame(Imaeda_vSW$statistic)
Imaeda_vSWp <- as.data.frame(Imaeda_vSW$p.value)

Imaeda_vtable <- merge(Imaeda_vSWp, Imaeda_vSWw)

logImaeda_vSWw <- as.data.frame(logImaeda_vSW$statistic)
logImaeda_vSWp <- as.data.frame(logImaeda_vSW$p.value)

Imaeda_vtablelog <- merge(logImaeda_vSWp, logImaeda_vSWw)

Imaeda_vtable <- merge(Imaeda_vtable, Imaeda_vtablelog)
Imaeda_vtable <- cbind(Imaeda_vn, Imaeda_vtable)
names(Imaeda_vtable) <- c("n", "p (norm)", "W (norm)", "p (log)", "W (log)")
Imaeda_vtable$Study <- "Imaeda (Validation)"
Imaeda_vtable <- Imaeda_vtable[, c("Study", "n", "W (norm)", "p (norm)", "W (log)", "p (log)")]


# MIDUS analysis and tables (note that MIDUS uses MSD)

MIDUSn <- nrow(IL6MIDUS)

IL6MIDUS <- as.numeric(IL6MIDUS$da29282.0001.B4BMSDIL6)
logIL6MIDUS <- log(IL6MIDUS)

histIL6MIDUS <- histIL6MIDUS <- hist(IL6MIDUS, main = "Raw data", xlab = "IL-6 (pg/ml)")
histlogIL6MIDUS <- hist(logIL6MIDUS, main = "Log-transformed", xlab = "IL-6 (pg/ml)")

qqIL6MIDUS <- qqnorm(IL6MIDUS)
qqlogIL6MIDUS <- qqnorm(logIL6MIDUS)

MIDUSSW <- shapiro.test(IL6MIDUS)
logMIDUSSW <- shapiro.test(logIL6MIDUS)

pdf("histIL6MIDUS.pdf", width = 3, height = 3) 
plot(histIL6MIDUS, main = "Raw data", xlab = "IL-6 (pg/ml)") 
dev.off()
pdf("histlogIL6MIDUS.pdf", width = 3, height = 3) 
plot(histlogIL6MIDUS, main = "Log-transformed", xlab = "IL-6 (pg/ml)") 
dev.off()
pdf("qqIL6MIDUS.pdf", width = 3, height = 3) 
plot(qqIL6MIDUS, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
dev.off() 
pdf("qqlogIL6MIDUS.pdf", width = 3, height = 3) 
plot(qqlogIL6MIDUS, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
dev.off()

MIDUSSWw <- as.data.frame(MIDUSSW$statistic)
MIDUSSWp <- as.data.frame(MIDUSSW$p.value)

MIDUStable <- merge(MIDUSSWp, MIDUSSWw)

logMIDUSSWw <- as.data.frame(logMIDUSSW$statistic)
logMIDUSSWp <- as.data.frame(logMIDUSSW$p.value)

MIDUStablelog <- merge(logMIDUSSWp, logMIDUSSWw)

MIDUStable <- merge(MIDUStable, MIDUStablelog)
MIDUStable <- cbind(MIDUSn, MIDUStable)

names(MIDUStable) <- c("n", "p (norm)", "W (norm)", "p (log)", "W (log)")
MIDUStable$Study <- "MIDUS (MSD)"
MIDUStable <- MIDUStable[, c("Study", "n", "W (norm)", "p (norm)", "W (log)", "p (log)")]

MIDUS2n <- nrow(IL6MIDUS2)

IL6MIDUS2 <- as.numeric(IL6MIDUS2$da29282.0001.B4BIL6)
logIL6MIDUS2 <- log(IL6MIDUS2)

histIL6MIDUS2 <- hist(IL6MIDUS2, main = "Raw data", xlab = "IL-6 (pg/ml)")
histlogIL6MIDUS2 <- hist(logIL6MIDUS2, main = "Log-transformed", xlab = "IL-6 (pg/ml)")

qqIL6MIDUS2 <- qqnorm(IL6MIDUS2)
qqlogIL6MIDUS2 <- qqnorm(logIL6MIDUS2)

pdf("histIL6MIDUS2.pdf", width = 3, height = 3) 
plot(histIL6MIDUS2, main = "Raw data", xlab = "IL-6 (pg/ml)") 
dev.off()
pdf("histlogIL6MIDUS2.pdf", width = 3, height = 3) 
plot(histlogIL6MIDUS2, main = "Log-transformed", xlab = "IL-6 (pg/ml)") 
dev.off()
pdf("qqIL6MIDUS2.pdf", width = 3, height = 3) 
plot(qqIL6MIDUS2, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
dev.off() 
pdf("qqlogIL6MIDUS2.pdf", width = 3, height = 3) 
plot(qqlogIL6MIDUS2, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
dev.off()

MIDUS2SW <- shapiro.test(IL6MIDUS2)
logMIDUS2SW <- shapiro.test(logIL6MIDUS2)

MIDUS2SWw <- as.data.frame(MIDUS2SW$statistic)
MIDUS2SWp <- as.data.frame(MIDUS2SW$p.value)

MIDUS2table <- merge(MIDUS2SWp, MIDUS2SWw)

logMIDUS2SWw <- as.data.frame(logMIDUS2SW$statistic)
logMIDUS2SWp <- as.data.frame(logMIDUS2SW$p.value)

MIDUS2tablelog <- merge(logMIDUS2SWp, logMIDUS2SWw)

MIDUS2table <- merge(MIDUS2table, MIDUS2tablelog)

MIDUS2table <- cbind(MIDUS2n, MIDUS2table)

names(MIDUS2table) <- c("n", "p (norm)", "W (norm)", "p (log)", "W (log)")
MIDUS2table$Study <- "MIDUS"
MIDUS2table <- MIDUS2table[, c("Study", "n", "W (norm)", "p (norm)", "W (log)", "p (log)")]


# MIDJA analysis and tables

MIDJAn <- nrow(IL6MIDJA)

IL6MIDJA <- as.numeric(da34969.0001$J2BIL6)
logIL6MIDJA <- log(IL6MIDJA)

histIL6MIDJA <- hist(IL6MIDJA, main = "Raw data", xlab = "IL-6 (pg/ml)")
histlogIL6MIDJA <- hist(logIL6MIDJA, main = "Log-transformed", xlab = "IL-6 (pg/ml)")

qqIL6MIDJA <- qqnorm(IL6MIDJA)
qqlogIL6MIDJA <- qqnorm(logIL6MIDJA)

MIDJASW <- shapiro.test(IL6MIDJA)
logMIDJASW <- shapiro.test(logIL6MIDJA)

pdf("histIL6MIDJA.pdf", width = 3, height = 3) 
plot(histIL6MIDJA, main = "Raw data", xlab = "IL-6 (pg/ml)") 
dev.off()
pdf("histlogIL6MIDJA.pdf", width = 3, height = 3) 
plot(histlogIL6MIDJA, main = "Log-transformed", xlab = "IL-6 (pg/ml)") 
dev.off()
pdf("qqIL6MIDJA.pdf", width = 3, height = 3) 
plot(qqIL6MIDJA, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
dev.off() 
pdf("qqlogIL6MIDJA.pdf", width = 3, height = 3) 
plot(qqlogIL6MIDJA, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
dev.off()

MIDJASWw <- as.data.frame(MIDJASW$statistic)
MIDJASWp <- as.data.frame(MIDJASW$p.value)

MIDJAtable <- merge(MIDJASWp, MIDJASWw)

logMIDJASWw <- as.data.frame(logMIDJASW$statistic)
logMIDJASWp <- as.data.frame(logMIDJASW$p.value)

MIDJAtablelog <- merge(logMIDJASWp, logMIDJASWw)

MIDJAtable <- merge(MIDJAtable, MIDJAtablelog)

MIDJAtable <- cbind(MIDJAn, MIDJAtable)
names(MIDJAtable) <- c("n", "p (norm)", "W (norm)", "p (log)", "W (log)")
MIDJAtable$Study <- "MIDJA"
MIDJAtable <- MIDJAtable[, c("Study", "n", "W (norm)", "p (norm)", "W (log)", "p (log)")]

# Sothern analysis and tables

Sothernn <- nrow(IL6Sothern)

IL6Sothern <- as.numeric(IL6Sothern$Sothern_df.values)
logIL6Sothern <- log(IL6Sothern)

histIL6Sothern <- hist(IL6Sothern, main = "Raw data", xlab = "IL-6 (pg/ml)")
histlogIL6Sothern <- hist(logIL6Sothern, main = "Log-transformed", xlab = "IL-6 (pg/ml)")

qqIL6Sothern <- qqnorm(IL6Sothern)
qqlogIL6Sothern <- qqnorm(logIL6Sothern)

pdf("histIL6Sothern.pdf", width = 3, height = 3) 
plot(histIL6Sothern, main = "Raw data", xlab = "IL-6 (pg/ml)") 
dev.off()
pdf("histlogIL6Sothern.pdf", width = 3, height = 3) 
plot(histlogIL6Sothern, main = "Log-transformed", xlab = "IL-6 (pg/ml)") 
dev.off()
pdf("qqIL6Sothern.pdf", width = 3, height = 3) 
plot(qqIL6Sothern, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
dev.off() 
pdf("qqlogIL6Sothern.pdf", width = 3, height = 3) 
plot(qqlogIL6Sothern, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
dev.off()

SothernSW <- shapiro.test(IL6Sothern)
logIL6SothernSW <- shapiro.test(logIL6Sothern)

SothernSWw <- as.data.frame(SothernSW$statistic)
SothernSWp <- as.data.frame(SothernSW$p.value)

Sotherntable <- merge(SothernSWp, SothernSWw)

logIL6SothernSWw <- as.data.frame(logIL6SothernSW$statistic)
logIL6SothernSWp <- as.data.frame(logIL6SothernSW$p.value)

Sotherntablelog <- merge(logIL6SothernSWp, logIL6SothernSWw)

Sotherntable <- merge(Sotherntable, Sotherntablelog)
Sotherntable <- cbind(Sothernn, Sotherntable)
names(Sotherntable) <- c("n", "p (norm)", "W (norm)", "p (log)", "W (log)")
Sotherntable$Study <- "Sothern"
Sotherntable <- Sotherntable[, c("Study", "n", "W (norm)", "p (norm)", "W (log)", "p (log)")]

# Abhimanyu analysis and tables

Abhimanyun <- nrow(IL6Abhimanyu)

IL6Abhimanyu <- as.numeric(IL6Abhimanyu$Abhimanyu_df.IL.6)
logIL6Abhimanyu <- log(IL6Abhimanyu)

histIL6Abhimanyu <- hist(IL6Abhimanyu, main = "Raw data", xlab = "IL-6 (pg/ml)")
histlogIL6Abhimanyu <- hist(logIL6Abhimanyu, main = "Log-transformed", xlab = "IL-6 (pg/ml)")

qqIL6Abhimanyu <- qqnorm(IL6Abhimanyu)
qqlogIL6Abhimanyu <- qqnorm(logIL6Abhimanyu)

pdf("histIL6Abhimanyu.pdf", width = 3, height = 3) 
plot(histIL6Abhimanyu, main = "Raw data", xlab = "IL-6 (pg/ml)") 
dev.off()
pdf("histlogIL6Abhimanyu.pdf", width = 3, height = 3) 
plot(histlogIL6Abhimanyu, main = "Log-transformed", xlab = "IL-6 (pg/ml)") 
dev.off()
pdf("qqIL6Abhimanyu.pdf", width = 3, height = 3) 
plot(qqIL6Abhimanyu, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
dev.off() 
pdf("qqlogIL6Abhimanyu.pdf", width = 3, height = 3) 
plot(qqlogIL6Abhimanyu, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
dev.off()

AbhimanyuSW <- shapiro.test(IL6Abhimanyu)
logAbhimanyuSW <- shapiro.test(logIL6Abhimanyu)

AbhimanyuSWw <- as.data.frame(AbhimanyuSW$statistic)
AbhimanyuSWp <- as.data.frame(AbhimanyuSW$p.value)

Abhimanyutable <- merge(AbhimanyuSWp, AbhimanyuSWw)

logAbhimanyuSWw <- as.data.frame(logAbhimanyuSW$statistic)
logAbhimanyuSWp <- as.data.frame(logAbhimanyuSW$p.value)

Abhimanyutablelog <- merge(logAbhimanyuSWp, logAbhimanyuSWw)

Abhimanyutable <- merge(Abhimanyutable, Abhimanyutablelog)
Abhimanyutable <- cbind(Abhimanyun, Abhimanyutable)

names(Abhimanyutable) <- c("n", "p (norm)", "W (norm)", "p (log)", "W (log)")
Abhimanyutable$Study <- "Abhimanyu"
Abhimanyutable <- Abhimanyutable[, c("Study", "n", "W (norm)", "p (norm)", "W (log)", "p (log)")]

# Lekander analysis and tables

Lekandern <- nrow(IL6Lekander)

IL6Lekander <- as.numeric(IL6Lekander$Lekander.il6m)
logIL6Lekander <- log(IL6Lekander)

histIL6Lekander <- hist(IL6Lekander, main = "Raw data", xlab = "IL-6 (pg/ml)")
histlogIL6Lekander <- hist(logIL6Lekander, main = "Log-transformed", xlab = "IL-6 (pg/ml)")

qqIL6Lekander <- qqnorm(IL6Lekander)
qqlogIL6Lekander <- qqnorm(logIL6Lekander)

pdf("histIL6Lekander.pdf", width = 3, height = 3) 
plot(histIL6Lekander, main = "Raw data", xlab = "IL-6 (pg/ml)") 
dev.off()
pdf("histlogIL6Lekander.pdf", width = 3, height = 3) 
plot(histlogIL6Lekander, main = "Log-transformed", xlab = "IL-6 (pg/ml)") 
dev.off()
pdf("qqIL6Lekander.pdf", width = 3, height = 3) 
plot(qqIL6Lekander, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
dev.off() 
pdf("qqlogIL6Lekander.pdf", width = 3, height = 3) 
plot(qqlogIL6Lekander, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
dev.off()

LekanderSW <- shapiro.test(IL6Lekander)
logLekanderSW <- shapiro.test(logIL6Lekander)

LekanderSWw <- as.data.frame(LekanderSW$statistic)
LekanderSWp <- as.data.frame(LekanderSW$p.value)

Lekandertable <- merge(LekanderSWp, LekanderSWw)

logLekanderSWw <- as.data.frame(logLekanderSW$statistic)
logLekanderSWp <- as.data.frame(logLekanderSW$p.value)

Lekandertablelog <- merge(logLekanderSWp, logLekanderSWw)

Lekandertable <- merge(Lekandertable, Lekandertablelog)
Lekandertable <- cbind(Lekandern, Lekandertable)

names(Lekandertable) <- c("n", "p (norm)", "W (norm)", "p (log)", "W (log)")
Lekandertable$Study <- "Lekander"
Lekandertable <- Lekandertable[, c("Study", "n", "W (norm)", "p (norm)", "W (log)", "p (log)")]

# 3. Merge datasets and write table

goftable <- rbind(Abhimanyutable, Imaedatable, Imaeda_vtable, Lekandertable, MIDJAtable, MIDUStable, MIDUS2table, Sotherntable, Wandtable)

write.csv2(goftable, "Tabell GOF.csv")
