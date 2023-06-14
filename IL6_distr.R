# Script to analyse distributions of cytokine data

# Written by Jimmy Hermansson
# Modified by Love Ahnström

# Clear environment
rm(list=ls())

# 1. Read data

# Require packages
require(foreign)
require(tidyverse)
require(fitdistrplus)
require(haven)
require(moments)

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

# Stack data in a 2 column df 
IL6_only <- bind_rows(IL6Abhimanyu, IL6Wand_stacked, IL6Imaeda, IL6Imaeda_v, IL6Sothern, IL6Lekander, IL6MIDUS, IL6MIDUS2, IL6MIDJA, IL6Lasselin, IL6MIDJA2)
colnames(IL6_only) <- c("Abhimanyu", "Wand", "Imaeda", "Imaeda_v", "Sothern", "Lekander", "IL6MIDUS", "IL6MIDUS2", "IL6MIDJA", "IL6Lasselin", "IL6MIDJA2")
IL6_only_stacked <- data.frame(stack(IL6_only)) %>% filter(!is.na(values))

# To simplify the syntax when writing files
setwd("../output")

# Group by dataset, for each group run the model and export csv file, once for dnorm and once for dexp
IL6_only_stacked %>% 
  group_by(ind) %>%
  group_walk(~{
    
    # For readability
    IL6_raw <- .x$values
    IL6_log <- log(IL6_raw)
    dataset_name <- .y$ind[1]
    
    # Fitdistr
    test_dnorm <- fitdist(IL6_raw, distr = dnorm)
    model_summary_dnorm <- summary(test_dnorm)
    
    # Scaled by a factor of 1/1000 due to model limitations
    test_dexp <- fitdist(IL6_raw/1000, distr = dexp)
    model_summary_dexp <- summary(test_dexp)
    
    # Shapiro-wilk normality test
    shapiro_test <- shapiro.test(IL6_raw)
    shapiro_test_log <- shapiro.test(IL6_log)
    
    # Skewness
    skewness_norm <- skewness(IL6_raw)
    skewness_log <- skewness(IL6_log)
    
    #group_walk silences output to console, going around this by writing csv files instead
    df <- data.frame(
      dataset=.y$ind[1], 
      n=model_summary_dnorm$n, 
      w_raw = round(shapiro_test$statistic, digits=3), 
      p_raw = signif(shapiro_test$p.value, digits=4), 
      w_log = round(shapiro_test_log$statistic, digits=3) , 
      p_log = signif(shapiro_test_log$p.value, digits=4), 
      loglik_norm=round(model_summary_dnorm$loglik, digits=1), 
      aic_norm=round(model_summary_dnorm$aic, digits=1), 
      bic_norm=round(model_summary_dnorm$bic, digits=1), 
      loglik_exp=round(model_summary_dexp$loglik, digits=1), 
      aic_exp=round(model_summary_dexp$aic, digits=1), 
      bic_exp=round(model_summary_dexp$bic, digits=1),
      skew_norm=round(skewness_norm, digits=3),
      skew_log=round(skewness_log, digits=3)
      )
    write.csv(df, sprintf("%s-fitdistr.csv", .y$ind[1]), row.names = F)
    
    # Create plots
    histIL6 <- hist(IL6_raw)
    histlogIL6 <- hist(IL6_log)
    qqIL6 <- qqnorm(IL6_raw)
    qqlogIL6 <- qqnorm(IL6_log)
    
    # Write to pdf
    setwd("../pdf")
    
    pdf(sprintf("%s.pdf", dataset_name), width = 5, height = 5) 
    
    par(mfrow=c(2,2), mar=c(3,3,3,3))
    plot(histIL6, main = "Raw data", xlab = "IL-6 (pg/µl)") 
    plot(histlogIL6, main = "Log-transformed", xlab = "IL-6 (pg/µl)") 
    plot(qqIL6, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", frame.plot = F) 
    abline(0,1)
    plot(qqlogIL6,  xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", frame.plot = F) 
    abline(0,1)
    
    dev.off()
    
    # Write to image
    setwd("../png")
    
    png(sprintf("%s.png", dataset_name), width = 5, height = 5, units="in", res=150) 
    
    par(mfrow=c(2,2), mar=c(3,3,3,3))
    plot(histIL6, main = "Raw data", xlab = "IL-6 (pg/µl)") 
    plot(histlogIL6, main = "Log-transformed", xlab = "IL-6 (pg/µl)") 
    plot(qqIL6, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", frame.plot = F) 
    abline(0,1)
    plot(qqlogIL6,  xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", frame.plot = F) 
    abline(0,1)
    
    dev.off()
    
    
    setwd("../output")
    
  })


#Read and merge the CSV files again
outputs <- list.files(path='.') %>% 
  lapply(read_csv) %>% 
  bind_rows 
write.csv(outputs, "../all_outputs.csv", row.names = F)

setwd("..")
png("sample-size-norm.png", width = 5, height = 5, units="in", res=150) 
plot(n ~ skew_norm, data = outputs, main = "Raw data", xlab="Skewness")
abline(v=1)
dev.off()

png("sample-size-log.png", width = 5, height = 5, units="in", res=150) 
plot(n ~ skew_log, data = outputs, main = "Log-transformed", xlab="Skewness")
abline(v=1)
dev.off()

