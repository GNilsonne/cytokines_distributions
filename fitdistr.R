# Script to test distributions of IL-6 data

# By Gustav Nilsonne
# This is an example of how to run fitdistrplus to test distributions of IL-6 in one dataset

require(fitdistrplus)

Abhimanyu_derivation <- read.csv2("C:/Users/gusta/Box Sync/Gustavs_arbete/Pek/Cytokines should be log transformed/cytokines_distributions/Abhimanyu_derivation.csv")
Abhimanyu_derivation <- Abhimanyu_derivation[1:216, ] # Remove empty lines
Abhimanyu_validation <- read.csv2("C:/Users/gusta/Box Sync/Gustavs_arbete/Pek/Cytokines should be log transformed/cytokines_distributions/Abhimanyu_validation.csv")
Abhimanyu_followup <- read.csv2("C:/Users/gusta/Box Sync/Gustavs_arbete/Pek/Cytokines should be log transformed/cytokines_distributions/Abhimanyu_followup.csv")

IL6 <- data.frame(IL6 = Abhimanyu_derivation$IL.6[Abhimanyu_derivation$Class.Labels == "HC"], dataset = "Abhimanyu_derivation")
IL6$logIL6 <- log(IL6$IL6)
IL6$logIL6[IL6$logIL6 == -Inf] <- log(0.3) # Set 0 values to reasonable lower bound

test <- fitdist(IL6$IL6, distr = dnorm)
summary(test)
plot(test)

test2 <- fitdist(IL6$IL6, distr = dexp)
summary(test2)
plot(test2)

test3 <- fitdist(IL6$IL6, distr = dlnorm) # error

test4 <- fitdist(IL6$IL6, distr = dpois) # error


