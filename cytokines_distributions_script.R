# Script to analyse distributions of cytokine data

# 1. Read data

# Read data from files
Abhimanyu_derivation <- read.csv2("C:/Users/Gustav Nilsonne/Box Sync/Gustavs_arbete/Pek/Cytokines should be log transformed/cytokines_distributions/Abhimanyu_derivation.csv")
Abhimanyu_derivation <- Abhimanyu_derivation[1:216, ] # Remove empty lines
Abhimanyu_validation <- read.csv2("C:/Users/Gustav Nilsonne/Box Sync/Gustavs_arbete/Pek/Cytokines should be log transformed/cytokines_distributions/Abhimanyu_validation.csv")
Abhimanyu_followup <- read.csv2("C:/Users/Gustav Nilsonne/Box Sync/Gustavs_arbete/Pek/Cytokines should be log transformed/cytokines_distributions/Abhimanyu_followup.csv")


# Make data frames for each cytokine
IL1b <- data.frame(IL1b = Abhimanyu_derivation$IL.1b[Abhimanyu_derivation$Class.Labels == "HC"], dataset = "Abhimanyu_derivation")
IL1b$logIL1b <- log(IL1b$IL1b)

IL6 <- data.frame(IL6 = Abhimanyu_derivation$IL.6[Abhimanyu_derivation$Class.Labels == "HC"], dataset = "Abhimanyu_derivation")
IL6$logIL6 <- log(IL6$IL6)
IL6$logIL6[IL6$logIL6 == -Inf] <- log(0.3) # Set 0 values to reasonable lower bound

# 2. Analyse distributions

# IL-1b
hist(IL1b$IL1b)
hist(log(IL1b$IL1b))

d_IL1b <- density(IL1b$IL1b)
d_logIL1b <- density(IL1b$logIL1b)
plot(d_IL1b, xaxs = "i", yaxs = "i", frame.plot= F, type = "n", ylim = c(0, 1))
polygon(d_IL1b, col="gray", border = "gray")
lines(d_logIL1b, lwd = 2)
shapiro.test(IL1b$IL1b)
shapiro.test(IL1b$logIL1b)




# IL-6
hist(IL6$IL6)
hist(IL6$logIL6))
d_IL6 <- density(IL6$IL6)
d_logIL6 <- density(IL6$logIL6)
plot(d_IL6, xlim = c(-2, 3), xaxs = "i", yaxs = "i", frame.plot= F, type = "n")
polygon(d_IL6, col="gray", border = "gray")
lines(d_logIL6, lwd = 2)
shapiro.test(IL6$IL6)
shapiro.test(IL6$logIL6)

qqIL6 <- qqnorm(IL6$IL6)
plot(qqIL6, frame.plot = F, xlab = "Theoretical quantile", ylab = "Sample quantile", main = "IL-6, raw", type = "n")
y <- quantile(IL6$IL6, c(0.25, 0.75), names = FALSE, type = 7, na.rm = TRUE)
x <- qnorm(c(0.25, 0.75))
slope <- diff(y)/diff(x)
int <- y[1L] - slope * x[1L]
abline(int, slope, col = "gray")
points(qqIL6)

qqlogIL6 <- qqnorm(IL6$logIL6)
plot(qqlogIL6, frame.plot = F, xlab = "Theoretical quantile", ylab = "Sample quantile", main = "IL-6, log-transformed", type = "n")
y <- quantile(IL6$logIL6, c(0.25, 0.75), names = FALSE, type = 7, na.rm = TRUE)
x <- qnorm(c(0.25, 0.75))
slope <- diff(y)/diff(x)
int <- y[1L] - slope * x[1L]
abline(int, slope, col = "gray")
points(qqlogIL6)




require(lattice)
myDist<-function(x) {
  qexp(x, 5)
}

set.seed(15)
x <- rexp(100, 5)
qqmath(~x, distribution=myDist, main="qqmath")


exp.x <- myDist(ppoints(length(x)))
xyplot(sort(x)~exp.x, main="xyplot")


plot(sort(IL6$IL6) ~ sort(rnorm(mean = 0, n = length(IL6$IL6))))


plot(quantile(IL6$IL6, probs = c(1:100)/100) ~ qnorm(p = c(1:100)/100))

ppoints(100)

mean()

