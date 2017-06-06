# Script to analyse distributions of cytokine data

# 1. Read data

# Require packages
require(psych) 

# Define functions
skewd <- function(x, y, B = 1000){ 
  sx <- replicate(B, skew(sample(x, replace = TRUE))) 
  sy <- replicate(B, skew(sample(y, replace = TRUE))) 
  sx - sy  
} 

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

pdf("IL1b_1.pdf", height = 5, width = 5)
par(mar = c(0,0,0,0))
plot(d_IL1b, xaxs = "i", yaxs = "i", frame.plot= F, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "")
polygon(d_IL1b, col="gray", border = "gray")
plot(d_logIL1b, xaxs = "i", yaxs = "i", frame.plot= F, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "")
polygon(d_logIL1b, col="black", border = "black")

par(mar = c(5.1, 4.1, 4.1, 2.1))


dev.off()

# IL-6
d_IL6_1 <- density(IL6$IL6)
d_IL6_1$y[1] <- 0
d_IL6_1$y[length(d_IL6_1$y)] <- 0
d_logIL6_1 <- density(IL6$logIL6)
d_logIL6_1$y[1] <- 0
d_logIL6_1$y[length(d_IL6_1$y)] <- 0
shapiro_IL6_1 <- shapiro.test(IL6$IL6)
shapiro_logIL6_1 <- shapiro.test(IL6$logIL6)
qqIL6 <- qqnorm(IL6$IL6)
qqlogIL6 <- qqnorm(IL6$logIL6)

pdf("IL6_1.pdf", height = 5, width = 5)
par(mar = c(0,0,0,0))
plot(d_IL6_1, xaxs = "i", yaxs = "i", frame.plot= F, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "", ylim = c(0, max(d_IL6$y)))
polygon(d_IL6_1, col="black", border = "black")
plot(d_logIL6_1, xaxs = "i", yaxs = "i", frame.plot= F, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "")
polygon(d_logIL6_1, col="black", border = "black")

plot(qqIL6, frame.plot = F, xlab = "", ylab = "", main = "", type = "n", xaxt = "n", yaxt = "n")
y <- quantile(IL6$IL6, c(0.25, 0.75), names = FALSE, type = 7, na.rm = TRUE)
x <- qnorm(c(0.25, 0.75))
slope <- diff(y)/diff(x)
int <- y[1L] - slope * x[1L]
abline(int, slope, col = "gray", lwd = 3)
points(qqIL6, pch = 16, cex = 2)
legend(min(qqIL6$x)-1, max(qqIL6$y)+10, legend = paste("W =", round(shapiro_IL6_1$statistic, 2)), cex = 4, bty = "n")

plot(qqlogIL6, frame.plot = F, xlab = "", ylab = "", main = "", type = "n", xaxt = "n", yaxt = "n")
y <- quantile(IL6$logIL6, c(0.25, 0.75), names = FALSE, type = 7, na.rm = TRUE)
x <- qnorm(c(0.25, 0.75))
slope <- diff(y)/diff(x)
int <- y[1L] - slope * x[1L]
abline(int, slope, col = "gray", lwd = 3)
points(qqlogIL6, pch = 16, cex = 2)
legend(min(qqlogIL6$x)-0.5, max(qqlogIL6$y)+0.2, legend = paste("W =", round(shapiro_logIL6_1$statistic, 2)), cex = 4, bty = "n")
dev.off()

IL61boot <- skewd(IL6$logIL6, IL6$IL6, B = 10000) 
mean(IL61boot > 0) 
hist(IL61boot, breaks = 50) 





qqIL6 <- qqnorm(IL6$IL6)
plot(qqIL6, frame.plot = F, xlab = "Theoretical quantile", ylab = "Sample quantile", main = "IL-6, raw", type = "n")
y <- quantile(IL6$IL6, c(0.25, 0.75), names = FALSE, type = 7, na.rm = TRUE)
x <- qnorm(c(0.25, 0.75))
slope <- diff(y)/diff(x)
int <- y[1L] - slope * x[1L]
abline(int, slope, col = "gray")
points(qqIL6)


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

