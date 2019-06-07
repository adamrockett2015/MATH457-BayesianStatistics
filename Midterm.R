### MATH457: Midterm 2 ###

data = read.table("Y12.dat", header = TRUE)
attach(data)

plot(density(y), ylim = c(0, 0.025), main = "Histogram of Data", xlab = "y", ylab = "P(Y=y)")
hist(y, breaks = 40, add = TRUE, probability = TRUE)

# Set Prior Values
a = 1
b = 1
mu0 = 140
t20 = 625
sigma20 = 625
nu0 = 5
B = 100000
n = length(y)

THETA1 = numeric(B)
THETA2 = numeric(B)
SIGMA1 = numeric(B)
SIGMA2 = numeric(B)
Pr = numeric(B)
YPRED = numeric(B)


# initial values
P = 0.5
theta1 = 110
theta2 = 150
s21 = 100
s22 = 64

dist1 = 0
dist2 = 0

for (i in 1:B) {
  # Update X
  p1 = P*dnorm(y, theta1, sqrt(s21))
  p2 = (1-P)*dnorm(y, theta2, sqrt(s22))
  ptemp = p1 / (p1+p2)
  X = rbinom(n, 1, ptemp)
  
  # Get values that will be used later
  n1 = sum(X)
  n2 = n - n1
  y1 = y[X==1]
  y2 = y[X==0]
  ybar1 = mean(y1)
  ybar2 = mean(y2)
  yvar1 = var(y1)
  yvar2 = var(y2)
  
  # Update P
  P = rbeta(1, a+n1, b+n2)
  
  # Update theta1
  t2n1 = 1 / (1/t20 + n1/s21)
  mun1 = t2n1 * (1/t20 * mu0 + n1/s21 * ybar1)
  theta1 = rnorm(1, mun1, sqrt(t2n1))
  
  # Update theta2
  t2n2 = 1 / (1/t20 + n2/s22)
  mun2 = t2n2 * (1/t20 * mu0 + n2/s22 * ybar2)
  theta2 = rnorm(1, mun2, sqrt(t2n2))
  
  # Update sigma1
  nun1 = nu0 + n1
  s2n1 = 1/nun1 * (nu0*sigma20 + (n1-1)*yvar1 + n1*(ybar1-theta1)^2)
  s21 = 1 / rgamma(1, nun1/2, (nun1*s2n1)/2)
  
  # Update sigma2
  nun2 = nu0 + n2
  s2n2 = 1/nun2 * (nu0*sigma20 + (n2-1)*yvar2 + n2*(ybar2-theta2)^2)
  s22 = 1 / rgamma(1, nun2/2, (nun2*s2n2)/2)
  
  # Posterior Predictive
  xpred = runif(1) < P
  if (xpred) {
    ypred = rnorm(1, theta1, sqrt(s21))
    if (ypred > 119 && ypred < 121) {
      dist1 = dist1 + 1
    }
  } else {
    ypred = rnorm(1, theta2, sqrt(s22))
    if (ypred > 119 && ypred < 121) {
      dist2 = dist2 + 1
    }
  }
  
  # Save Values
  THETA1[i] = theta1
  THETA2[i] = theta2
  SIGMA1[i] = s21
  SIGMA2[i] = s22
  Pr[i] = P
  YPRED[i] = ypred
}

#Diagnostics
plot(density(THETA1), xlim = c(102,117),main = expression(paste("Density of ",theta[1])), xlab = expression(paste(theta[1])), ylab = expression(paste(p(theta[1]))), col = "black")
abline(v=quantile(THETA1, 0.025))
abline(v=quantile(THETA1, 0.975))
plot(density(THETA2), xlim = c(137,152), main = expression(paste("Density of ",theta[2])), xlab = expression(paste(theta[2])), ylab = expression(paste(p(theta[2]))), col = "black")
abline(v=quantile(THETA2, 0.025))
abline(v=quantile(THETA2, 0.975))


plot(density(y), ylim = c(0, 0.025), main = "Kernel Densities",lwd = 2, xlab = "y", ylab = "P(Y=y)", col = "black")
lines(density(YPRED), col = "black", lty = 2, lwd = 2)
legend(160, 0.0225, c("Original", "Mixture Model"), lty = c(1, 2), lwd = 2, cex = 1, bty = "n")

#
mean(Pr)
quantile(THETA1, c(0.025, 0.5, 0.975))
quantile(THETA2, c(0.025, 0.5, 0.975))
quantile(sqrt(SIGMA1), c(0.025, 0.5, 0.975))
quantile(sqrt(SIGMA2), c(0.025, 0.5, 0.975))
quantile(Pr, c(0.025, 0.5, 0.975))

# Autocorrelation
library(coda)
acf(THETA1)
effectiveSize(THETA1)

acf(THETA2)
effectiveSize(THETA2)

dist1 / (dist1 + dist2)

THETA16250 = THETA1
THETA26250 = THETA2

# Sensitivity
plot(density(THETA1), xlim = c(102,117),main = expression(theta[1] ~ tau[0]^2 ~ paste("=625")), xlab = expression(paste(theta[1])), ylab = expression(paste(p(theta[1]))), col = "black")
abline(v=quantile(THETA1, 0.025))
abline(v=quantile(THETA1, 0.975))
plot(density(THETA2), xlim = c(137,152), main = expression(theta[2] ~ tau[0]^2 ~ paste("=625")), xlab = expression(paste(theta[2])), ylab = expression(paste(p(theta[2]))), col = "black")
abline(v=quantile(THETA2, 0.025))
abline(v=quantile(THETA2, 0.975))
plot(density(THETA16250), xlim = c(102,117),main = expression(theta[1] ~ tau[0]^2 ~ paste("=6250")), xlab = expression(paste(theta[1])), ylab = expression(paste(p(theta[1]))), col = "black")
abline(v=quantile(THETA16250, 0.025))
abline(v=quantile(THETA16250, 0.975))
plot(density(THETA26250), xlim = c(137,152), main = expression(theta[2] ~ tau[0]^2 ~ paste("=6250")), xlab = expression(paste(theta[2])), ylab = expression(paste(p(theta[2]))), col = "black")
abline(v=quantile(THETA26250, 0.025))
abline(v=quantile(THETA26250, 0.975))







