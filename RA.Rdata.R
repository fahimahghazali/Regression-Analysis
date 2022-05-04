env_data = read.table("Environmental_data.txt", sep = " ", header = TRUE)
radiation = env_data$radiation
temp = env_data$temperature
wind = env_data$wind
ozone = env_data$ozone
#(A)
pairs(env_data)

#(B)(ii)
library(matlib)
intercept = rep(1, length(wind))
X = matrix(c(intercept, radiation, temp, wind), ncol = 4, byrow = F)
Xt = t(X)
G = Xt %*% X
solve(Xt %*% X)

Y = env_data$ozone
Xt %*% Y

t(Y) %*% Y

#(B)(iii),(vi)
beta.hat = solve(Xt %*% X) %*% Xt %*% Y
fit.model = lm(ozone~radiation + temp + wind); fit.model
fit.model.sum = summary(fit.model)

#B(iv)
anova(fit.model)

#B(v)
y.bar = rep(mean(ozone), times = 30)
Ftest = ((t(Y -y.bar) %*% (Y- y.bar) -SSEr)/3)/(SSEr/(30-4))
pf(Ftest, 3, 30-4)

#B(vi)
sstc = t(Y -y.bar) %*% (Y- y.bar) - 30*mean(ozone)^2
R.sq = (sstc- SSEr)



#B(vii)
f0 = c(1, 100-50, 70-80, 10-10)
f0 %*% beta.hat + pt(0.025/2, 30-4)*sigma(fit.model)*sqrt(1 + t(f0)%*%solve(G)%*%f0)
f0 %*% beta.hat - pt(0.025/2, 30-4)*sigma(fit.model)*sqrt(1 + t(f0)%*%solve(G)%*%f0)

#B(viii)
par(mfrow=c(2,2), mai = c(0.8, 0.8, 0.5, 0.3))
plot(fitted.values(fit.model), residuals(fit.model), xlab = "standardised residuals", ylab = "fitted values")
abline(h = 0, col = "red")

plot(radiation, residuals(fit.model), xlab = "standardised residuals", ylab = "radiation")
abline(h = 0, col = "red")

plot(temp, residuals(fit.model), xlab = "standardised residuals", ylab = "temperature in degrees fahrenheit")
abline(h = 0, col = "red")

plot(wind, residuals(fit.model), xlab = "standardised residuals", ylab = "wind speed in miles per hour (mph)")
abline(h = 0, col = "red")

#C(ii)
tstat = qt(0.95/2, 30-4)
beta.hat[3] + tstat* sigma(fit.model)* sqrt(G[2,2])
beta.hat[3] - tstat* sigma(fit.model)* sqrt(G[2,2])


#D(i)
Xf = matrix(c(intercept, radiation, temp, wind, radiation*temp, radiation*wind, temp*wind), ncol = 7, byrow = F); Xf
Xft = t(Xf); Xft
Gf = Xft %*% Xf
solve(Xft %*% Xf)

Xft %*% Y

betaf.hat = solve(Xft %*% Xf) %*% Xft %*% Y

ffit.model = lm(ozone~radiation + temp + wind + radiation*temp + radiation*wind + temp*wind)
ffit.model.sum = summary(ffit.model)

#D(iii)
ub = rep(0, times = 7)
for (i in 1:7) {
  ub[i] = betaf.hat[i,1] + pt(0.95/2, 30-7) * sigma(ffit.model)*sqrt(Gf[i,i])
}

lb = rep(0, times = 7)
for (i in 1:7) {
  lb[i] = betaf.hat[i,1] - pt(0.95/2, 30-7) * sigma(ffit.model)*sqrt(Gf[i,i])
}
lb
ci = cbind(lb, ub); ci
confint(ffit.model)

#D(iv)
SSEr = t(residuals(fit.model)) %*% residuals(fit.model)
SSEf = t(residuals(ffit.model)) %*% residuals(ffit.model)
Fstat = ((SSEr[1,1]- SSEf[1,1])/3)/sigma(ffit.model)^2





