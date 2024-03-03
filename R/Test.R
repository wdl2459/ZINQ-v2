library(quantreg)
n = 300
p <- function(x0, gam0=0.75, gam1=-0.15){
   lc = gam0 + gam1*x0
   exp(lc) / (1 + exp(lc))
}
x = c(rep(0, n), rep(1, n))
w = 0.5 + 1.5*x + (1+0.15*x)*rchisq(2*n,df=1)
b = rbinom(2*n, 1, p(x))
y = w*b
dat = data.frame(y, x)


ZINQ_tests(formula.logistic=y~x, formula.quantile=y~x, C="x", data=dat)


formula.logistic=y~x
formula.quantile=y~x
C="x"
data=dat





library(quantreg)
n = 300
p <- function(x0, x1, gam0=0.75, gam1=-0.15, gam2=0.10){
  lc = gam0 + gam1*x0 + gam2*x1
  exp(lc) / (1 + exp(lc))
}
x0 = c(rep(0, n), rep(1, n))
x1 = rnorm(2*n)
w = 0.5 + 1.5*x0 + 0.1*x1 + (1+0.15*x0)*rchisq(2*n,df=1)
b = rbinom(2*n, 1, p(x0, x1))
y = w*b
dat = data.frame(y, x0, x1)


ZINQ_tests(formula.logistic=y~x, formula.quantile=y~x, C="x", data=dat)


formula.logistic=y~x0 + x1
formula.quantile=y~x0 + x1
C=c("x0", "x1")
data=dat

