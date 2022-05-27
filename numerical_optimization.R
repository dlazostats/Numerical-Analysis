# Numerical Optimization
#----------------------
# the golden section search
golden <- function (f, a, b, tol = 0.0000001){
  ratio <- 2 / (sqrt(5) + 1)
  x1 <- b - ratio*(b - a)
  x2 <- a + ratio*(b - a)
  f1 <- f(x1)
  f2 <- f(x2)
  while(abs(b - a) > tol) {
    if (f2 > f1) {
      b <- x2
      x2 <- x1
      f2 <- f1
      x1 <- b - ratio*(b - a)
      f1 <- f(x1)
    } else {
      a <- x1
      x1 <- x2
      f1 <- f2
      x2 <- a + ratio*(b - a)
      f2 <- f(x2)
    }
  }
  return((a + b)/2)
}
### f1
f<-function(x){
  abs(x-3.5)+(x-2)^2
}
curve(f,from=1,to=5)
golden(f,1,5)

g<-function(x){
  -x^3+4*x+1
}
curve(g,from=-6,to=6)
golden(g,-6,6)

w<-function(x){
  -x^2+4*x+1
}
curve(w,from=-6,to=8)
golden(g,-6,8)

# Newton-Raphson
f <- function(x) exp(-x) + x^4
curve(f, from = -1, to = 4)

f <- function(x) exp(-x) + x^4
fprime <- function(x) -exp(-x) + 4 * x^3
fprimeprime <- function(x) exp(-x) + 12 * x^2
x <- c(0.6, rep(NA, 6))
fval <- rep(NA, 7)
fprimeval <- rep(NA, 7)
fprimeprimeval <- rep(NA, 7)
for (i in 1:6) {
  fval[i] <- f(x[i])
  fprimeval[i] <- fprime(x[i])
  fprimeprimeval[i] <- fprimeprime(x[i])
  x[i + 1] <- x[i] - fprimeval[i] / fprimeprimeval[i]
}
data.frame(x, fval, fprimeval, fprimeprimeval)

x<-seq(-1,4,0.01)
y<-f(x)
min(y)
df<-data.frame(x,y)
plot(y,type="l")

## Gradient descent for local minima
graddsc <- function (fp , x, h = 1e-3, tol = 1e-4, m = 1e3){
  iter <- 0
  oldx <- x
  x = x - h*fp(x)
  while (abs(x - oldx) > tol) {
    iter <- iter + 1
    if( iter > m)
      stop ("No solution found ")
    oldx <- x
    x = x - h*fp(x)
  }
  return (x)
}
fp <- function(x) { x^3 + 3*x^2 - 1 }
curve(fp,from = -10,to=10)
graddsc(fp, x = 1, h = 1/100)
