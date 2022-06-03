euler <- function (f, x0 , y0 , h, n) {
  x <- x0
  y <- y0
  for(i in 1:n) {
    y0 <- y0 + h*f(x0 , y0)
    x0 <- x0 + h
    x <- c(x, x0)
    y <- c(y, y0)
  }
  return ( data.frame (x = x, y = y))
}