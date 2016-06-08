## ------------------------------------------------------------------------
mu  <- function(X,t){-alpha*X*sin(X*pi)}
sig <- function(X,t){sigma}

## ------------------------------------------------------------------------
library(DiffusionRimp)

alpha <- 1
sigma <- 0.5
res   <- MOL.density(0.1, 0, 5, c(-4, 4), 101, 0.01)

## ----fig.align = 'center'------------------------------------------------
persp(x = res$Xt, y = res$time, z = res$density, col = 'white', xlab = 'State (X_t)',ylab='Time (t)', 
      zlab='Density f(X_t|X_s)', border = NA, shade = 0.5, theta = 145)

## ----eval=FALSE----------------------------------------------------------
#  browseVignettes('DiffusionRimp')

