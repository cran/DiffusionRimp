## ---- fig.align='center'-------------------------------------------------
library("DiffusionRimp")

mu  <- function(X){(0.5*X*(1-X^2))}
sig <- function(X) {0.5}
res <- MOL.passage(Xs = 0.5,t = 10, limits=c(-0.5, 2), N = 51, delt = 1/250)
  
persp(res$Xt, res$time, res$surface, col = 'white', xlab = 'X_s',
      ylab = 'Time (t)', zlab = 'Survival Prob.', border = NA, shade = 0.5,
      theta = 145, phi = 30, ticktype = 'detailed')
  
plot(res$density~res$time, col = '#222299', type='l', main = 'First passage time density')

#-------------------------------------------------------------------------------

## ---- fig.align='center'-------------------------------------------------
#===============================================================================
# Bi-cubic diffusion with concentration in the even quadrants.
#===============================================================================
 
 # Define drift and diffusion terms.
 mu1   <- function(X,Y){0.5*(1 - X^2)*X - Y}
 mu2   <- function(X,Y){0.5*(1 - Y^2)*Y - X}
 sig11 <- function(X,Y){1}
 sig22 <- function(X,Y){1}

 # Parameters of the problem.
 Xs    <- 0.5      # Starting X-coordinate
 Ys    <- 0.5      # Starting Y-coordinate
 t     <- 10       # Final horizon time
 Xbar  <- c(-2, 2) # Limits in X-dim
 Ybar  <- c(-2, 2) # Limits in Y-dim
 Nodes <- 51       # How many nodes x nodes (incl. ends)
 delt  <- 1/250    # Time step size

 res <- BiMOL.passage(Xs, Ys, t, c(Xbar, Ybar), Nodes, delt)
 
 time.sequence <- c(0.5, 2, 4, 6, 8, 10)
 k = 0
 for(i in time.sequence)
 {
  k = k+1
  persp(res$Xt, res$Yt, res$surface[,,i/delt], col='white',
        xlab = 'State (X_0)', ylab = 'State (Y_0)', 
        zlab = 'Survival Probability', border = NA, 
        shade = 0.5, theta = 145 + 180*k/6)
 }
 
  plot(res$density~res$time, col = '#222299', type = 'l',
       main = 'First passage time density')
#-------------------------------------------------------------------------------

## ----fig.align = 'center'------------------------------------------------
# Define drift and diffusion terms:
mu1   <- function(X,Y){-1*X}
mu2   <- function(X,Y){-1*Y}
sig11 <- function(X,Y){0.5}
sig22 <- function(X,Y){0.5}

Nodes <- 101
delt  <- 1/500

# Define a region in the x-y plane:
region <- function(x,y){sqrt((x-0.5)^2+(y-0.5)^2)<1}
res    <- BiMOL.passage(0.5, 0.5, 10, c(-1,2,-1,2), Nodes, delt, Phi = region)

# Draw a nice plot of the evolution of the survival prob.:
library("colorspace")
colpal=function(n){rev(sequential_hcl(n,power=0.8,l=c(40,100)))}
time.sequence <- c(0.1, 0.25, 0.5, 1)
k = 0
for(i in time.sequence)
{
  k = k+1
  filled.contour(res$Xt, res$Yt, res$surface[,,i/delt], color.palette = colpal,
                main = paste0('Survival probability \n (t = ', i,')'),
                xlab = expression(X[0]), ylab = expression(Y[0]),
                plot.axis = {
                th =seq(0,1,1/100)
                lines(0.5+sin(2*pi*th),0.5+cos(2*pi*th))
                abline(h = 0, v = 0, lty = 'dashed')
                })

}


## ---- fig.align= 'center'------------------------------------------------
sim  = function(N = 50000, plt = FALSE)
{
   set.seed(200707881)
   delt = 1/2000
   times = rep(0, N)
   X = rep(0.5, N)
   Y = rep(0.5, N)
   k = 1
   d = 0
   while(N >1)
   {
     X = X + mu1(X,Y)*delt + sig11(X,Y)*rnorm(N,sd = sqrt(delt))
     Y = Y + mu2(X,Y)*delt + sig22(X,Y)*rnorm(N,sd = sqrt(delt))
     d = d + delt
     wh = (sqrt((X-0.5)^2+(Y-0.5)^2)>=1)
     if(plt)
     {
       plot(X,Y,pch = 20,xlim = c(-1,2), ylim =c(-1,2))
                     th =seq(0,1,1/100)
                lines(0.5+sin(2*pi*th),0.5+cos(2*pi*th))
     }
     if(any(wh))
     {
        X = X[which(!wh)]
        Y = Y[which(!wh)]
        m = N - length(X)
        N = length(X)
        times[k:(k+m-1)] =  d
        k = k+m
     }
   }
   return(times)
}

sm = sim()

hist(sm, freq = FALSE, border = 'white', col= 'gray75',
     ylim = range(res$density), xlim = c(0, 6))
lines(res$density~res$time, col = '#222299', lwd = 2)

## ----eval=FALSE----------------------------------------------------------
#  browseVignettes('DiffusionRimp')

