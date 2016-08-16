## ----fig.align='center'--------------------------------------------------
library("DiffusionRimp")
# Set the parameters of the problem:
Xs   <- 1
s    <- 0
t    <- 5
lims <- c(-3, 3)
delt <- 0.0025
nodes<- 101

# Define drift and diffusion terms:
mu  <- function(X,t){-X^3}
sig <- function(X,t){1}
# Approximate the transitional density:
res <- MOL.density(Xs, s, t, lims, nodes, delt)

# Plot the density:
persp(res$Xt, res$time, pmin(res$density, 1), col = 'white', xlab = 'State (X_t)',ylab='Time (t)',
    zlab='Density f(X_t|X_s)', border = NA, shade = 0.5, theta = 145)

## ----fig.align='center'--------------------------------------------------
# Define drift and diffusion terms:
mu  <- function(X,t){X-X^3}
sig <- function(X,t){1}
# Approximate the transitional density:
res <- MOL.density(Xs, s, t, lims, nodes, delt)

# Plot the density:
persp(res$Xt, res$time, pmin(res$density, 1), col = 'white', xlab = 'State (X_t)',ylab='Time (t)',
    zlab='Density f(X_t|X_s)', border = NA, shade = 0.5, theta = 145)

## ----fig.align='center'--------------------------------------------------
# Define drift and diffusion terms:
mu  <- function(X,t){X*((1+cos(2*pi*t))-X^2)}
sig <- function(X,t){(1+0.25*sin(3*pi*t))}
# Approximate the transitional density:
res <- MOL.density(1, 0, 5, c(-3, 3), 101, 0.0025)

# Plot the density:
for(i in 0:1)
{
  persp(res$Xt, res$time, pmin(res$density, 1), col = 'white', xlab = 'State (X_t)',ylab='Time (t)',
    zlab='Density f(X_t|X_s)', border = NA, shade = 0.5, theta = 145+i*60)
}

## ----fig.align='center'--------------------------------------------------
 # Set the parameters of the problem:
 Xs     <- 1            # Starting X-coordinate
 Ys     <- 1            # Starting Y-coordinate
 s      <- 0            # Starting time
 t      <- 10           # Final horizon time
 Xlim   <- c(-2.2,2.2)  # Lattice endpoints in X dim
 Ylim   <- c(-2.2,2.2)  # Lattice endpoints in Y dim
 Nodes  <- 51           # How many nodes (incl. ends)
 delt   <- 1/100        # Time stepsize
 
 # Define drift and diffusion terms:
 mu1   <- function(X,Y,t){X*(1 - X^2) + sin(0.5*pi*t)*Y}
 mu2   <- function(X,Y,t){Y*(1 - Y^2) - sin(0.5*pi*t)*X}
 sig11 <- function(X,Y,t){0.5}
 sig22 <- function(X,Y,t){0.5}
 
 # Run the Method of Lines:
 res <- BiMOL.density(Xs, Ys, s, t, Xlim, Ylim, Nodes, delt)
 
 library("colorspace")
 colpal=function(n){rev(c(sequential_hcl(n-1,power=0.8,l=c(40,100)),'white'))}
  
 for(i in c(251,501,751,1001))
 {
   filled.contour(res$Xt, res$Yt, res$density[,,i], color.palette = colpal,
                  main = paste0('Transition Density \n (t = ', res$time[i],')'),
                  xlab='Xt',ylab='Yt')
 }
 

## ----fig.align='center'--------------------------------------------------
# Define drift and diffusion terms:
mu1   <- function(X,Y,t){-Y*sin(X*pi)} 
mu2   <- function(X,Y,t){-X*sin(Y*pi)}
sig11 <- function(X,Y,t){0.5}
sig22 <- function(X,Y,t){0.5}

# Parameters of the problem:
X0  <- 0        # Initial X-coordinate
Y0  <- 0        # Initial Y-coordinate
s   <- 0        # Starting time
t   <- 2.5      # Final horizon time
Xlim <- c(-3,3) # Lattice endpoints in X dim
Ylim <- c(-3,3) # Lattice endpoints in Y dim
Nodes <- 121    # How many nodes per dimension (incl. ends)
delt  <- 1/200  # Step size

# Run the Method of Lines:
res <- BiMOL.density(X0, Y0, s, t, Xlim, Ylim, Nodes, delt)

time.sequence <- c(0.5,1,1.5,2,2.5)
k = 0
for(i in time.sequence)
{
  k = k+1
  persp(res$Xt,res$Yt,res$density[,,i/delt],col='white',xlab='State (X_t)',ylab='State (Y_t)',zlab='Density',border=NA,shade=0.5,theta=145+180*k/5)
}

## ----eval=FALSE----------------------------------------------------------
#  browseVignettes('DiffusionRimp')

