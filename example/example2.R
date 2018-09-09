# Examples used in Bioinformatics Applications Note, Supplemental Material
# MBQN: R package for mean balanced quantile-normalization, Bioinformatics, 2018

# Generate a simulated data matrix and apply MBQN
# including 2 conditions, 10 replicates, and dep and ndep proteins

# A. Schad, April 2018

# Model 1: all features are RI feature

N <- 20
M <- 10
scale <- 0.45

#dat0 <- matrix(NA, nrow = N, ncol = M)

x <- 1:N
dat0 <- matrix(rep(x,each=M), ncol=M, byrow=TRUE)
m <- (matrix(runif(N*M, min = 0, max = 1),N,M)-0.5)*scale

expMatrix <- dat0 + m

# add batch effects, i.e. transform to mp <- m*s+c
c <- runif(M, min = -1, max = 1)
s <- runif(M, min = 0.5, max = 5)
mp <- sweep(m, 2, s, "*")
mp <- sweep(mp, 2, c, "+")

boxplot(m)
boxplot(mp)

filename <- "example_RI_features"

library("MBQN")
devtools::document()
res <- mbqn.check_saturation(expMatrix, FUN = median, flag_show_fig = TRUE, low_thr = 0.2, filename = filename)



# Model 2: RI features, NRI features, and well mixing features

filename <- ""
#expMatrix <- matrix(xxxx)
expMatrix <- log2(expMatrix)

library("MBQN")
devtools::document()
res <- mbqn.check_saturation(expMatrix, FUN = median, flag_show_fig = TRUE, low_thr = 0.2, filename = filename, feature_index = 124)

# get index of nri/ri feature with largest rank invariance frequency
nri_max <-as.numeric(names(which.max(res$nri)))

