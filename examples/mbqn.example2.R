# Examples used in Bioinformatics Applications Note, Supplemental Material
# MBQN: R package for mean balanced quantile-normalization, Bioinformatics, 2018

# Simulate a data matrix, add batch effects and apply MBQN
# Matrix includes 2 conditions, 10 replicates, and 10 % dep ( = diff. expr. proteins) and
# ndep ( = not diff. expr. proteins)

# A. Schad, April 2018

###############################################################################
# Model: Matrix with RI features, NRI features, and well mixing features

N <- 200
M <- 10
x <- 1:N
#scale <- 0.45
#dat0 <- matrix(rep(x,each=M), ncol=M, byrow=TRUE)
#set.seed(1234)
#m <- (matrix(runif(N*M, min = 0, max = 1),N,M)-0.5)*scale

#mtx <- dat0 + m

# simulate batch effects, i.e. transform to mp <- m*s+c
#set.seed(2345)
#c <- runif(M, min = -1, max = 1)
#set.seed(3456)
#s <- runif(M, min = 0.5, max = 5)
#mp <- sweep(m, 2, s, "*")
#mp <- sweep(mp, 2, c, "+")


filename <- "example_mixed_features"
dat0 <- matrix(rep(x,each=M), ncol=M, byrow=TRUE)/200

scale <- .45

set.seed(1234)
m <- (matrix(runif(N*M, min = 0, max = 1),N,M)-0.5)*scale

scale <- 1
m[1:5,] <- matrix(rep(x[1:5],each=M), ncol=M, byrow=TRUE)/sqrt(2)
m[1:5,] <- m[1:5,] + (matrix(runif(5*M, min = 0, max = 1),5,M)-0.5)*scale
mtx <- dat0 + m

res <- MBQN::mbqn.check_saturation(mtx, FUN = median,
                                   show_fig = TRUE, low_thr = 0.5,
                                   filename = filename, feature_index = 124)

res

# get index of nri/ri feature with largest rank invariance frequency
nri_max <-as.numeric(names(which.max(res$nri)))
print(paste("Index of maximum rank invariant feature:", nri_max))
