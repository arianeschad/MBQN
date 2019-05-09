#' Demonstrate mean/median-balanced quantile normalization
#'
#' @description Visualize the intensity distribution of unnormalized,
#' quantile and mean/median-balanced quantile normalized data.
#' @param x a matrix where rows represent features, e.g. of protein
#' abundances/intensities, and columns represent samples. Default
#' \code{NULL} - a simple matrix is generated.
#' @importFrom stats median
#' @return Various graphics illustrating the effect of normalization on
#' rank mixing and rank invariant intensity features.
#' @family example
#' @examples
#' # Simple example:
#' mbqnDemo()
#' # Example with a single RI feature:
#' set.seed(1234)
#' x <- mbqnSimuData(model = "omics.dep")
#' mbqnDemo(x)
#' @details Normalize matrix and return boxplots of quantile normalized and mean balanced normalized
#' data.
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean/median-balanced quantile normalization. In prep. 2019
#' @author Ariane Schad
# Aug. 2017
#' @export mbqnDemo

mbqnDemo <- function(x = NULL){

  # if no matrix is given, create a simple dummy matrix
  if(is.null(x)){
    x <- matrix(c(5,2,3,NA,2,4,1,4,2,3,1,4,6,NA,1,3,NA,1,4,3,NA,1,2,3),ncol=4)

    print(x)
    #      [,1] [,2] [,3] [,4]
    # [1,]    5    1    6    4
    # [2,]    2    4   NA    3
    # [3,]    3    2    1   NA
    # [4,]   NA    3    3    1
    # [5,]    2    1   NA    2
    # [5,]    4    4    1    3

  }

  # perform qn, median balanced qn, and qn with median balanced nri feature
  qn_x <- mbqn(x,FUN=NULL, verbose = FALSE)
  mbqn_x <- mbqn(x,FUN = "median", verbose = FALSE)
  qn_nri_x <- mbqnNRI(x,FUN = "median", low_thr = 0.5, verbose = FALSE)

  # check saturation i.e. for rank invariance
  res <- mbqnGetNRIfeatures(x,verbose = FALSE)

  mtx <- matrix(c(2, 3, 3, 6,
                  4, 1, 1, 5), byrow=TRUE, nrow=2)
  nf <- layout(mtx, heights=c(1,1), widths=c(6,5,1,0.5))

  ylim <- range(x, na.rm =TRUE)
  # create a boxplot for qn-data with nri features median balanced
  mbqnBoxplot(qn_nri_x, irow = res$ip,
              add.leg = TRUE, main ="QN data with MBQN NRI feature",
              cex.lab = 1.2,
              cex = 1.2,
              ylim = ylim)
  # create a boxplot for x
  mbqnBoxplot(x,
              add.leg = FALSE, main = "data",
              cex.lab = 1.2,
              cex = 1.2,
              ylim = ylim)
  # create a boxplot for qn-data
  mbqnBoxplot(qn_x, irow = res$ip,
              add.leg = FALSE, main = "QN data",
              cex.lab = 1.2,
              cex = 1.2,
              ylim = ylim)
  # create a boxplot for mbqn-data
  mbqnBoxplot(mbqn_x,
              add.leg = FALSE, main = "MBQN data",
              cex.lab = 1.2,
              cex = 1.2,
              ylim = ylim)

}

