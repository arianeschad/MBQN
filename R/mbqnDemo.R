#' Demonstration of mean/median-balanced quantile normalization
#'
#' @description Visualize the intensity distribution of unnormalized, quantile and mean/median-balanced quantile normalized data
#' @param dat a matrix where rows represent features, e.g. of protein
#' abundances/intensities, and columns represent samples.
#' Default \code{NULL} - a simple matrix is generated.
#' @importFrom stats median
#' @return Various graphics illustrating the effect of normalization on
#' rank mixing and rank invariant intensity features
#' @family data
#' @examples
#' mbqnDemo()
#' dat <- mbqnSimuData()
#' mbqnDemo(dat)
#' @details Normalize matrix and return boxplots of quantile normalized and mean balanced normalized
#' data.
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean balanced quantile normalization. In prep. 2019
#' @author Ariane Schad
# Aug. 2017
#' @export mbqnDemo

mbqnDemo <- function(dat = NULL){

  # if no matrix is given, create a simple dummy matrix
  if(is.null(dat)){
    dat <- matrix(c(5,2,3,NA,2,4,1,4,2,3,1,4,6,NA,1,3,NA,1,4,3,NA,1,2,3),ncol=4)

    print(dat)
    #      [,1] [,2] [,3] [,4]
    # [1,]    5    1    6    4
    # [2,]    2    4   NA    3
    # [3,]    3    2    1   NA
    # [4,]   NA    3    3    1
    # [5,]    2    1   NA    2
    # [5,]    4    4    1    3

  }

  # perform qn, median balanced qn, and qn with median balanced nri feature
  qn_dat <- mbqn(dat,FUN=NULL, verbose = FALSE)
  mbqn_dat <- mbqn(dat,FUN = median, verbose = FALSE)
  qn_nri_dat <- mbqnNRI(dat,FUN = "median", low_thr = 0.5, verbose = FALSE)

  # check saturation i.e. for rank invariance
  res <- mbqnCheckSaturation(dat, save_fig = FALSE, verbose = FALSE, show_fig = FALSE)

  dev.off()
  plot.new()
  frame()
  mtx <- matrix(c(2, 3, 3, 6,
                  4, 1, 1, 5), byrow=TRUE, nrow=2)
  nf <- layout(mtx, heights=c(1,1), widths=c(6,5,1,0.5))

  ylim <- range(dat, na.rm =T)
  # create a boxplot for qn-data with nri features median balanced
  mbqnBoxplot(qn_nri_dat, irow = res$ip, filename = NULL,
              add.leg = T, main ="QN data with MBQN NRI feature", ylim = ylim)
  # create a boxplot for dat
  mbqnBoxplot(dat, filename = NULL,
              add.leg = F, main = "data", ylim = ylim)
  # create a boxplot for qn-data
  mbqnBoxplot(qn_dat, filename = NULL, irow = res$ip,
              add.leg = F, main = "QN data", ylim = ylim)
  # create a boxplot for mbqn-data
  mbqnBoxplot(mbqn_dat, filename = NULL,
              add.leg = F, main = "MBQN data", ylim = ylim)

}

