# Script file to create example file for NA pattern
# A. Schad

MBQN:::mbqnLoadExample(which.example = 1)
mtx <- !is.na(mtx$mtx[,c(1:9,19:27)])
colnames(mtx) <- rownames(mtx) <- NULL
six <- sort(apply(mtx,1,sum),decreasing = TRUE, index.return = TRUE)
mtx <- ifelse(mtx[six$ix,], 1,0)
example_NApattern <- mtx
#usethis::use_data(example_NApattern, internal = TRUE, overwrite = TRUE)
usethis::use_data(example_NApattern, overwrite = TRUE)
