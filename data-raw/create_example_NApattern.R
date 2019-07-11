# Script file to create example file for NA pattern
# A. Schad

pxd_id <- "PXD001584"

# Load file
out <- mbqnLoadFile(pxd_id, file.pattern = "proteinGroups.txt")

# filter for potential contaminants and identified only by site features
out <- out[!rowData(out)[["ixs"]],] 

# extract data and feature annotation
mtx <- assays(out)[["data"]]

mtx <- !is.na(mtx[,c(1:9,19:27)])
colnames(mtx) <- rownames(mtx) <- NULL
six <- sort(apply(mtx,1,sum),decreasing = TRUE, index.return = TRUE)
mtx <- ifelse(mtx[six$ix,], 1,0)
example_NApattern <- mtx
#usethis::use_data(example_NApattern, internal = TRUE, overwrite = TRUE)
usethis::use_data(example_NApattern, overwrite = TRUE)
