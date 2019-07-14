# Internal functions used to passing and modifying optional arguments
# within functions
# 1. replace must be a list of the form
# list(A = newval)
# 2. remove is a list of variables of the form c("x","y")
# Ariane Schad, Feb 2019

.optargsReplace <- function(..., replace){
    opt.args <- list(...)
    opt.args[[names(replace)]] <- replace[[1]]
    return(opt.args)
}

.optargsRemove <- function(..., remove){
    opt.args <- list(...)
    opt.args[remove] <- NULL
    return(opt.args)
}
