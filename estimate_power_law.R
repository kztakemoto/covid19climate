estimate.power.law <- function(count){
    #https://mathworld.wolfram.com/LeastSquaresFittingPowerLaw.html
    y <- count[count > 0]
    x <- as.integer(as.Date(names(y)) - as.Date(names(y)[[1]]) + 1)
    b <- (length(x) * sum(log(x) * log(y)) - sum(log(x)) * sum(log(y))) / (length(x) * sum(log(x)**2) - sum(log(x))**2)
    return(b)
}