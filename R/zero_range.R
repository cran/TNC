# function by Hadley Wickham
# https://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-vector

zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}
