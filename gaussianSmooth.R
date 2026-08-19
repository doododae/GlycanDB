gaussianSmooth <- function(x, window, sd = window / 4) {
  # window should be odd
  half <- floor(window / 2)
  kernel <- dnorm(-half:half, mean = 0, sd = sd)
  kernel <- kernel / sum(kernel)  # normalize so it sums to 1
  
  stats::filter(x, kernel, sides = 2)
}