# wrappers for commonly used visualization tools

sim_hist <- function(counts, observed, xlab = "Simulation Outcomes", main = "Comparing Observation to Simulation", col = "cyan"){
  
  # makes a histogram of simulation outcomes and plots observed counts as dotted black line
  
  par(mar=c(5,4,5,1) + 0.5)   # extra large bottom margin
  breaks = seq(min(counts) - 0.5,max(max(counts) + 0.5, observed),1)
  h = hist(counts, xlab=xlab, main=main, breaks = breaks, col=col, xaxt="n")
  axis(side=1, at = breaks+0.5, labels=breaks + 0.5)
  abline(v=observed, col="black", lty=3, lw=5)
}