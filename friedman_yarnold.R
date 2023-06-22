source("/Users/lilygebhart/Downloads/Summer 23 R Code/yarnold_implementation.R")

# Calculates moments for the Friedman test statistic distribution under the null hypothesis. 
friedman_mom_calc <- function(num_groups, num_blocks, num_per){
  moms <- central_moments_calc(rep(num_per, num_groups))
  return(list(m1=moms$m1 * num_blocks, m2=moms$m2 * num_blocks,
              m3=moms$m3 * num_blocks, m4=moms$m4 * num_blocks))
}

# Yarnold simulation for different c values given cumulants. 
calc_yarnold_sim_F <- function(num_groups, num_blocks, num_per, m2, c1, c2, c3, c4, c_val, n = 1) {
  len <- num_groups
  k <- len - 1
  center <- rep(0, num_groups)
  for (jj in 1:len) {
    center[ii] <- num_blocks*(num_blocks*num_groups*num_per+1)*num_per/2
  }
  yarnold_stat <- rep(NA, length(c_val))
  chisq_stats <- pchisq(c_val,df=num_treatments - 1)
  for (ii in 1:length(c_val)) { 
    c <- c_val[ii]
    num_per_ellipse <- countl(center=center, nn=k, kk=k, rsq=c, a=chol(c2), v=rep(0, k))
    vol_ellipse <- vol_ellipsoid(apply(m2,2,c), n, c, k)
    term2_ellipse <- num_per_ellipse - vol_ellipse
    term2 <- term2_ellipse * (exp(-c/2)/((2*pi*n)^(k/2) * det(as.matrix(c2,2,c))^(1/2)))
    term3_d1 <- (calc_delta1(c2, c4)/ n) *
      (choose(2, 0)* dchisq(c, k) - choose(2, 1)* dchisq(c, k+2) + choose(2, 2)* dchisq(c, k+4))
    term3_d2 <- (calc_delta2(c2, c3)/2) * 
      (-choose(3, 0)*dchisq(c, k) + choose(3, 1)*dchisq(c, k+2) - choose(3, 2)*dchisq(c, k+4) +
         choose(3, 3)*dchisq(c, k+6))
    yarnold_stat[ii] <- chisq_stats[ii] + term2 - (term3_d1 + term3_d2)
  }
  return(yarnold_stat)
}

# Produces a plot of Yarnold approximation for Friedman test vs. Friedman test
# statistic vs. the Chi Square distribution. 
yarnold_Friedman_chisq_comp <- function(num_groups, num_blocks, num_per, c_val) {
  Friedman_stats <- rep(NA, 100000)
  for (jj in seq(length(Friedman_stats))) {
    Friedman_stats[jj] <- friedman.test(matrix(rnorm(num_blocks*num_groups*num_per),
                                               nrow=num_blocks, byrow=TRUE))$statistic
    
  }
  moments <- friedman_mom_calc(num_groups, num_blocks, num_per)
  cumu <- calc_cumulants(moments$m1, moments$m2, moments$m3, moments$m4)
  yarnold_stats <- calc_yarnold_sim_F(num_groups, num_blocks, num_per, moments$m2,
                                      cumu$c1, cumu$c2, cumu$c3, cumu$c4, c_val)
  chisq_stats <- pchisq(c_val,df=num_treatments - 1)
  plot(c_val, yarnold_stats, col="brown", type="l")
       #xlim=c(0, min(c(max(Friedman_stats), max(chisq_stats)))), ylim = c(0, 1))
  lines(c_val, chisq_stats, col="blue")
  lines(sort(Friedman_stats), seq(length(Friedman_stats))/length(Friedman_stats))
  legend(x=0, y=1, legend=c("Friedman", "Chi-Square", "Yarnold"), fill= c("black", "blue", "brown"), cex=0.6)
}


# Testing...
yarnold_Friedman_chisq_comp(4, 4, 1, seq(from=1, to=1000, by=5)/ 100)

