source("/Users/lilygebhart/Downloads/Summer 23 R Code/yarnold_implementation.R")

# Calculate approximate chi square distribution from Kruskal-Wallis Test Statistic, assumes
# same number of entries in each group, not used elsewhere.
calc_grpsum_chisq <- function(nv, c2) {
  len <- length(nv)
  avg_rank <- (length(nv) + 1)/2
  t <- rep(0, len)
  for (i in 1:len){
    t[i] <- nv[i]*avg_rank
  }
  S <- t[2:len]
  tr_S <- t(S)
  inv_c2 <- solve(c2)
  return(list(grpsum = S, chisq = tr_S %*% inv_c2 %*% S))
}

# Yarnold simulation for different c values given cumulants. 
calc_yarnold_sim <- function(m2, c1, c2, c3, c4, nv, c_val, n = 1) {
  len <- length(nv)
  k <- len - 1
  center <- rep(0, len)
  for (jj in 1:len){
    center[jj] <- (sum(nv)+1)*nv[jj]/2
  }
  yarnold_stat <- rep(NA, length(c_val))
  chisq_stats <- pchisq(c_val,df=length(nv) - 1)
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


# Produces a plot of Yarnold approximation for KW test vs. KW test
# statistic vs. the Chi Square distribution. 
# Uses loopme to calculate moments
yarnold_KW_chisq_comp <- function(nv, c_val) {
  KW_stats <- rep(NA, 100000)
  for (jj in seq(length(KW_stats))) {
    KW_stats[jj] <- kruskal.test(rnorm(length(nv)*nv[1]), c(rep(1:length(nv), nv[1])))$statistic
  }
  moments <- loopme(nv)
  cumu <- calc_cumulants(moments$k1, moments$k2, moments$k3, moments$k4)
  yarnold_stats <- calc_yarnold_sim(moments$k2, cumu$c1, cumu$c2, cumu$c3, cumu$c4, nv, c_val)
  chisq_stats <- pchisq(c_val,df=length(nv) - 1)
  plot(c_val, chisq_stats, type="l", col="blue",
       xlim=c(0, max(KW_stats)), ylim = c(0, 1))
  lines(c_val, yarnold_stats, col="brown")
  lines(sort(KW_stats), seq(length(KW_stats))/length(KW_stats))
  legend(x=0, y=1, legend=c("Kruskal-Wallis", "Chi-Square", "Yarnold"),
         fill= c("black", "blue", "brown"), cex=0.6)
}


# Produces a plot of Yarnold approximation for KW test vs. KW test
# statistic vs. the Chi Square distribution. 
# Uses central_moments_calc() to calculate moments, which uses truetrue() function.
yarnold_KW_chisq_comp2 <- function(nv, c_val) {
  KW_stats <- rep(NA, 100000)
  for (jj in seq(length(KW_stats))) {
    KW_stats[jj] <- kruskal.test(rnorm(sum(nv)), c(rep(1:length(nv), nv[1])))$statistic
    }
  moments <- central_moments_calc(nv)
  cumu <- calc_cumulants(moments$m1, moments$m2, moments$m3, moments$m4)
  yarnold_stats <- calc_yarnold_sim(moments$m2, cumu$c1, cumu$c2, cumu$c3, cumu$c4, nv, c_val)
  chisq_stats <- pchisq(c_val,df=length(nv) - 1)
  plot(c_val, chisq_stats, type="l", col="blue",
       xlim=c(0, max(KW_stats)), ylim = c(0, 1))
  lines(c_val, yarnold_stats, col="brown")
  lines(sort(KW_stats), seq(length(KW_stats))/length(KW_stats))
  legend(x=0, y=1, legend=c("Kruskal-Wallis", "Chi-Square", "Yarnold"),
         fill= c("black", "blue", "brown"), cex=0.6)
}

#Testing...
yarnold_KW_chisq_comp2(c(5, 5, 5, 5), (seq(from=1, to=100))/10)






