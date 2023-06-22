source("/Users/lilygebhart/Downloads/Summer 23 R Code/loop3.R")
source("/Users/lilygebhart/Downloads/Summer 23 R Code/countl.R")

# Calculates 1st - 4th cumulants from 1st - 4th moments
calc_cumulants<-function(m1, m2, m3, m4) {
  len <- dim(m2)[1]
  m1 <- m1[2:len]
  m2 <- m2[2:len, 2:len]
  m3 <- m3[2:len, 2:len, 2:len]
  m4 <- m4[2:len, 2:len, 2:len, 2:len]
  len <- len - 1
  c1 <- array(0,rep(len,1))
  c2 <- array(0,rep(len,2))
  c3 <- array(0,rep(len,3))
  c4 <- array(0,rep(len,4))
  for (ii in 1:len) {
    for (jj in 1:len) {
      c2[ii, jj] <- m2[ii, jj] - m1[ii]*m1[jj]
      for (kk in 1:len) {
        term21 <- m2[ii, jj]*m1[kk] + m2[ii, kk]*m1[jj] + m2[jj, kk]*m1[ii]
        term111 <- m1[ii] * m1[jj] * m1[kk]
        c3[ii, jj, kk] <- m3[ii, jj, kk] - term21 + 2*term111
        for (ll in 1:len) {
          term31 <- m3[ii, jj, kk]*m1[ll] + m3[ii, jj, ll]*m1[kk] +
            m3[ii, kk, ll]*m1[jj] + m3[jj, kk, ll]*m1[ii]
          term22 <- m2[ii, jj]*m2[kk, ll] + m2[ii, kk]*m2[jj, ll] +
            m2[ii, ll]*m2[jj, kk]
          term211 <- m2[ii, jj]*m1[kk]*m1[ll] + m2[ii, kk]*m1[jj]*m1[ll] +
            m2[ii, ll]*m1[jj]*m1[kk] + m2[jj, kk]*m1[ii]*m1[ll] +
            m2[jj, ll]*m1[ii]*m1[kk] + m2[kk, ll]*m1[ii]*m1[jj]
          term1111 <- m1[ii]*m1[jj]*m1[kk]*m1[ll]
          c4[ii, jj, kk, ll] <- m4[ii, jj, kk, ll] - term31 - term22 + 2*term211 - 6*term1111
        }
      }
    }
  }
  return(list(c1=c1, c2=c2, c3=c3, c4=c4))
}

# Calculating the deltas from Yarnold 1972, 4.4
calc_delta1 <- function(c2, c4){
  delta1<- 0
  len <- dim(c2)[1]
  inv_c2<-solve(apply(c2,2,c))
  for (i_1 in 1:len) {
    for (i_2 in 1:len) {
      for (i_3 in 1:len) {
        for (i_4 in 1:len) {
          delta1<- delta1 + (inv_c2[i_1,i_2] * inv_c2[i_3,i_4] * c4[i_1, i_2, i_3, i_4])
        }
      }
    }
  }
  return(delta1*(1/8))
}
calc_delta2 <- function(c2, c3) {
  delta2<-0
  inv_c2<-solve(apply(c2,2,c))
  len <- dim(c2)[1]
  for (i_1 in 1:len) {
    for (i_2 in 1:len) {
      for (i_3 in 1:len) {
        for (i_4 in 1:len) {
          for (i_5 in 1:len) {
            for (i_6 in 1:len) {
              delta2<- delta2 + (((inv_c2[i_1, i_2] * inv_c2[i_3, i_4] * inv_c2[i_5, i_6])/8 +
                                    (inv_c2[i_1, i_4] * inv_c2[i_2, i_5] * inv_c2[i_3, i_6])/12) *
                                   c3[i_1,i_2,i_3] *c3[i_4, i_5, i_6])
            }
          }
        }
      }
    }
  }
  return(delta2)
}

# Calculate volume of ellipsoid
vol_ellipsoid <- function(m2, n, c, k) {
  m2 <- as.matrix(m2)
  v <- ((pi * n * c)^(k/2) * det(m2)^(1/2))/ gamma(0.5*k + 1)
  return(v)
}

# 4.4 Yarnold Equation Implementation (First Three Terms)
calc_yarnold <- function(nv, c, n = 1) {
  len <- length(nv)
  # Set n = 1 for now.
  k <- len - 1
  moments <- central_moments_calc(nv)
  cumulants <- calc_cumulants(moments$k1, moments$k2, moments$k3, moments$k4)
  avg_rank <- (length(nv) + 1)/2
  #cat("avg_rank", avg_rank, "\n")
  t <- rep(0, len)
  for (i in 1:len){
    t[i] <- nv[i]*avg_rank
  }
  S <- t[2:len]
  cat("S", S, "\n")
  S_means <- rep(0, k)
  for (i in 1:k) {
    S_means[i] <- S[i]/nv[i + 1]
  }
  term1 <- t(S) %*% solve(cumulants$c2) %*% S
  num_per_ellipse <- countl(S_means, k, k, c, chol(cumulants$c2), rep(0, k))
  vol_ellipse <- vol_ellipsoid(apply(moments$m2,2,c), n, c, k)
  term2_ellipse <- num_per_ellipse - vol_ellipse
  term2 <- term2_ellipse * (exp(-c/2)/((2*pi*n)^(k/2) * det(as.matrix(cumulants$c2,2,c))^(1/2)))
  term3_d1 <- (calc_delta1(cumulants$c2, cumulants$c4)/ n) *
    (choose(2, 0)* dchisq(c, k) - choose(2, 1)* dchisq(c, k+2) + choose(2, 2)* dchisq(c, k+4))
  term3_d2 <- (calc_delta2(cumulants$c2, cumulants$c3)/2) * 
    (-choose(3, 0)*dchisq(c, k) + choose(3, 1)*dchisq(c, k+2) - choose(3, 2)*dchisq(c, k+4) +
       choose(3, 3)*dchisq(c, k+6))
  return(term1 + term2 - (term3_d1 + term3_d2))
}

# Yarnold function that takes in nv vector and moments arrays
calc_yarnold_mom <- function(m1, m2, m3, m4, nv, c, n = 1) {
  len <- length(nv)
  k <- len - 1
  cumulants <- calc_cumulants(m1, m2, m3, m4)
  avg_rank <- (length(nv) + 1)/2
  t <- rep(0, len)
  for (i in 1:len){
    t[i] <- nv[i]*avg_rank
  }
  S <- t[2:len]
  cat("S", S, "\n")
  S_means <- rep(0, k)
  for (i in 1:k) {
    S_means[i] <- S[i]/nv[i + 1]
  }
  term1 <- t(S) %*% solve(cumulants$c2) %*% S
  num_per_ellipse <- countl(S_means, k, k, c, chol(cumulants$c2), rep(0, k))
  vol_ellipse <- vol_ellipsoid(apply(m2,2,c), n, c, k)
  term2_ellipse <- num_per_ellipse - vol_ellipse
  term2 <- term2_ellipse * (exp(-c/2)/((2*pi*n)^(k/2) * det(as.matrix(cumulants$c2,2,c))^(1/2)))
  term3_d1 <- (calc_delta1(cumulants$c2, cumulants$c4)/ n) *
    (choose(2, 0)* dchisq(c, k) - choose(2, 1)* dchisq(c, k+2) + choose(2, 2)* dchisq(c, k+4))
  term3_d2 <- (calc_delta2(cumulants$c2, cumulants$c3)/2) * 
    (-choose(3, 0)*dchisq(c, k) + choose(3, 1)*dchisq(c, k+2) - choose(3, 2)*dchisq(c, k+4) +
       choose(3, 3)*dchisq(c, k+6))
  return(term1 + term2 - (term3_d1 + term3_d2))
}

# Testing...


