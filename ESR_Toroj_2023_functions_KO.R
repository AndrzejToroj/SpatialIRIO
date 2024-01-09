#FUNCTIONS

level2trans <- function(theta) {
  if(drivetime == FALSE) {
    theta_t <- c(log(theta[1:(2*N_S_reg)]),
                 qnorm((theta[(2*N_S_reg+1):(2*N_S_reg+N_S_reg*(N_S_reg-1)/2)]+1)/2),
                 log(theta[(2*N_S_reg+N_S_reg*(N_S_reg-1)/2+1):(3*N_S_reg+N_S_reg*(N_S_reg-1)/2)]))
  } else {
    theta_t <- c(log(theta[1:(2*N_S_reg)]),
                 theta[(2*N_S_reg+1):(3*N_S_reg)],
                 qnorm((theta[(3*N_S_reg+1):(3*N_S_reg+N_S_reg*(N_S_reg-1)/2)]+1)/2),
                 log(theta[(3*N_S_reg+N_S_reg*(N_S_reg-1)/2+1):(4*N_S_reg+N_S_reg*(N_S_reg-1)/2)]))
  }
  return(theta_t)
}

trans2level <- function(theta_t) {
  if(drivetime == FALSE) {
    theta <- c(exp(theta_t[1:(2*N_S_reg)]),
               pnorm(theta_t[(2*N_S_reg+1):(2*N_S_reg+N_S_reg*(N_S_reg-1)/2)])*2-1,
               exp(theta_t[(2*N_S_reg+N_S_reg*(N_S_reg-1)/2+1):(3*N_S_reg+N_S_reg*(N_S_reg-1)/2)]))
  } else {
    theta <- c(exp(theta_t[1:(2*N_S_reg)]),
               theta_t[(2*N_S_reg+1):(3*N_S_reg)],
               pnorm(theta_t[(3*N_S_reg+1):(3*N_S_reg+N_S_reg*(N_S_reg-1)/2)])*2-1,
               exp(theta_t[(3*N_S_reg+N_S_reg*(N_S_reg-1)/2+1):(4*N_S_reg+N_S_reg*(N_S_reg-1)/2)]))
  }
  return(theta)
}

expected_scale <- function(shape) {
  l_shape <- log(shape)
  X <- c(1, l_shape, l_shape^2, l_shape^3, l_shape^4, l_shape^5, l_shape^6, l_shape^7, l_shape^8, l_shape^9, l_shape^10, l_shape^11)
  l_scale <- rep(NA, N_S_reg)
  for (sector in 1:N_S_reg) {
    l_scale[sector] <- sum(X * prior_builder$EV[[sectors$prior_sector[sector]]]) #KO1
  }
  return(exp(l_scale))
}

#Prior: one sector parameters' contribution to log-density
log_prior_1S_PL <- function(sector, shape, scale) {
  l_shape <- shape
  logD_shape <- dunif(l_shape, min = minlsh, max = maxlsh, log = TRUE)
  if (logD_shape == -Inf) {
    logD_prior <- -Inf
  } else {
    X <- c(1, l_shape, l_shape^2, l_shape^3, l_shape^4, l_shape^5, l_shape^6, l_shape^7, l_shape^8, l_shape^9, l_shape^10, l_shape^11)
    first <- sum(X * prior_builder$first)
    last <- sum(X * prior_builder$last) 
    EV.ss <- sum(X * prior_builder$EV[[sectors$prior_sector[sector]]])#KO2
    SD.ss <- sum(X * prior_builder$SD[[sectors$prior_sector[sector]]])#KO3
    l_scale <- scale
    D_scale <- dtruncnorm(l_scale, a=first, b=last, mean=EV.ss, sd=SD.ss)
    if (D_scale == 0) {
      logD_prior <- -Inf
    } else {
      logD_prior <- logD_shape + log(D_scale)
    }
  }
  return(logD_prior)
}

#Prior log-density
log_prior <- function(theta, hyper) {
  
  shapes <- theta[1:N_S_reg]
  scales <- theta[(N_S_reg+1):(2*N_S_reg)]
  variances <- exp(theta[(length(theta)-N_S_reg+1):length(theta)])
  correlations <- pnorm(theta[(2*N_S_reg+1):(length(theta)-N_S_reg)])*2 - 1

  Sd.small <- variances^0.5
  Corr.small <- diag(N_S_reg)
  Corr.small[upper.tri(Corr.small, diag = FALSE)] <- correlations
  Corr.small <- Corr.small + t(Corr.small) - diag(N_S_reg)
  Sigma.small <- diag(Sd.small) %*% Corr.small %*% diag(Sd.small)
  Sigma.small <- round(Sigma.small, 3)
  if (is.symmetric.matrix(Sigma.small)) {
    Sigma.small <- Sigma.small
  } else {
    Sigma.small <- round(Sigma.small, 2)
    if (is.symmetric.matrix(Sigma.small)) {
      Sigma.small <- Sigma.small
    } else {
      Sigma.small <- forceSymmetric(round(Sigma.small, 1))
      Sigma.small <- matrix(Sigma.small@x, nrow = N_S_reg)
    }
  }
  
  lp_sector <- rep(0, N_S_reg)
  for (sector in 1:N_S_reg) {
    lp_sector[sector] <- log_prior_1S_PL(sector, shapes[sector], scales[sector])
  }
  
  hyper_var <- hyper[1:N_S_reg]
  hyper_cor <- hyper[N_S_reg+1]
  
  df <- N_S_reg
  H <- diag(length(hyper_var))
  H[upper.tri(H, diag = FALSE)] <- hyper_cor
  H <- H + t(H) - diag(length(hyper_var))
  hyper_sd <- hyper_var^0.5
  H <- diag(hyper_sd) %*% H %*% diag(hyper_sd) * (df+N_S_reg+1)
  H <- round(H, 7)
  H <- forceSymmetric(H)
  H <- matrix(H@x, N_S_reg, N_S_reg)
  
  if (is.positive.semi.definite(Sigma.small)) {
    scaling2 <- 1 #1/mean(Sigma.small) #10^(-12)*1.95
    lp_stochastic <- log(diwish(scaling2*Sigma.small, df, scaling2*H))
  } else {
    lp_stochastic <- -Inf
  }
  
  logD_prior <- sum(lp_sector) + lp_stochastic
  return(logD_prior)
}


frames <- function(theta, dist, va_reg) {

  vertical.frame <- list()
  diagonal.frame <- list()
  square.frame <- list()
  
  W <- construct_W_est(theta, dist, va_reg)
  n <- N_S_reg * N_R
  
  for (yy in 1:N_T) {
    vertical.frame[[yy]] <- matrix(NA, nrow = n, ncol = N_R)
    diagonal.frame[[yy]] <- matrix(0, nrow = n, ncol = n)
    for(ss in 1:N_S_reg) {
      vertical.frame[[yy]][((ss-1)*N_R+1):(ss*N_R), 1:N_R] <- as.matrix(W$W[[ss]])
      diagonal.frame[[yy]][((ss-1)*N_R+1):(ss*N_R), ((ss-1)*N_R+1):(ss*N_R)] <- as.matrix(W$W[[ss]])
    }
    square.frame[[yy]] <- kronecker(t(as.matrix(rep(1,N_S_reg))), vertical.frame[[yy]])
  }
  
  frames <- list()
  frames$diagonal.frame <- diagonal.frame
  frames$vertical.frame <- vertical.frame
  frames$square.frame <- square.frame
  return(frames)
  
}

construct_W_est <- function(theta, dist, va_reg) {
  
  W2 <- list()
  #W <- list()
  if (length(dim(va_reg)) > 2) {
    years <- dim(va_reg)[3]
  } else {
    years <- 1
  }
  
  W1 <- lapply(1:N_S_reg, construct_W1, regions = N_R, sectors = N_S_reg, dist = dist, theta = theta)
  for (yy in 1:years) {
    if (years > 1) {
      VA_vector <- as.matrix(c(va_reg[,,yy]))
    } else {
      VA_vector <- as.matrix(c(va_reg))
    }
    W2[[yy]] <- lapply(1:N_S_reg, construct_W2, regions = N_R, VA_vector = VA_vector)
  }
  
  WW <- list(W=W1,W2=W2)
  return(WW)
}

construct.B <- function(beta0, regions, diagonal.frame) {
  if (length(dim(beta0)) > 2) {
    years.est <- dim(beta0)[3]
  } else {
    years.est <- 1
  }
  B <- list()
  for (yy in 1:years.est) {
    if (years.est > 1) {
      b0 <- beta0[, , yy]
    } else {
      b0 <- beta0
    }
    B[[yy]] <- kronecker(diag(as.vector(b0)), matrix(1, nrow = regions, ncol = regions)) * diagonal.frame[[yy]]
  }
  return(B)
}

construct.A <- function(beta, regions, square.frame) {
  if (length(dim(beta)) > 2) {
    years.est <- dim(beta)[3]
  } else {
    years.est <- 1
  }
  A <- list()
  for (yy in 1:years.est) {
    if (years.est > 1) {
      b <- beta[, , yy]
    } else {
      b <- beta[,,1]
    }
    A[[yy]] <- kronecker(as.matrix(b), matrix(1, nrow = regions, ncol = regions)) * square.frame[[yy]]
    A[[yy]] <- diag(nrow(A[[yy]])) - as.matrix(A[[yy]])
  }
  return(A)
}

construct_W1 <- function(jj, regions, sectors, dist, theta) {
  dist.utri.v <- dist[upper.tri(dist, diag = FALSE)]
  W1.v <- rep(1, length(dist.utri.v)) - sapply(X = dist.utri.v, FUN = pgamma, 
                                               shape = theta[jj], scale = theta[sectors + jj])
  W1 <- diag(regions)
  W1[upper.tri(W1, diag = FALSE)] <- W1.v
  W1 <- W1 + t(W1) - diag(regions)
  sum.col <- kronecker(as.matrix(rep(1, regions)), t(as.matrix(rowSums(W1))))
  W1 <- W1 / sum.col
  return(W1)
}

construct_W2 <- function(jj, regions, VA_vector) {
  W2 <- kronecker(t(as.matrix(rep(1,regions))), VA_vector[((jj-1)*regions+1):(jj*regions)])
  W2[W2 < 0] <- 0
  sum.col <- kronecker(as.matrix(rep(1,regions)), t(as.matrix(colSums(W2))))
  W2 <- W2 / sum.col
  return(W2)
}

#Log-likelihood
loglik <- function(theta, va_reg, y, beta, beta0, dist) {
  
  theta <- trans2level(theta)
  if (length(dim(va_reg)) > 2) {
    years <- dim(va_reg)[3]
  } else {
    years <- 1
  }
  n <- N_S_reg * N_R
  
  Sd.small <- diag(theta[(length(theta) - (N_S_reg - 1)):length(theta)]) ^ 0.5
  Corr.small <- diag(N_S_reg)
  Corr.small[upper.tri(Corr.small, diag = FALSE)] <- theta[(length(theta) - ((N_S_reg*(N_S_reg-1)/2)-1)-N_S_reg):(length(theta) - N_S_reg)]
  Corr.small <- Corr.small + t(Corr.small) - diag(N_S_reg)
  Sigma.small <- Sd.small %*% Corr.small %*% Sd.small
  
  inv.Corr.small <- tryCatch({
    solve(Corr.small)
  }, error = function(err) {
    matrix(NA, nrow = N_S_reg, ncol = N_S_reg)
  })
  
  inv.Sd.small <- diag(1 / diag(Sd.small))
  det.Sigma.small <- det(Sigma.small)
  if (det.Sigma.small <= 0) {
    return(-Inf)
  } else {
  
    if (sum(as.integer(is.nan(Sd.small))) > 0) {
      det.Sigma.small <- -1
      inv.Sd.small <- diag(c(Inf, N_S_reg))
    }
    
    inv.Sigma.small <- inv.Sd.small %*% inv.Corr.small %*% inv.Sd.small
    Sigma <- kronecker(Sigma.small, diag(N_R))
    inv.Sigma <- kronecker(inv.Sigma.small, diag(N_R))
    
    VA_vector <- list()
    for (yy in 1:years) {
      if (years > 1) {
        VA_vector[[yy]] <- as.matrix(c(va_reg[,,yy]))
      } else {
        VA_vector[[yy]] <- as.matrix(c(va_reg))
      }
    }
    
    frames.list <- frames(theta, dist, va_reg)
    B <- construct.B(beta0, N_R, frames.list$diagonal.frame)
    A <- construct.A(beta, N_R, frames.list$square.frame)
    
    ll <- 0
    for (yy in 1:years){
      ll <- ll - n / 2 * log(2 * pi) - (1 / 2) * (t(A[[yy]] %*% as.matrix(VA_vector[[yy]]) - B[[yy]] %*% y[[yy]]) %*% inv.Sigma %*% (A[[yy]] %*% as.matrix(VA_vector[[yy]]) - B[[yy]] %*% y[[yy]])) + log(abs(det(A[[yy]]))) - ((1 / 2) * N_R * log(det.Sigma.small))
    }
    #ll <- -ll
    return(ll)
  }
}

log_posterior <- function(theta, hyper, data) {
  LL <- loglik(theta, data$va_reg, data$y, data$beta, data$beta0, data$dist)
  prior <- log_prior(theta, hyper)
  return(LL+prior)
}

proposalfunction <- function(param, scale){
  relatives <- diag(c(rep(5,N_S_reg), rep(25, N_S_reg), rep(2, N_S_reg*(N_S_reg-1)/2), rep(25,N_S_reg)) / 1000)
  return(rmvnorm(1, mean = param, sigma = scale * relatives))
}

run_metropolis_MCMC <- function(theta0 = theta0, S = S, scale = scale){
  chain = matrix(NA, nrow = S+1, ncol = length(theta0))
  chain[1, ] = theta0
  for (ss in 2:(S+1)){
    print(paste0("Iteration ", ss))
    proposal <- proposalfunction(chain[ss-1, ], scale)
    previous <- t(as.matrix(chain[ss-1, ]))
    probab = exp(log_posterior(proposal, hyper, data) - log_posterior(previous, hyper, data))
    if (is.na(probab) | is.nan (probab) | !is.numeric(probab)) {
      print(log_posterior(proposal, hyper, data))
      print(log_posterior(previous, hyper, data))
    }
    if (runif(1) < probab){
      chain[ss, ] <- proposal
    } else {
      chain[ss, ] <- previous
    }
  }
  return(chain)
}

draw_prior_posterior_all <- function(minlsh, maxlsh, prior_builder, results_list, theta0){
  png(filename = paste0("output/priors_posteriors_", mode, country, ".png"), width = 2000, height = 1000)
  par(mfrow = c(2, N_S_reg), mar = c(4, 4, 4, 4), omd = c(0.05, 0.95, 0.05, 0.95), xpd = FALSE)
  colors <- c(rgb(169,169,169,100, maxColorValue = 255), rgb(0,139,139,100, maxColorValue = 255), "darkred", rgb(0, 0, 1, 0.5))
  temp <- rbind(results_list[[1]], results_list[[2]])
  results_list <- list()
  results_list[[1]] <- temp

  maxYshape <- 0
  minXshape <- minlsh
  maxXshape <- maxlsh
  for(ss in 1:N_S_reg) {
    temp_hist <- list()
    for (ii in 1:length(results_list)) {
      temp_hist[[ii]] <- hist(results_list[[ii]][(burnIn+1):S,ss], breaks = 30, plot = FALSE)
      maxYshape <- max(c(maxYshape, temp_hist[[ii]]$density))
      minXshape <- min(c(minXshape, min(results_list[[ii]][(burnIn+1):S,ss])))
      maxXshape <- max(c(maxXshape, max(results_list[[ii]][(burnIn+1):S,ss])))
    }
  }
  #### by hand!#
  ##############
  maxXshape <- 3
  ##############

  for (ss in 1:N_S_reg) {
    draw_pp_shape_1s(ss, minlsh, maxlsh, results_list, colors, theta0, maxYshape, minXshape, maxXshape)
  }
  
  maxYscale <- 0
  temp_prior <- list()
  minXscale <- 0
  maxXscale <- 0

  for (ss in 1:N_S_reg) {
    temp_prior[[ss]] <- marginal_prior_scale(ss)
    minXscale <- min(c(minXscale, temp_prior[[ss]]))
    maxXscale <- max(c(maxXscale, temp_prior[[ss]]))
    temp_prior[[ss]] <- hist(temp_prior[[ss]], breaks = 30, plot = FALSE)
    maxYscale <- max(maxYscale, max(temp_prior[[ss]]$density))
    temp_hist <- list()
    for (ii in 1:length(results_list)) {
      temp_hist[[ii]] <- hist(results_list[[ii]][(burnIn+1):S,ss+N_S_reg], breaks = 30, plot = FALSE)
      maxYscale <- max(c(maxYscale, temp_hist[[ii]]$density))
      minXscale <- min(c(minXscale, min(results_list[[ii]][(burnIn+1):S,ss+N_S_reg])))
      maxXscale <- max(c(maxXscale, max(results_list[[ii]][(burnIn+1):S,ss+N_S_reg])))
    }
  }
  #### by hand!#
  ##############
  minXscale <- 0
  ##############
  
  for (ss in 1:N_S_reg) {
    draw_pp_scale_1s(ss, prior_builder, results_list, colors, theta0, maxYscale, minXscale, maxXscale)
  }
  dev.off()
}

draw_pp_shape_1s <- function(ss, minlsh, maxlsh, results_list, colors, theta0, maxYshape, minXshape, maxXshape) {
  domain <- seq(from = minlsh, to = maxlsh, length.out = 1000)
  temp_hist <- list()
  maxY <- maxYshape#0
  minX <- minXshape
  maxX <- maxXshape
  for (ii in 1:length(results_list)) {
    temp_hist[[ii]] <- hist(results_list[[ii]][(burnIn+1):S,ss], breaks = 30, plot = FALSE)
  }
  plot(x = domain, y = dunif(domain, min = minlsh, max = maxlsh, log = FALSE), col = colors[1], type = "l", lwd = 3,
       ylim = c(0, 0.7*maxY), xlim = c(minX, maxX), ylab = "", xlab = "", main = paste0("log of Shape, sec. ", codes_S_reg[ss]),
       cex.main = 2.5, cex.axis=2.5, cex.lab = 2.5)
  polygon(x = c(domain, rev(domain)), y = c(dunif(domain, min = minlsh, max = maxlsh, log = FALSE), rep(0, length(domain))), col = colors[1], border = NA)
  for (ii in 1:length(temp_hist)) {
    lo <- loess(temp_hist[[ii]]$density ~ temp_hist[[ii]]$mids, span = 0.65)
    polygon(x = c(temp_hist[[ii]]$mids, rev(temp_hist[[ii]]$mids)), y = c(predict(lo), rep(0, length(temp_hist[[ii]]$mids))), col = colors[ii+1], border = NA)
    abline(v = theta0[ss], col = colors[ii+1])
  }
}

draw_pp_scale_1s <- function(ss, prior_builder, results_list, colors, initial_pts_shape0, maxYscale, minXscale, maxXscale) {
  temp_prior <- marginal_prior_scale(ss)
  temp_hist <- list()
  temp_hist[[1]] <- hist(temp_prior, breaks = 30, plot = FALSE)
  maxY <- maxYscale
  minX <- minXscale
  maxX <- maxXscale
  for (ii in 1:length(results_list)) {
    temp_hist[[ii+1]] <- hist(results_list[[ii]][(burnIn+1):S,ss+7], breaks = 30, plot = FALSE)
  }
  
  plot(x = temp_hist[[1]]$mids, y = temp_hist[[1]]$density, type = "l", col = colors[1], 
       ylim = c(0, 0.85*maxY), xlim = c(minX, maxX), ylab = "", xlab = "", lwd = 3,
       main = paste0("log of Scale, sec. ", codes_S_reg[ss]), cex.main = 2.5, cex.axis=2.5, cex.lab = 2.5)
  polygon(x = c(temp_hist[[ii]]$mids, rev(temp_hist[[ii]]$mids)), y = c(temp_hist[[1]]$density, rep(0, length(temp_hist[[1]]$density))), col = colors[ii], border = NA)
  for (ii in 2:length(temp_hist)) {
    lo <- loess(temp_hist[[ii]]$density ~ temp_hist[[ii]]$mids, span = 0.65)
    polygon(x = c(temp_hist[[ii]]$mids, rev(temp_hist[[ii]]$mids)), y = c(predict(lo), rep(0, length(temp_hist[[ii]]$mids))), col = colors[ii], border = NA)
    abline(v = theta0[ss+7], col = colors[ii])
  }
}

marginal_prior_scale <- function(ss) {
  NN <- 100000
  l_shape <- runif(NN, min = minlsh, max = maxlsh)
  X <- as.matrix(data.frame(rep(1,NN), l_shape, l_shape^2, l_shape^3, l_shape^4, l_shape^5, l_shape^6, l_shape^7, l_shape^8, l_shape^9, l_shape^10, l_shape^11))
  first <- X %*% as.matrix(prior_builder$first)
  last <- X %*% as.matrix(prior_builder$last)
  EV.ss <- X %*% as.matrix(prior_builder$EV[[sectors$prior_sector[ss]]]) #KO4
  SD.ss <- X %*% as.matrix(prior_builder$SD[[sectors$prior_sector[ss]]]) #KO5
  prior_scale <- rtruncnorm(1, a=first, b=last, mean=EV.ss, sd=SD.ss)
  prior_scale <- prior_scale[!is.na(prior_scale)]
  return(prior_scale)
}

log_posterior_wrap <- function(theta0, data2) {
  data <- data2$data
  hyper <- data2$hyper
  return(log_posterior(theta0, hyper, data))
}

RAS_blocks <- function(x_R, X_R, VA_Sio, stacked_VA_shares, N_R, N_S_io, x_Sio_times_Sio, proj_CS, proj_RS) {
  
  VA_NUTS3 <- rep(as.vector(VA_Sio), each = N_R) * stacked_VA_shares
  
  x_temp <- x_R
  A <- x_R / kronecker(as.matrix(rep(1,length(X_R))), X_R)
  A[is.nan(A)] <- 0
  RABSB <- A
  stop_criterion <- 0.1
  iter <- 1
  print("Starting blocksum-augmented tri-proportional RAS procedure")
  
  while (stop_criterion > 10^(-2)) {
    
    r <- proj_RS / rowSums(x_temp)
    r[is.nan(r)] <- 0
    #R <- diag(as.vector(r))
    #RA <- R %*% RABSB
    RA <- kronecker(t(as.matrix(rep(1,N_R*N_S_io))), as.matrix(r)) * RABSB #quicker version
    x_temp <- RA * kronecker(as.matrix(rep(1,length(X_R))), X_R)
    
    stop_criterion_1 <- max(abs(r[r!=0]-1))
    
    xs_temp <- matrix(NA, nrow = N_S_io, ncol = N_S_io)
    for (ii in 1:N_S_io) {
      for (jj in 1:N_S_io) {
        xs_temp[ii,jj] <- sum(x_temp[((ii-1)*N_R+1):(ii*N_R), ((jj-1)*N_R+1):(jj*N_R)])
      }
    }
    multipliers <- kronecker((x_Sio_times_Sio / xs_temp), matrix(1, nrow = N_R, ncol = N_R))
    multipliers[is.nan(multipliers)] <- 0
    x_temp <- x_temp * multipliers
    RAB <- x_temp / kronecker(as.matrix(rep(1,length(X_R))), X_R)
    RAB[is.nan(RAB)] <- 0
    
    stop_criterion_2 <- max(abs(multipliers[multipliers!=0]-1))
    
    s <- proj_CS / colSums(x_temp)
    s[is.nan(s)] <- 0
    #S <- diag(as.vector(s))
    #RABS <- RAB %*% S
    RABS <- RAB * kronecker(as.matrix(rep(1,N_R*N_S_io)), t(as.matrix(s)))
    x_temp <- RABS * kronecker(as.matrix(rep(1,length(X_R))), X_R)
    
    stop_criterion_3 <- max(abs(s[s!=0]-1))
    
    xs_temp <- matrix(NA, nrow = N_S_io, ncol = N_S_io)
    for (ii in 1:N_S_io) {
      for (jj in 1:N_S_io) {
        xs_temp[ii,jj] <- sum(x_temp[((ii-1)*N_R+1):(ii*N_R), ((jj-1)*N_R+1):(jj*N_R)])
      }
    }
    multipliers <- kronecker((x_Sio_times_Sio / xs_temp), matrix(1, nrow = N_R, ncol = N_R))
    multipliers[is.nan(multipliers)] <- 0
    x_temp <- x_temp * multipliers
    RABSB <- x_temp / kronecker(as.matrix(rep(1,length(X_R))), X_R)
    RABSB[is.nan(RABSB)] <- 0
    
    stop_criterion_4 <- max(abs(multipliers[multipliers!=0]-1))
    print(paste0(iter, " | ", round(stop_criterion_1,4), " | ", round(stop_criterion_2,4), " | ", round(stop_criterion_3,4), " | ", round(stop_criterion_4,4)))
    stop_criterion <- max(stop_criterion_1, stop_criterion_2, stop_criterion_3, stop_criterion_4)
    iter <- iter + 1
    
  }
  
  x_NUTS3_proj <- x_temp
  return(x_NUTS3_proj)
  
}

MRIO <- function(theta) {
  W <- construct_W_est(theta, dist, va_reg[,,N_T])
  W_ras0 <- list()
  for (ss in 1:N_S_reg) {
    W_ras0[[ss]] <- (W$W[[ss]]) * (W$W2[[1]][[ss]] ^ 0.98)
    scaling <- kronecker(matrix(1, nrow = N_R, ncol = 1),
                         t(as.matrix(colSums(W_ras0[[ss]]))))
    W_ras0[[ss]] <- W_ras0[[ss]] / scaling
  }
  stacked_VA_shares <- matrix(NA,
                              nrow = 0,
                              ncol = 1)
  stacked_W <- matrix(NA,
                      nrow = 0,
                      ncol = N_R)
  for(ss in 1:N_S_io){
    stacked_VA_shares <- rbind(stacked_VA_shares,
                               as.matrix(W$W2[[1]][[sectors$id_vareg[ss]]][,1]))
    stacked_W <- rbind(stacked_W,
                       W_ras0[[sectors$id_vareg[ss]]])
  }
  X_R <- t(stacked_VA_shares) * kronecker(t(sectors$X), matrix(1, nrow = 1, ncol = N_R))
  interpolation_grid <- kronecker(matrix(1, nrow=N_S_io*N_R, 1), t(stacked_VA_shares)) * kronecker(matrix(1, nrow=1, ncol=N_S_io), stacked_W)
  x_R <- interpolation_grid * kronecker(x_Sio_times_Sio, matrix(1, nrow = N_R, ncol = N_R))
  rm(interpolation_grid)
  
  #3-RAS balanced MRIO
  proj_CS <- colSums(x_R)
  if (country == "PL") {
    proj_RS <- rep(rowSums(x_Sio_times_Sio), each = N_R) * stacked_VA_shares #t(X_R) - kronecker(as.matrix(sectors$Y), as.matrix(Y_shares[,ncol(Y_shares)]))#
  }
  if (country == "KO") {
    aux <- sectors_SR[,c("region_code", "sector_code", "Y", "X")]
    aux <- aux[order(aux$sector_code), ]
    proj_RS <- aux$X - aux$Y
  }
  x_R <- RAS_blocks(x_R, X_R, VA_Sio, stacked_VA_shares, N_R, N_S_io, x_Sio_times_Sio, proj_CS, proj_RS)
  
  #Create header and save
  x_R_header <- data.frame(sector_id = rep(sectors$id, each = N_R),
                           sector = rep(sectors$name_S_io, each = N_R),
                           region = rep(names_R, N_S_io),
                           stringsAsFactors = FALSE)
  save(x_R, x_R_header, X_R, file = paste0("MRIO_", country, ".RData"))
  
  if(country == "KO") {
    print("Matrix x saved")
    return(NULL)
  }
  if(country == "PL") {
    
    A_R <- x_R / kronecker(matrix(1, nrow = N_S_io * N_R, ncol = 1), X_R)
    A_R[is.nan(A_R)] <- 0
    A_R_ind_block_12 <- kronecker(as.matrix(sectors$induced_col), matrix(1, nrow = N_R, ncol = N_R)) * stacked_W
    all_local_induced <- c(6, 37) #6 - food, 37 - retail trade
    for(ll in 1:length(all_local_induced)) {
      A_R_ind_block_12[((all_local_induced[ll]-1) * N_R + 1) : (all_local_induced[ll] * N_R), ] <- diag(N_R) * sectors$induced_col[all_local_induced[ll]]
    }
    rm(all_local_induced, ll)
    
    if (useDatarino == FALSE) {
      block21maker <- diag(N_R)
    } else {
      block21maker <- as.matrix(read.csv( file = 'commuting.csv', 
                                          sep = ";",
                                          dec = ".",
                                          header = FALSE,
                                          stringsAsFactors = FALSE))
    }
    A_R_ind_block_21 <- kronecker(t(as.matrix(sectors$induced_row)), block21maker)
    A_R_ind_block_22 <- matrix(0, nrow = N_R, ncol = N_R)
    A_R_ind <- rbind(cbind(A_R, A_R_ind_block_12),cbind(A_R_ind_block_21, A_R_ind_block_22))
    #rm(A_R_ind_block_12, A_R_ind_block_21, A_R_ind_block_22)
    return(A_R_ind)
  }
}
