#FUNCTIONS

#Prior: one sector parameters' contribution to log-density
log_prior_1S_PL <- function(sector, shape, scale, alpha) {
  l_shape <- shape
  logD_shape <- dunif(l_shape, min = minlsh, max = maxlsh, log = TRUE)
  logD_alpha <- dunif(alpha, min = 0, max = 1, log = TRUE)
  if (logD_shape == -Inf | logD_alpha == -Inf) {
    logD_prior <- -Inf
  } else {
    X <- c(1, l_shape, l_shape^2, l_shape^3, l_shape^4, l_shape^5, l_shape^6, l_shape^7, l_shape^8, l_shape^9, l_shape^10, l_shape^11)
    first <- sum(X * prior_builder$first)
    last <- sum(X * prior_builder$last) 
    EV.ss <- sum(X * prior_builder$EV[[sector]])
    SD.ss <- sum(X * prior_builder$SD[[sector]])
    l_scale <- scale
    D_scale <- dtruncnorm(l_scale, a=first, b=last, mean=EV.ss, sd=SD.ss)
    if (D_scale == 0) {
      logD_prior <- -Inf
    } else {
      logD_prior <- logD_shape + log(D_scale) + logD_alpha
    }
  }
  return(logD_prior)
}

#Prior log-density
log_prior <- function(theta, hyper) {

  shapes <- theta[1:N_S_reg]
  scales <- theta[(N_S_reg+1):(2*N_S_reg)]
  alphas <- theta[(2*N_S_reg+1):(3*N_S_reg)]
  variances <- exp(theta[(length(theta)-N_S_reg+1):length(theta)])
  correlations <- pnorm(theta[(3*N_S_reg+1):(length(theta)-N_S_reg)])*2 - 1
  
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
      Sigma.small <- round(Sigma.small, 1)
    }
  }
  
  lp_sector <- rep(0, N_S_reg)
  for (sector in 1:N_S_reg) {
    lp_sector[sector] <- log_prior_1S_PL(sector, shapes[sector], scales[sector], alphas[sector])
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
    lp_stochastic <- log(diwish(Sigma.small, df, H))
  } else {
    lp_stochastic <- -Inf
  }
  
  logD_prior <- sum(lp_sector) + lp_stochastic
  return(logD_prior)
}


frames <- function(theta, dist_r, dist_t, va_reg) {

  vertical.frame <- list()
  diagonal.frame <- list()
  square.frame <- list()
  
  W <- construct_W_est(theta, dist_r, dist_t, va_reg)
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

construct_W_est <- function(theta, dist_r, dist_t, va_reg) {
  
  W2 <- list()
  #W <- list()
  if (length(dim(va_reg)) > 2) {
    years <- dim(va_reg)[3]
  } else {
    years <- 1
  }
  
  W1 <- lapply(1:N_S_reg, construct_W1, regions = N_R, sectors = N_S_reg, dist_r = dist_r, dist_t = dist_t, theta = theta)
  for (yy in 1:years) {
    if (years > 1) {
      VA_vector <- as.matrix(c(va_reg[,,yy]))
    } else {
      VA_vector <- as.matrix(c(va_reg))
    }
    W2[[yy]] <- lapply(1:N_S_reg, construct_W2, regions = N_R, VA_vector = VA_vector)
    #W[[yy]] <- lapply(1:N_S_reg, aggregate_W, W1 = W1, W2 = W2[[yy]], theta = theta, regions = N_R)
  }
  
  
  WW <- list(W=W1,W2=W2)
  return(WW)
}

construct_W1 <- function(jj, regions, sectors, dist_r, dist_t, theta) {
  dist.r.utri.v <- dist_r[upper.tri(dist_r, diag = FALSE)]
  dist.t.utri.v <- dist_t[upper.tri(dist_t, diag = FALSE)]
  dist.utri.v.ratio <- dist.t.utri.v / dist.r.utri.v
  dist.utri.v <- dist.r.utri.v * (dist.utri.v.ratio^theta[2*sectors+jj])
  W1.v <- rep(1, length(dist.utri.v)) - sapply(X = dist.utri.v, FUN = pgamma, 
                                               shape = theta[jj], scale = theta[sectors + jj])
  W1 <- diag(regions)
  W1[upper.tri(W1, diag = FALSE)] <- W1.v
  W1 <- W1 + t(W1) - diag(regions)
  sum.col <- kronecker(as.matrix(rep(1, regions)), t(as.matrix(rowSums(W1))))
  W1 <- W1 / sum.col
  return(W1)
}

#Log-likelihood
loglik <- function(theta, va_reg, y, beta, beta0, dist_r, dist_t) {
  
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
    
    frames.list <- frames(theta, dist_r, dist_t, va_reg)
    B <- construct.B(beta0, N_R, frames.list$diagonal.frame)
    A <- construct.A(beta, N_R, frames.list$square.frame)
    
    ll <- 0
    for (yy in 1:years){
      ll <- ll - n / 2 * log(2 * pi) - (1 / 2) * (t(A[[yy]] %*% as.matrix(VA_vector[[yy]]) - B[[yy]] %*% y[[yy]]) %*% inv.Sigma %*% (A[[yy]] %*% as.matrix(VA_vector[[yy]]) - B[[yy]] %*% y[[yy]])) + log(abs(det(A[[yy]]))) - ((1 / 2) * N_R * log(det.Sigma.small))
    }
    return(ll)
  }
}

log_posterior <- function(theta, hyper, data) {
  LL <- loglik(theta, data$va_reg, data$y, data$beta, data$beta0, data$dist_r, data$dist_t)
  prior <- log_prior(theta, hyper)
  return(LL+prior)
}

proposalfunction <- function(param, scale){
  relatives <- diag(c(rep(5,N_S_reg), rep(25, N_S_reg), rep(10, N_S_reg), rep(2, N_S_reg*(N_S_reg-1)/2), rep(25,N_S_reg)) / 1000)
  return(rmvnorm(1, mean = param, sigma = scale * relatives))
}

draw_prior_posterior_time <- function(results_list){
  png(filename = paste0("output/priors_posteriors_gamma_", mode, country,".png"), width = 2000, height = 500)
  par(mfrow = c(1, N_S_reg), mar = c(4, 4, 4, 4), omd = c(0.05, 0.95, 0.05, 0.95), xpd = FALSE)
  colors <- c(rgb(169,169,169,100, maxColorValue = 255), rgb(0,139,139,100, maxColorValue = 255), "darkred", rgb(0, 0, 1, 0.5))
  temp <- rbind(results_list[[1]], results_list[[2]])
  results_list <- list()
  results_list[[1]] <- temp
  theta0 <- rep(0.5, N_S_reg)
  
  maxYshape <- 0
  minXshape <- 0.5
  maxXshape <- 0.5

  for(ss in 1:N_S_reg) {
    for (ii in 1:length(results_list)) {
      temp_hist <- list()
      temp_hist[[ii]] <- hist(results_list[[ii]][,ss+(2*N_S_reg)], breaks = 30, plot = FALSE)
      maxYshape <- max(c(maxYshape, temp_hist[[ii]]$density))
      minXshape <- min(c(minXshape, min(results_list[[ii]][,ss+(2*N_S_reg)])))
      maxXshape <- max(c(maxXshape, max(results_list[[ii]][,ss+(2*N_S_reg)])))
    }
  }
  #### by hand!#
  ##############
  maxYshape <- 3000
  ##############
  
  for (ss in 1:N_S_reg) {
    draw_pp_time_1s(ss, results_list, colors, theta0, maxYshape, minXshape, maxXshape)
  }
  dev.off()
}

draw_pp_time_1s <- function(ss, results_list, colors, theta0, maxYshape, minXshape, maxXshape) {
  domain <- seq(from = 0, to = 1, length.out = 1000)
  temp_hist <- list()
  maxY <- maxYshape

  for (ii in 1:length(results_list)) {
    temp_hist[[ii]] <- hist(results_list[[ii]][(burnIn+1):S,ss+(2*N_S_reg)], breaks = 30, plot = FALSE)
  }
  if (length(temp_hist[[1]]$density) == 1) {
    minX <- 0.99
    maxX <- 1.001
    maxY <- 2000
    temp_hist[[1]]$density <- 2000
    temp_hist[[1]]$mids <- 1
  } else {
    minX <- min(temp_hist[[1]]$mids)
    maxX <- max(temp_hist[[1]]$mids)
    maxY <- max(temp_hist[[1]]$density)
  }
  
  plot(x = domain, y = dunif(domain, min = 0, max = 1, log = FALSE), col = colors[1], type = "l", lwd = 3,
       ylab = "", xlab = "", main = paste0("gamma, sec. ", codes_S_reg[ss]),
       ylim = c(0, maxY), xlim = c(minX, maxX), 
       cex.main = 2.5, cex.axis=2.5, cex.lab = 2.5)
  polygon(x = c(domain, rev(domain)), y = c(dunif(domain, min = 0, max = 1, log = FALSE), rep(0, length(domain))), col = colors[1], border = NA)
  for (ii in 1:length(temp_hist)) {
    lo <- loess(temp_hist[[ii]]$density ~ temp_hist[[ii]]$mids, span= 1)
    polygon(x = c(temp_hist[[ii]]$mids, rev(temp_hist[[ii]]$mids)), y = c(predict(lo), rep(0, length(temp_hist[[ii]]$mids))), col = colors[ii+1], border = NA)
  }
}

MRIO <- function(theta) {
  W <- construct_W_est(theta, dist_r, dist_t, va_reg[,,N_T])
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
  proj_RS <- rep(rowSums(x_Sio_times_Sio), each = N_R) * stacked_VA_shares #t(X_R) - kronecker(as.matrix(sectors$Y), as.matrix(Y_shares[,ncol(Y_shares)]))#
  x_R <- RAS_blocks(x_R, X_R, VA_Sio, stacked_VA_shares, N_R, N_S_io, x_Sio_times_Sio, proj_CS, proj_RS)
  
  #Create header and save
  x_R_header <- data.frame(sector_id = rep(sectors$id, each = N_R),
                           sector = rep(sectors$name_S_io, each = N_R),
                           region = rep(names_R, N_S_io),
                           stringsAsFactors = FALSE)
  save(x_R, x_R_header, X_R, file = paste0("MRIO_", country, "_", mode, ".RData"))
  
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


