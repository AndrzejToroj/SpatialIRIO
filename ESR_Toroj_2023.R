############## DIRECTORY, LIBRARIES, SETTINGS ################################################

setwd("C:/Andrzej/OneDrive - SGH/regio_leontief/2019_Sonata/ESR_replication_pack/")

if(!require("xlsx")) {install.packages("xlsx"); library(xlsx)}
if(!require("mvtnorm")) {install.packages("mvtnorm"); library(mvtnorm)}
if(!require("Matrix")) {install.packages("Matrix"); library(Matrix)}
if(!require("truncnorm")) {install.packages("truncnorm"); library(truncnorm)}
if(!require("MCMCpack")) {install.packages("MCMCpack"); library(MCMCpack)}
if(!require("matrixcalc")) {install.packages("matrixcalc"); library(matrixcalc)}
if(!require("coda")) {install.packages("coda"); library(coda)}
if(!require("bridgesampling")) {install.packages("bridgesampling"); library(bridgesampling)}
if(!require("sf")) {install.packages("sf"); library(sf)}
if(!require("ggspatial")) {install.packages("ggspatial"); library(ggspatial)}
if(!require("ggplot2")) {install.packages("ggplot2"); library(ggplot2)}
if(!require("reshape2")) {install.packages("reshape2"); library(reshape2)}
if(!require("spdep")) {install.packages("spdep"); library(spdep)}

rm(list = ls())
if(!is.null(dev.list())) dev.off()
cat("\014")

mode <- "base"
h_value <- 0.4 #sensit: 0.3, 0.5
drivetime <- FALSE
if (h_value != 0.4) {mode <- paste0("h", gsub("\\.", "" , as.character(h_value)))}
if (drivetime == TRUE) {mode <- "t"}
country <- "KO"  # "PL", "KO"
estimation <- FALSE #if FALSE, previous estimation results are loaded
useDatarino <- FALSE

#Load functions
source("ESR_Toroj_2023_functions.R")
#...overwrite some of them in "drivetime" mode
if (drivetime == TRUE) {
  source("ESR_Toroj_2023_functions_time.R")
}


############## DATA IMPORTS ################################################

#Poland:
if (country == "PL") {
  #Regional VA data
  ### LOADS: 1 va_reg, 2 codes_S_reg
  va_reg_temp <- array(data=NA, dim=c(73,7,5))
  for (tt in 1:dim(va_reg_temp)[3]) {
    temp <- xlsx::read.xlsx(file = paste0("GVA_PL.xlsx"), 
                            sheetName = paste0("PLN_", 2012 + tt), 
                            colIndex = 1:9,
                            rowIndex = 1:74,
                            stringsAsFactors = FALSE, 
                            header = TRUE, 
                            encoding = "UTF-8")
    colnames(temp) <- c("region", "TOTAL", "A", "B_EexC", "C", "F", "G_J", "K_N", "O_U")
    temp$B_EexC <- temp$B_EexC - temp$C
    temp <- temp[,3:9]
    va_reg_temp[,,tt] <- as.matrix(temp)
  }
  va_reg <- va_reg_temp
  rm(temp, va_reg_temp, tt)
  codes_S_reg <- c("A", "B-EexC", "C", "F", "G-J", "K-N", "O-U")
  
  #I/O data PL & sector aggregation key
  ### LOADS: 3 sectors, 4 x, 5 names_R
  load(paste0("IO_data_", country, ".RData"))
  sectors$id_vareg <- match(sectors$code_S_reg, codes_S_reg)
  
  #Population data
  ### LOADS: 6 popul
  popul <- xlsx::read.xlsx(paste0("Y_indicator_", country, ".xlsx"), sheetName = "Data", header = TRUE, encoding = "UTF-8")
  
  #Distance data: 
  ### LOADS: 7 dist (in 2 versions for Poland if "drivetime" mode on)
  load(paste0("distance_", country, ".RData"))
  if (drivetime == FALSE & country == "PL") {
    dist <- mean_dist
  } else {
    dist_r <- mean_dist
    dist_t <- mean_time
  }
  if (country == "PL") { 
    rm(mean_dist, mean_time) 
  } else {
    drivetime <- FALSE
  }
}

#Korea
if (country == "KO") {
  #x_SR, sectors_SR
  scaling <- 10^3
  x_SR <- xlsx::read.xlsx(file = "KO.xlsx", 
                          sheetName = paste0("transaction"), 
                          colIndex = 1:561,
                          rowIndex = 1:561,
                          stringsAsFactors = FALSE, 
                          header = FALSE, 
                          encoding = "UTF-8") / scaling
  sectors_SR <- xlsx::read.xlsx(file = "KO.xlsx", 
                                sheetName = paste0("header"), 
                                colIndex = 1:8,
                                rowIndex = 1:562,
                                stringsAsFactors = FALSE, 
                                header = TRUE, 
                                encoding = "UTF-8")
  sectors_SR[,c("X", "Y", "VA", "M_abroad")] <- sectors_SR[,c("X", "Y", "VA", "M_abroad")] / scaling
  sum(colSums(x_SR)+sectors_SR$VA <= sectors_SR$X)
  (rowSums(x_SR)+sectors_SR$Y) / sectors_SR$X
  x_SR[x_SR<0] <- 0
  sectors_SR$VA[sectors_SR$VA<0] <- 0
  sectors_SR$X[sectors_SR$X<0] <- 0
  sum(colSums(x_SR)+sectors_SR$VA <= sectors_SR$X)
  sectors_SR$X <- sectors_SR$X + pmax(colSums(x_SR)+sectors_SR$VA-sectors_SR$X,0)
  sum(colSums(x_SR)+sectors_SR$VA <= sectors_SR$X)
  sectors_SR$Y[sectors_SR$Y<0] <- 0
  delta <- (rowSums(x_SR)+sectors_SR$Y) - sectors_SR$X
  sectors_SR$X[delta>0] <- sectors_SR$X[delta>0]+delta[delta>0]
  sectors_SR$Y[delta<0] <- sectors_SR$Y[delta<0]-delta[delta<0]
  sum(colSums(x_SR)+sectors_SR$VA <= sectors_SR$X)
  (rowSums(x_SR)+sectors_SR$Y) / sectors_SR$X
  sum(x_SR<0)
  sum(sectors_SR$VA<0)
  sum(sectors_SR$X<0)
  sum(sectors_SR$Y<0)
  rm(delta)
  
  #sectors
  sectors.aux <- sectors_SR[,-c(1:2,4)]
  codes_S_reg <- c("A", "B-EexC", "C", "F", "G-J", "K-N", "O-U")
  names_S_aux <- unique(sectors_SR$sector)
  sectors <- data.frame(id=unique(sectors_SR$sector_code), name_S_io = names_S_aux,
                        stringsAsFactors = FALSE)
  #sectors$id_vareg <- sectors$id; sectors$name_S_reg <- sectors$name_S_io #no sectorial aggregation for sectorial-regional VA
  sectors.aux <- aggregate(sectors.aux[,2:4], by = list(sector_code = sectors.aux$sector_code), FUN = sum)
  sectors <- merge(x=sectors, y=sectors.aux, by.x = "id", by.y = "sector_code")
  rm(names_S_aux, sectors.aux)
  sectors$prior_sector <- sectors$id_vareg <- c(1, 2, rep(3,14), rep(2,2), 4, rep(5,4), rep(6,4), rep(7,6))
  sectors$code_S_reg <- codes_S_reg[sectors$id_vareg]
  
  #va_reg
  sectors.aggkey <- sectors[,c("id", "code_S_reg")]
  sectors_SR_va <- sectors_SR[, c("region_code", "sector_code", "VA")]
  sectors_SR_va <- merge(x = sectors_SR_va, y = sectors.aggkey, by.x = "sector_code", by.y = "id")
  sectors_SR_va$year <- 1
  va_reg <- acast(data = sectors_SR_va, formula = region_code~code_S_reg~year, fun.aggregate = sum, value.var = "VA") 
  rm(sectors_SR_va)
  
  #x
  x.aux <- aggregate(x_SR, by = list(sector_code = sectors_SR$sector_code), FUN = sum)
  x.aux <- t(as.matrix(x.aux[2:ncol(x.aux)]))
  x.aux <- aggregate(x.aux, by = list(sector_code = sectors_SR$sector_code), FUN = sum)
  x <- t(x.aux[2:ncol(x.aux)])
  rm(x.aux)
  
  #dist
  map <- read_sf("gadm41_KOR_shp", "gadm41_KOR_1")
  map <- st_transform(map, "+proj=longlat")
  names_R <- unique(sectors_SR$region)
  map$NAME_1[12] <- "Jeju-do"
  map$key <- as.numeric(match(map$NAME_1, names_R))
  map <- map[order(map$key),]
  centroids <- st_centroid(map$geometry)
  dist <- st_distance(centroids)
  units(dist) <- "km"
  dist <- matrix(dist, nrow = length(names_R))
  rm(centroids)
  
  source("ESR_Toroj_2023_functions_KO.R")   
  
  #less section T
  #codes_S_reg <- codes_S_reg[sectors$code_S_reg!="T"]
  #va_reg <- array(va_reg[,sectors$code_S_reg!="T",],dim=c(nrow(dist),length(codes_S_reg),1))
  sectors[sectors$id=="S", c("Y", "X", "VA")] <- sectors[sectors$id=="S", c("Y", "X", "VA")]+sectors[sectors$id=="T", c("Y", "X", "VA")]
  sectors <- sectors[!sectors$id=="T",]
  x_SR[sectors_SR$sector_code=="S",] <- x_SR[sectors_SR$sector_code=="S",] + x_SR[sectors_SR$sector_code=="T",]
  x_SR[,sectors_SR$sector_code=="S"] <- x_SR[,sectors_SR$sector_code=="S"] + x_SR[,sectors_SR$sector_code=="T"]
  x_SR <- x_SR[sectors_SR$sector_code!="T",sectors_SR$sector_code!="T"]
  sectors_SR[sectors_SR$sector_code=="S", c("Y", "X", "VA")] <- sectors_SR[sectors_SR$sector_code=="S", c("Y", "X", "VA")]+sectors_SR[sectors_SR$sector_code=="T", c("Y", "X", "VA")]
  sectors_SR <- sectors_SR[!sectors_SR$sector_code=="T",]
  x[nrow(x)-1,] <- x[nrow(x)-1,] + x[nrow(x),] 
  x[,ncol(x)-1] <- x[,ncol(x)-1] + x[,ncol(x)] 
  x <- x[-nrow(x),-ncol(x)]
  
}
  
#Imports complete. Validity checks:
if (drivetime == FALSE) {
  if (dim(dist)[1] != dim(dist)[2]) warning("Distance matrix not square!")
  if (dim(dist)[1] != length(names_R)) warning("Region names vector (names_R) of non-comformable size with distance matrix!")
  if (country == "PL") { if (dim(dist)[1] != nrow(popul)) warning("Regional final demand indicator of non-comformable size with distance matrix!") }
  if (dim(dist)[1] != dim(va_reg)[1]) warning("Regional VA data non conformable with distance matrix!")
  
} else {
  if (dim(dist_r)[1] != dim(dist_r)[2]) warning("Distance matrix not square!")
  if (dim(dist_r)[1] != length(names_R)) warning("Region names vector (names_R) of non-comformable size with distance matrix!")
  if (dim(dist_t)[1] != dim(dist_t)[2]) warning("Distance matrix not square!")
  if (dim(dist_t)[1] != length(names_R)) warning("Region names vector (names_R) of non-comformable size with distance matrix!")
  if (country == "PL") { if (dim(dist_r)[1] != nrow(popul)) warning("Regional final demand indicator of non-comformable size with distance matrix!") }
  if (dim(dist_r)[1] != dim(va_reg)[1]) warning("Regional VA data non conformable with distance matrix!")
}
if (dim(x)[1] != dim(x)[2]) warning("Intermediate transaction matrix not square!")
if (dim(x)[1] != nrow(sectors)) warning("Sectors database (sectors) of non-comformable size with intermediate transaction matrix!")
if (length(unique(sectors$code_S_reg)) != dim(va_reg)[2]) warning("Sector aggregation key inconsistent with regional value added data!")
if (country == "PL") { if (ncol(popul)-1 != dim(va_reg)[3]) warning("Regional final demand indicator inconsistent with VA panel in terms of time!") }
if (sum(is.na(sectors$id_vareg)) > 0) warning("Some national IO sectors unmatched to regional VA data!")


############## DATA PROCESSING ################################################

#More variables extracted from previous imports
N_S_io <- nrow(sectors)
N_S_reg <- dim(va_reg)[2]
N_R <- dim(va_reg)[1]
N_T <- dim(va_reg)[3]

x_Sio_times_Sio <- x
VA_Sio <- sectors$VA

#Aggregate sectors in x; create aggregated VA, Y, X vectors
sec.names <- paste0("s", sectors$id)
sectors.ord <- data.frame(sectors$code_S_reg, 
                          x, 
                          sectors$X, 
                          sectors$Y,
                          sectors$VA,
                          stringsAsFactors = FALSE)
colnames(sectors.ord) <- c("code_S_reg", sec.names, "X", "Y", "VA")
sectors.ord.agg <- aggregate(sectors.ord[,-1], 
                                      by = list(Sector = sectors.ord$code_S_reg), 
                                      FUN = sum)
X <- sectors.ord.agg$X
Y <- sectors.ord.agg$Y
VA <- sectors.ord.agg$VA
sectors.ord.agg.t <- data.frame(sectors.ord$code_S_reg, 
                                         t(sectors.ord.agg[,2:(nrow(x)+1)]), 
                                         stringsAsFactors = FALSE)
sectors.ord.agg.t.agg <- aggregate(sectors.ord.agg.t[,-1], 
                                             by = list(Sector = sectors.ord.agg.t$sectors.ord.code_S_reg), 
                                             FUN = sum)

x <- t(sectors.ord.agg.t.agg[,
                             2:(N_S_reg+1)])
colnames(x) <- codes_S_reg
rownames(x) <- codes_S_reg
rm(sectors.ord, sectors.ord.agg, sectors.ord.agg.t, sectors.ord.agg.t.agg, sec.names)

#Create regional final demand y, remove popul
if (country == "PL") {
  Y_shares <- popul[,-1] / kronecker(t(as.matrix(colSums(popul[,-1]))),as.matrix(rep(1,nrow(popul))))
  y <- list()
  for (tt in 1:N_T) {
    y[[tt]] <- kronecker(Y, as.matrix(Y_shares[, tt]))
  }
  rm(popul, tt)
}
if (country == "KO") {
  sectors.yagg <- sectors_SR[,c("sector_code", "region", "VA")]
  sectors.yagg <- merge(x = sectors.yagg, y = sectors.aggkey, by.x = "sector_code", by.y = "id")
  sectors.yagg <- sectors.yagg[,c("code_S_reg", "region", "VA")]
  sectors.yagg <- aggregate(sectors.yagg$VA, by = list(code_S_reg = sectors.yagg$code_S_reg, region = sectors.yagg$region), FUN = sum)
  sectors.yagg$regord <- match(sectors.yagg$region, names_R)
  sectors.yagg <- sectors.yagg[order(sectors.yagg$code_S_reg, sectors.yagg$regord),]
  
  sectors.tbm <- sectors[,c("id", "Y", "VA")]
  sectors.tbm <- merge(x = sectors.tbm, y = sectors.aggkey, by.x = "id", by.y = "id")
  sectors.tbm <- aggregate(sectors.tbm[,c("VA", "Y")], by = list(code_S_reg = sectors.tbm$code_S_reg), FUN = sum)
  
  sectors.yagg <- merge(x = sectors.yagg, y = sectors.tbm, by = "code_S_reg", all.x = TRUE)
  sectors.yagg$Yr <- sectors.yagg$x/sectors.yagg$VA * sectors.yagg$Y
  
  y <- list()
  y[[1]] <- as.matrix(sectors.yagg$Yr)
  rm(sectors.tbm)
}

#Create VA2: another version of VA vector consistent with regional data (and with additional time dimension)
VA2 <- t(as.matrix(apply(va_reg, c(2, 3), sum)[, 1:N_T]))

#Create beta, beta0 (and A)
beta <- rep(NA, N_S_reg * N_S_reg * N_T)
dim(beta) <- c(N_S_reg, N_S_reg, N_T)
beta0 <- rep(NA, N_S_reg * 1 * N_T)
dim(beta0) <- c(N_S_reg, 1, N_T)
A <- x/kronecker(as.matrix(rep(1,nrow(x))),t(as.matrix(X)))
for (yy in 1:N_T) {
  beta[, , yy] <- kronecker(t(as.matrix(rep(1,N_S_reg))), as.matrix(VA2[yy, ])/as.matrix(X))*A*kronecker(as.matrix(rep(1,N_S_reg)),t(as.matrix(X))/t(as.matrix(VA2[yy, ])))
  beta0[, 1, yy] <- as.matrix(VA2[yy, ])/as.matrix(X)
}

#Complete list of estimation data
data <- list(va_reg = va_reg,
             y = y,
             beta = beta,
             beta0 = beta0,
             dist = dist)
if (drivetime == TRUE) {
  data <- list(va_reg = va_reg,
               y = y,
               beta = beta,
               beta0 = beta0,
               dist_r = dist_r,
               dist_t = dist_t)
}


############## PRIOR DISTRIBUTION ################################################

#Load data for building priors for 7 sectors (EU regional data)
#LOADED: prior_builder (list), shapes
load("prior.RData")
lsh <- log(shapes)
minlsh <- min(lsh)
maxlsh <- max(lsh)

#Hyperparameters: expected values of variances for Wishart
hyper <- rep(NA, N_S_reg+1)
for (ii in 1:N_S_reg) {
  hyper[ii] <- h_value * var(c(va_reg[,ii,]))
}
hyper[N_S_reg+1] <- 0.3


############## ESTIMATION ################################################

if (estimation == TRUE) {

  ###############
  initial_pts_shape <- c(exp(-4), exp(-2), 0.2817238, exp(1), exp(5))
  which_point0 <- 3 #only for start
  set.seed(2)
  S <- 500000
  burnIn <- 0.5*S
  scale <- 0.3 #PL 0.01
  continue <- TRUE
  chain_ID <- 1 #only to continue
  ###############
  
  if (continue == FALSE) {
    #a. start new chain
    selected_shape <- initial_pts_shape[which_point0]
    initial_pts_scale <- expected_scale(selected_shape)
    theta0 <- c(rep(selected_shape,N_S_reg),
                initial_pts_scale,
                rep(hyper[N_S_reg+1], N_S_reg*(N_S_reg-1)/2),
                hyper[1:N_S_reg])
    theta0 <- level2trans(theta0)
    if (drivetime == TRUE) {
      theta0 <- c(theta0[1:(2*N_S_reg)], rep(0.5, N_S_reg), theta0[(2*N_S_reg+1):length(theta0)])
    }
    rm(selected_shape, initial_pts_scale)
  } else {
    # b. continue old chain
    load(paste0("MCMCchains_", mode, country,".RData"))
    eval(parse(text = paste0("theta0 <- MCMCchain", chain_ID, "[nrow(MCMCchain", chain_ID, "),]")))
  }
    
  log_prior(theta0, hyper)
  log_posterior(theta0, hyper, data)
  
  
  ###########################################################
  sample <- run_metropolis_MCMC(theta0, S = S, scale = scale)
  ###########################################################
  
  acceptance <- 1 - mean(duplicated(sample))
  print(paste0("Acceptance probability: ", acceptance))
  if (continue == FALSE) {
    
    ###########################################################
    sample2 <- run_metropolis_MCMC(theta0, S = S, scale = scale)
    ###########################################################
    
    acceptance2 <- 1 - mean(duplicated(sample2))
    print(paste0("Acceptance probability 2: ", acceptance2))
  }
  
  if (continue == FALSE) {
    MCMCchain1 <- sample
    MCMCchain2 <- sample2
  } else {
    eval(parse(text = paste0("MCMCchain", chain_ID, " <- rbind(MCMCchain", chain_ID, ", sample[2:nrow(sample),])")))
  }
  filename <- paste0("MCMCchains_", mode, country, ".RData")
  save(MCMCchain1, MCMCchain2, file = filename)
  rm(sample, sample2)
  
}


############## POST-ESTIMATION CHECKS ################################################

#Load & acceptance probabilities
load(paste0("MCMCchains_", mode, country,".RData"))
theta0 <- MCMCchain1[1,]
S <- nrow(MCMCchain2)
burnIn = 0.5*S
acceptance_1 <- 1 - mean(duplicated(MCMCchain1[- (1 : burnIn), ]))
acceptance_2 <- 1 - mean(duplicated(MCMCchain2[- (1 : burnIn), ]))
MCMC_chainA <- MCMCchain1[(burnIn+1):S,]
MCMC_chainA <- MCMC_chainA[seq(from=1, to = nrow(MCMC_chainA), by = 5), ]
MCMC_chainB <- MCMCchain2[(burnIn+1):S,]
MCMC_chainB <- MCMC_chainB[seq(from=1, to = nrow(MCMC_chainB), by = 5), ]
combined.chains <- mcmc.list(mcmc(MCMC_chainA), 
                             mcmc(MCMC_chainB))

#Priors and posteriors
results_list <- list(MCMCchain1, MCMCchain2)
draw_prior_posterior_all(minlsh, maxlsh, prior_builder, results_list, theta0)
if(drivetime == TRUE) {
  draw_prior_posterior_time(results_list)
}

#Convergence checks
gelman.diag(combined.chains) #PSRF
geweke.diag(combined.chains, frac1 = 0.1, frac2 = 0.5)
autocorr(combined.chains, lags = c(0, 1, 5, 10, 50), relative = TRUE)
autocorr.diag(combined.chains)
autocorr.plot(combined.chains, lag.max = 50)
crosscorr(combined.chains)
crosscorr.plot(combined.chains)

#Posterior means
MCMC_chain <- rbind(MCMC_chainA, MCMC_chainB)
means <- apply(MCMC_chain[,1:(2*N_S_reg)], 2, FUN = mean)
shape_vec <- exp(means[1:N_S_reg])
scale_vec <- exp(means[(N_S_reg+1):(2*N_S_reg)])
params <- data.frame(sector = codes_S_reg, 
                     shape = shape_vec,
                     scale = scale_vec)
write.table(params, file = paste0('output/theta_', country, '.csv'), 
            sep = ';', dec = ',', row.names = FALSE)

#Distance plot
domain <- seq(from=0, to=800, by=5)
intv1 <- list()
intv2 <- list()
intv3 <- list()
str0 <- data.frame(codes_S_reg, avg_dist = rep(NA,N_S_reg))
png(filename = paste0("output/distance_profiles_", mode, country, ".png"), width = 1900, height = 600)
par(mfrow = c(2, 4), mar = c(6, 6, 6, 6), xpd = TRUE)
plot(1, type = "n", xlab = NULL, ylab = NULL, axes = FALSE, ann = FALSE, xpd = TRUE)
legend(x = "bottom", cex = 3, legend = c("post.mean", "50% HPDI", "90% HPDI"), fill = c("black", rgb(0.6, 0.6, 0.6), rgb(0.85, 0.85, 0.85)))
for (ss in 1:length(codes_S_reg)) {
  intv1[[ss]] <- pmax(as.matrix(MCMC_chain[,c(ss, ss+N_S_reg)]), matrix(0.001, nrow = nrow(MCMC_chain), ncol = 2))
  if (nrow(intv1[[ss]]) > 0) {
    intv2[[ss]] <- matrix(NA, nrow = nrow(intv1[[ss]]), ncol = length(domain))
    intv3[[ss]] <- matrix(NA, nrow = 5, ncol = length(domain))
    for (ii in 1:nrow(intv1[[ss]])) {
      intv2[[ss]][ii, ] <- matrix(1, nrow = 1, ncol = length(domain)) - apply(t(as.matrix(domain)), c(1, 2), pgamma, shape = exp(intv1[[ss]][ii, 1]), scale = exp(intv1[[ss]][ii, 2]))
    }
    for (jj in 1:ncol(intv2[[ss]])) {
      intv3[[ss]][1, jj] <- quantile(intv2[[ss]][, jj], probs = 0.05)
      intv3[[ss]][2, jj] <- mean(intv2[[ss]][, jj])
      intv3[[ss]][3, jj] <- quantile(intv2[[ss]][, jj], probs = 0.95)
      intv3[[ss]][4, jj] <- quantile(intv2[[ss]][, jj], probs = 0.25)
      intv3[[ss]][5, jj] <- quantile(intv2[[ss]][, jj], probs = 0.75)
    }
    str0$avg_dist[ss] <- round(t(as.matrix(domain)) %*% as.matrix(intv3[[ss]][2, ]) / sum(intv3[[ss]][2, ]), 0)
    plot(domain, intv3[[ss]][2, ], col = NULL, lwd = 2, type = "l", lty = "solid", ylim = c(0, 1), pch = 20, cex = 3, cex.axis = 3, cex.lab = 3, xlab = 'distance [km]', ylab = "", main = codes_S_reg[ss], cex.main = 3)
    polygon(c(domain, rev(domain)), c(intv3[[ss]][1, ], rev(intv3[[ss]][3, ])), col = rgb(0.85, 0.85, 0.85), border = NA)
    polygon(c(domain, rev(domain)), c(intv3[[ss]][4, ], rev(intv3[[ss]][5, ])), col = rgb(0.6, 0.6, 0.6), border = NA)
    lines(x = domain, y = intv3[[ss]][2, ], lwd = 4, col = "black", lty = "solid")
    text(x=700,y=0.9, paste0("avg = ", as.character(str0$avg_dist[ss])), cex = 3)
  }
}
dev.off()

#Balanced MRIO
theta <- c(shape_vec, scale_vec)
if (drivetime == TRUE) {theta <- c(theta, apply(MCMC_chain[,(2*N_S_reg+1):(3*N_S_reg)], 2, FUN = mean))}
A_R_ind <- MRIO(theta)
if (country == "KO") {
  load(paste0("MRIO_", country, "_base.RData")) #X_R, X_R_header, x_R
  sectors_SR$regord <- match(sectors_SR$region, names_R)
  
  heatmap(x_Sio_times_Sio, Rowv=NA, Colv=NA, revC=TRUE)
  indexing1 <- order(sectors_SR$sector_code, sectors_SR$regord)
  x_SR_i1 <- as.matrix(x_SR)[indexing1,indexing1] #R slow -> fast, S fast -> slow
  x_R_i1 <- x_R
  heatmap(x_SR_i1, Rowv=NA, Colv=NA, revC=TRUE)
  heatmap(x_R_i1, Rowv=NA, Colv=NA, revC=TRUE)
  
  x_SR_i2 <- as.matrix(x_SR)
  aux <- sectors_SR[order(sectors_SR$sector_code, sectors_SR$regord), c("sector_code", "regord")]
  indexing2 <- order(aux$regord, aux$sector_code)
  x_R_i2 <- x_R[indexing2,indexing2]
  heatmap(x_SR_i2, Rowv=NA, Colv=NA, revC=TRUE)
  heatmap(x_R_i2, Rowv=NA, Colv=NA, revC=TRUE)
  
  errors_i2 <- x_SR_i2 - x_R_i2
  d.block <- kronecker(diag(N_R), matrix(1, nrow = N_S_io, ncol = N_S_io))
  nd.block <- 1-d.block
  
  cname <- "comparison_spatial"
  eval(parse(text = paste0(cname, " <- list()")))
  
  eval(parse(text = paste0(cname, "$MSE_overall <- sum(errors_i2^2)/((N_R*N_S_io)^2)")))
  eval(parse(text = paste0(cname, "$RMSE_overall <- sqrt(",cname,"$MSE_overall)")))
  eval(parse(text = paste0(cname, "$MSE_diag <- sum(errors_i2^2 * d.block)/(N_R*(N_S_io^2))")))
  eval(parse(text = paste0(cname, "$RMSE_diag <- sqrt(",cname,"$MSE_diag)")))
  eval(parse(text = paste0(cname, "$MSE_nondiag <- sum(errors_i2^2 * nd.block)/((N_R*N_S_io)^2-N_R*(N_S_io^2))")))
  eval(parse(text = paste0(cname, "$RMSE_nondiag <- sqrt(",cname,"$MSE_nondiag)")))
  
  fname <- paste0("output/RMSE_", cname, ".RData")
  eval(parse(text = paste0("save(", cname, ", file = '", fname, "')")))
  
  save.image("javaerror.RData")
  #close RStudio, then re-run the preamble 1-20 and continue below
  
  #here estimate gravity model for JP
  load("javaerror.RData")
  tblsrcv <- "Icio" #"Jp" # "Icio"
  if (tblsrcv == "Jp") {
    x_SR_JP <- xlsx::read.xlsx(file = "JP.xlsx", 
                               sheetName = "table", 
                               colIndex = 1:658,
                               rowIndex = 1:576,
                               stringsAsFactors = FALSE, 
                               header = FALSE, 
                               encoding = "UTF-8")
  }
  if (tblsrcv == "Icio") {x_SR_JP <- read.csv("ICIO.csv", header = TRUE, sep = ",", dec = "."); x_SR_JP <- x_SR_JP[,-1]}
  if (tblsrcv == "Jp") {
    x_SR_JP_row <- xlsx::read.xlsx(file = "JP.xlsx", 
                                   sheetName = "rowHeader", 
                                   colIndex = 1:4,
                                   rowIndex = 1:577,
                                   stringsAsFactors = FALSE, 
                                   header = TRUE, 
                                   encoding = "UTF-8")
  }
  if (tblsrcv == "Icio") {
    x_SR_JP_row <- xlsx::read.xlsx(file = "ICIO.xlsx", 
                                   sheetName = "rowHeader", 
                                   colIndex = 1:4,
                                   rowIndex = 1:3469,
                                   stringsAsFactors = FALSE, 
                                   header = TRUE, 
                                   encoding = "UTF-8")
  }
  save.image("javaerror.RData")
  #close RStudio, then re-run the preamble 1-20 and continue below
  
  load("javaerror.RData")
  if (tblsrcv == "Jp") {
    x_SR_JP_col <- xlsx::read.xlsx(file = "JP.xlsx", 
                                   sheetName = "colHeader", 
                                   colIndex = 1:4,
                                   rowIndex = 1:659,
                                   stringsAsFactors = FALSE, 
                                   header = TRUE, 
                                   encoding = "UTF-8")
  }
  if (tblsrcv == "Icio") {
    x_SR_JP_col <- xlsx::read.xlsx(file = "ICIO.xlsx", 
                                   sheetName = "colHeader", 
                                   colIndex = 1:4,
                                   rowIndex = 1:3544,
                                   stringsAsFactors = FALSE, 
                                   header = TRUE, 
                                   encoding = "UTF-8")
  }
  
  if (tblsrcv == "Jp") {
    row.filter <- x_SR_JP_row$Region.Code..Row. %in% 1:9 &
                x_SR_JP_row$Sector.Code..Row. %in% seq(10,530, by=10)
    col.filter <- x_SR_JP_col$Region.Code..Column. %in% 1:9 &
                  x_SR_JP_col$Sector.Code..Column. %in% seq(10,530, by=10)
    x_SR_JP <- x_SR_JP[row.filter,]
    X_R_JP <- x_SR_JP[,ncol(x_SR_JP)]
    x_SR_JP <- x_SR_JP[,col.filter] 
    x_SR_JP_header <- x_SR_JP_row[row.filter,]
    colnames(x_SR_JP_header) <- c("region_code", "region_name", "sector_code", "sector_name")
    rm(x_SR_JP_row, x_SR_JP_col, col.filter, row.filter)
  }
  if (tblsrcv == "Icio") {
    X_R_JP <- x_SR_JP[1:3465,3543]
    x_SR_JP <- x_SR_JP[1:3465,1:3465]
    x_SR_JP_header <- x_SR_JP_row[1:3465,]
    colnames(x_SR_JP_header) <- c("region_code", "region_name", "sector_code", "sector_name")
    rm(x_SR_JP_row, x_SR_JP_col)
  }

  X_JP <- data.frame(sector = x_SR_JP_header$sector_code,
                     region = x_SR_JP_header$region_code,
                     X = X_R_JP,
                     stringsAsFactors = FALSE)
  flows_JP <- cbind(x_SR_JP_header[,c("region_code", "sector_code")], x_SR_JP)
  colnames(flows_JP) <- c("source_region", "source_sector", paste("recip", flows_JP$region_code, flows_JP$sector_code, sep = "_"))
  flows_JP <- melt(data = flows_JP, 
                   id.vars = c("source_region", "source_sector"),
                   variable.name = "destination",
                   value.name = "flows")
  flows_JP$destination <- as.character(flows_JP$destination)
  aux <- colsplit(flows_JP$destination, "_", c("recip", "dest_region", "dest_sector"))
  flows_JP <- cbind(flows_JP[, c("source_region", "source_sector")],
                    aux[, c("dest_region", "dest_sector")],
                    flows_JP$flows)
  rm(aux)
  colnames(flows_JP)[5] <- "flows"
  flows_JP$source <- paste0(flows_JP$source_region, "_", flows_JP$source_sector)
  flows_JP$dest <- paste0(flows_JP$dest_region, "_", flows_JP$dest_sector)
  flows_JP$source_dest <- paste0(flows_JP$source_region, "_", flows_JP$dest_region)
  X_JP$region_sector <- paste0(X_JP$region, "_", X_JP$sector)
  flows_JP <- merge(x=flows_JP, y=X_JP[,3:4], by.x="source", by.y="region_sector",
                    sort = FALSE, all = TRUE)
  colnames(flows_JP)[9] <- "source_X"
  flows_JP <- merge(x=flows_JP, y=X_JP[,3:4], by.x="dest", by.y="region_sector",
                    sort = FALSE, all = TRUE)
  colnames(flows_JP)[10] <- "dest_X"
  if (tblsrcv == "Jp") {
    flows_JP <- flows_JP[flows_JP$source_region != flows_JP$dest_region, ]
  }
  if (tblsrcv == "Icio") {
    flows_JP <- flows_JP[flows_JP$source_region != "ROW" & flows_JP$dest_region != "ROW", ]
  }
  
  #dist
  if (tblsrcv == "Jp") {
    map_JP <- read_sf("jpn_adm_2019_shp", "jpn_admbnda_adm1_2019")
    map_JP <- st_transform(map_JP, "+proj=longlat")
    names_R_JP <- unique(x_SR_JP_header$region_name)
    names_R_JP_map <- c("Hokkaido", "Miyagi", "Tokyo", "Aichi", "Osaka", "Hiroshima", "Ehime", "Fukuoka", "Okinawa")
    names_R_JP <- data.frame(names_IO = names_R_JP, names_map = names_R_JP_map, region_code = 1:9)
    map_JP$ADM1_EN <- substr(map_JP$ADM1_EN, 2, nchar(map_JP$ADM1_EN))
    map_JP <- merge(x=map_JP, y=names_R_JP, by.x="ADM1_EN", by.y="names_map")
    map_JP <- map_JP[order(map_JP$region_code),]
    centroids <- st_centroid(map_JP$geometry)
    dist_JP <- st_distance(centroids)
    units(dist_JP) <- "km"
    dist_JP <- matrix(dist_JP, nrow = nrow(names_R_JP))
    rm(map_JP, centroids)
    adj_JP <- matrix(0, nrow = nrow(names_R_JP), ncol = nrow(names_R_JP))
    adj_JP[2,3] <- 1; adj_JP[3,2] <- 1
    adj_JP[2,4] <- 1; adj_JP[4,2] <- 1
    adj_JP[3,4] <- 1; adj_JP[4,3] <- 1
    adj_JP[4,5] <- 1; adj_JP[5,4] <- 1
    adj_JP[5,6] <- 1; adj_JP[6,5] <- 1
    dist_JP <- cbind(names_R_JP[,c(1,3)], dist_JP)
    colnames(dist_JP)[3:11] <- paste0("r_", names_R_JP$region_code)
    adj_JP <- cbind(names_R_JP[,c(1,3)], adj_JP)
    colnames(adj_JP)[3:11] <- paste0("r_", names_R_JP$region_code)
  }
  if (tblsrcv == "Icio") {
    map_JP <- read_sf("world", "world-administrative-boundaries")
    map_JP <- st_transform(map_JP, "+proj=longlat")
    unique(flows_JP$source_region) %in% unique(map_JP$iso3)
    map_JP <- map_JP[order(map_JP$iso3),]
    map_JP <- map_JP[map_JP$iso3 %in% unique(flows_JP$source_region), ]
    map_JP <- map_JP[-61,]
    (names_R_JP <- unique(x_SR_JP_header$region_name)[1:76])
    map_JP$iso3
    centroids <- st_centroid(map_JP$geometry)
    dist_JP <- st_distance(centroids)
    units(dist_JP) <- "km"
    dist_JP <- as.matrix(data.frame(dist_JP))#, nrow = nrow(names_R_JP))
    rm(centroids)
    
    adj_JP <- nb2mat(poly2nb(map_JP, queen = T), style = "B", zero.policy = TRUE)
    rm(map_JP)

    dist_JP <- cbind(as.matrix(names_R_JP), as.data.frame(dist_JP))
    colnames(dist_JP) <- c("region_code", paste0("r_", names_R_JP))
    adj_JP <- cbind(as.matrix(names_R_JP), as.data.frame(adj_JP))
    colnames(adj_JP) <- c("region_code", paste0("r_", names_R_JP))
  }
  
  if (tblsrcv == "Jp") {
    dist_JP <- melt(data = dist_JP, 
                     id.vars = c("names_IO", "region_code"),
                     variable.name = "region_code_2",
                     value.name = "distance")
    aux <- colsplit(as.character(dist_JP$region_code_2), "_", c("r", "region_code_2"))
    dist_JP$region_code_2 <- aux$region_code_2
    dist_JP$source_dest <- paste0(dist_JP$region_code, "_", dist_JP$region_code_2)
    adj_JP <- melt(data = adj_JP, 
                    id.vars = c("names_IO", "region_code"),
                    variable.name = "region_code_2",
                    value.name = "adjacency")
    aux <- colsplit(as.character(adj_JP$region_code_2), "_", c("r", "region_code_2"))
    adj_JP$region_code_2 <- aux$region_code_2
    adj_JP$source_dest <- paste0(adj_JP$region_code, "_", adj_JP$region_code_2)
    rm(aux)
    adj_JP <- adj_JP[,4:5]
    dist_JP <- dist_JP[,4:5]
  }
  
  if (tblsrcv == "Icio") {
    dist_JP <- melt(data = dist_JP, 
                    id.vars = c("region_code"),
                    variable.name = "region_code_2",
                    value.name = "distance")
    aux <- colsplit(as.character(dist_JP$region_code_2), "_", c("r", "region_code_2"))
    dist_JP$region_code_2 <- aux$region_code_2
    dist_JP$source_dest <- paste0(dist_JP$region_code, "_", dist_JP$region_code_2)
    adj_JP <- melt(data = adj_JP, 
                   id.vars = c("region_code"),
                   variable.name = "region_code_2",
                   value.name = "adjacency")
    aux <- colsplit(as.character(adj_JP$region_code_2), "_", c("r", "region_code_2"))
    adj_JP$region_code_2 <- aux$region_code_2
    adj_JP$source_dest <- paste0(adj_JP$region_code, "_", adj_JP$region_code_2)
    rm(aux)
    adj_JP <- adj_JP[,3:4]
    dist_JP <- dist_JP[,3:4]
  }
  

  dist_JP <- merge(x=dist_JP, y=adj_JP, by = "source_dest")
  rm(adj_JP)
  flows_JP <- merge(x = flows_JP, y = dist_JP, by = "source_dest", all.x = TRUE)
  flows_JP <- flows_JP[flows_JP$source_X != 0 & flows_JP$dest_X != 0,]
  
  if (tblsrcv == "Jp") {
    mapping <- xlsx::read.xlsx(file = "JP.xlsx", 
                                   sheetName = "mapping", 
                                   colIndex = 1:3,
                                   rowIndex = 1:54,
                                   stringsAsFactors = FALSE, 
                                   header = TRUE, 
                                   encoding = "UTF-8")[,c(1,3)]
  }
  if (tblsrcv == "Icio") {
    mapping <- xlsx::read.xlsx(file = "ICIO.xlsx", 
                               sheetName = "mapping", 
                               colIndex = 1:3,
                               rowIndex = 1:46,
                               stringsAsFactors = FALSE, 
                               header = TRUE, 
                               encoding = "UTF-8")[,c(1,3)]
  }
  flows_JP <- merge(x=flows_JP, y=mapping, by.x="source_sector", by.y = "code_JP")
  
  flows_JP$dist_A <- as.integer(flows_JP$code_S_reg == "A")*log(flows_JP$distance)
  flows_JP$dist_BEexC <- as.integer(flows_JP$code_S_reg == "B-EexC")*log(flows_JP$distance)
  flows_JP$dist_C <- as.integer(flows_JP$code_S_reg == "C")*log(flows_JP$distance)
  flows_JP$dist_F <- as.integer(flows_JP$code_S_reg == "F")*log(flows_JP$distance)
  flows_JP$dist_GJ <- as.integer(flows_JP$code_S_reg == "G-J")*log(flows_JP$distance)
  flows_JP$dist_KN <- as.integer(flows_JP$code_S_reg == "K-N")*log(flows_JP$distance)
  flows_JP$dist_OU <- as.integer(flows_JP$code_S_reg == "O-U")*log(flows_JP$distance)
  flows_JP$adj_A <- as.integer(flows_JP$code_S_reg == "A")*flows_JP$adjacency
  flows_JP$adj_BEexC <- as.integer(flows_JP$code_S_reg == "B-EexC")*flows_JP$adjacency
  flows_JP$adj_C <- as.integer(flows_JP$code_S_reg == "C")*flows_JP$adjacency
  flows_JP$adj_F <- as.integer(flows_JP$code_S_reg == "F")*flows_JP$adjacency
  flows_JP$adj_GJ <- as.integer(flows_JP$code_S_reg == "G-J")*flows_JP$adjacency
  flows_JP$adj_KN <- as.integer(flows_JP$code_S_reg == "K-N")*flows_JP$adjacency
  flows_JP$adj_OU <- as.integer(flows_JP$code_S_reg == "O-U")*flows_JP$adjacency
  
  if (tblsrcv == "Icio") {
    eur.smpl <- TRUE
    if (eur.smpl == TRUE) {
      eur <- c("AUT", "BEL", "BGR", "CYP", "CZE", "DEU", "DNK", "ESP", "EST", "FIN", "FRA", "GRC", "HRV",
               "HUN", "IRL", "ITA", "LTU", "LUX", "LVA",
               "MLT", "NLD", "POL", "PRT", "ROU", "SVK", "SVN", "SWE")
      flows_JP <- flows_JP[flows_JP$dest_region %in% eur & flows_JP$source_region %in% eur, ]
      tblsrcv <- "IcioEur"
    }
  }
  
  gravity_model_v1 <- glm(flows ~ log(source_X) + log(dest_X) + log(distance) + adjacency, data = flows_JP, family = poisson)
  gravity_model_v2 <- glm(flows ~ log(source_X) + log(dest_X) + 
                            dist_A + dist_BEexC + dist_C + dist_F + dist_GJ + dist_KN + dist_OU +
                            adj_A + adj_BEexC + adj_C + adj_F + adj_GJ + adj_KN + adj_OU, 
                          data = flows_JP, family = poisson)

  options(scipen = 100000)
  summary(gravity_model_v1)
  summary(gravity_model_v2)
  grav.het <- TRUE ###!!!###
  if (grav.het == TRUE) {gravv <- "Ghet"} else {gravv <- "Ghom"}
  
  ### start pasted to adjust for Korea
  delta_v1 <- rep(0.385, N_R)
      dom.imp <- 1 - kronecker(diag(N_R), matrix(1, nrow = N_S_io, ncol = N_S_io))
      dom.imp <- colSums(dom.imp * x_SR)
      sectors_SR$M_home <- dom.imp
      rm(dom.imp)
      regions <- aggregate(sectors_SR[, c("X", "M_abroad", "M_home")], by = list(region = sectors_SR$region_code), FUN = sum)
      regions$R <- regions$X / sum(regions$X) * 100
      regions$P <- regions$M_home / regions$X
      regions$P <- regions$P / mean(regions$P)
      regions$F <- regions$M_abroad / regions$X
      regions$F <- regions$F / mean(regions$F)
      regions$D <- 0; regions$D[4] <- 1
      regions$delta <- exp(-1.226 + 0.168*log(regions$R) + 0.325*log(regions$P) + 0.317*log(regions$F) + 0.577*regions$D)
      mean(regions$delta)
  delta_v2 <- regions$delta
  delta <- delta_v2 ###!!!###
  if (sum(delta == delta_v1) == length(delta)) {deltav <- "Dhom"} else {deltav <- "Dhet"}
  
  VA_SR <- as.matrix(sectors_SR$VA)
  VA_SR <- matrix(VA_SR, nrow = N_R, byrow = TRUE)
  X <- as.matrix(sectors$X)
  VA <- as.matrix(sectors$VA)
  A <- x_Sio_times_Sio / kronecker(as.matrix(rep(1,length(X))), t(X))
  A[is.nan(A)] <- 0
  
  SLQ_denominators <- VA / sum(VA)
  SLQ_numerators <- VA_SR / kronecker(as.matrix(rowSums(VA_SR)), t(as.matrix(rep(1, N_S_io))))
  SLQ <- SLQ_numerators / kronecker(as.matrix(rep(1, N_R)), t(as.matrix(SLQ_denominators)))
  SLQ[is.nan(SLQ)] <- 0
  A_SLQ_intrareg <- list()
  for (ii in 1:N_R) {
    A_SLQ_intrareg[[ii]] <- ifelse(kronecker(as.matrix(SLQ[ii, ]), t(as.matrix(rep(1,N_S_io)))) >= 1, A, A*kronecker(as.matrix(SLQ[ii, ]), t(as.matrix(rep(1,N_S_io)))))
  }
  
  CILQ_numerator <- list()
  CILQ_denominator <- list()
  CILQ <- list()
  A_CILQ_intrareg <- list()
  for (ii in 1:N_R) {
    CILQ_numerator[[ii]] <- kronecker(as.matrix(SLQ[ii, ]), t(as.matrix(rep(1, N_S_io))))
    CILQ_denominator[[ii]] <- kronecker(t(as.matrix(SLQ[ii, ])), as.matrix(rep(1, N_S_io)))
    CILQ[[ii]] <- CILQ_numerator[[ii]] / CILQ_denominator[[ii]]
    CILQ[[ii]][is.nan(CILQ[[ii]])] <- 0
    diag(CILQ[[ii]]) <- SLQ[ii, ]
    A_CILQ_intrareg[[ii]] <- ifelse(CILQ[[ii]] >= 1, A, A*CILQ[[ii]])
  }
  
  FLQ <- list()
  A_FLQ_intrareg <- list()
  lambda <- (log2(1+rowSums(VA_SR)/sum(VA)))^delta
  for (ii in 1:N_R) {
    FLQ[[ii]] <- lambda[ii] * CILQ[[ii]]
    FLQ[[ii]][is.nan(FLQ[[ii]])] <- 0
    diag(FLQ[[ii]]) <- SLQ[ii, ]
    A_FLQ_intrareg[[ii]] <- ifelse(FLQ[[ii]] >= 1, A, A*FLQ[[ii]])
  }
  
  A_reg <- A_FLQ_intrareg
  
  ord <- order(sectors_SR$sector_code, sectors_SR$region_code)
  stacked_VA_shares <- sectors_SR$VA[ord] / rep(VA, each = N_R)
  X_R_hat <- rep(as.vector(X), each = N_R) * stacked_VA_shares
  X_R_hat <- matrix(X_R_hat, byrow = FALSE, nrow = N_R)
  x_reg <- list()
  x_reg_sum <- matrix(0, nrow = N_S_io, ncol = N_S_io)
  for (ii in 1:N_R) {
    x_reg[[ii]] <- A_reg[[ii]] * (as.matrix(rep(1,N_S_io)) %*% t(as.matrix(X_R_hat[ii, ])))
    x_reg_sum <- x_reg_sum + x_reg[[ii]]
  }
  LQ_residuals <- x_Sio_times_Sio - x_reg_sum
  LQ_residuals[LQ_residuals < 0] <- 0
  
  h <- list()
  h_sum <- matrix(0, nrow = N_S_io, ncol = N_S_io)
  border <- nb2mat(poly2nb(map, queen = T), style = "B", zero.policy = TRUE)
  if(grav.het == FALSE) {
    beta1 <- gravity_model_v1$coefficients["log(distance)"]
    beta2 <- gravity_model_v1$coefficients["log(source_X)"]
    beta3 <- gravity_model_v1$coefficients["log(dest_X)"]
    beta4 <- gravity_model_v1$coefficients["adjacency"]
  } else {
    beta1.aux <- as.numeric(gravity_model_v2$coefficients[c("dist_A", "dist_BEexC", "dist_C", "dist_F", "dist_GJ", "dist_KN", "dist_OU")])
    beta4.aux <- as.numeric(gravity_model_v2$coefficients[c("adj_A", "adj_BEexC", "adj_C", "adj_F", "adj_GJ", "adj_KN", "adj_OU")])
    beta.aux <- data.frame(beta1 = beta1.aux, beta4 = beta4.aux, sec = codes_S_reg)
    beta2 <- gravity_model_v2$coefficients["log(source_X)"]
    beta3 <- gravity_model_v2$coefficients["log(dest_X)"]
    sectors.aux <- merge(x = sectors[,c("id", "code_S_reg")], y = beta.aux, by.x = "code_S_reg", by.y = "sec")
    beta1 <- kronecker(matrix(1, nrow = 1, ncol = N_S_io), as.matrix(sectors.aux$beta1))
    beta4 <- kronecker(matrix(1, nrow = 1, ncol = N_S_io), as.matrix(sectors.aux$beta4))
  }
  for (ii in 1:N_R) {
    h[[ii]] <- list()
    for (jj in 1:N_R) {


        
        h[[ii]][[jj]] <- dist[ii, jj]^beta1 * ((as.matrix(X_R_hat[ii, ])^beta2) %*% t((as.matrix(X_R_hat[jj, ])^beta3))) * exp(border[ii,jj]*beta4)
        h[[ii]][[jj]][is.nan(h[[ii]][[jj]]) | h[[ii]][[jj]]<0] <- 0

      
      if (ii == jj) { h[[ii]][[jj]] <- matrix(0, nrow = N_S_io, ncol = N_S_io) }
      h_sum <- h_sum + h[[ii]][[jj]]

    }
  }
  
  g <- list()
  z <- list()
  for (ii in 1:N_R) {
    g[[ii]] <- list()
    z[[ii]] <- list()
    for (jj in 1:N_R) {
      g[[ii]][[jj]] <- h[[ii]][[jj]] / h_sum
      g[[ii]][[jj]][is.nan(g[[ii]][[jj]])] <- 0
      z[[ii]][[jj]] <- g[[ii]][[jj]] * LQ_residuals
    }
  }
  
  x_R_iriolq <- matrix(NA, nrow = N_S_io*N_R, ncol = N_S_io*N_R)
  for (ii in 1:N_R) {
    for (jj in 1:N_R) {
      if (ii == jj) {
        x_R_iriolq[((ii-1)*N_S_io+1):(ii*N_S_io), ((jj-1)*N_S_io+1):(jj*N_S_io)] <- x_reg[[ii]]
      } else {
        x_R_iriolq[((ii-1)*N_S_io+1):(ii*N_S_io), ((jj-1)*N_S_io+1):(jj*N_S_io)] <- z[[ii]][[jj]]
      }
    }
  }
  
  #save.image("optimization.RData")
  
            load("optimization.RData")
            #transform to final version!
            cond1_LHS <- matrix(colSums(x_R_iriolq), nrow = N_R, byrow=TRUE)
            sum(x_Sio_times_Sio)
            sum(cond1_LHS)
            sum(cond1_LHS <= X_R_hat)
            sum(cond1_LHS > X_R_hat)
            test <- cond1_LHS - X_R_hat
            test[test>0]
            
            library(nloptr)
            x0 <- as.vector(x_R_iriolq) #reg slow, sec fast
            lb <- rep(0, length(x0))
            x_init <- as.vector(x_R_iriolq)
            x_init[x_init==0] <- 0.0001
            f_opt <- function(x) {
              f <- sum(((x-x_init)^2)/x_init)
              g <- (2*(x-x_init)) / (x_init)
              return(list("objective" = f, "gradient" = g))
            }
            f_opt(x0)
            f_constr_ineq <- function(x) {
              lhs1 <- colSums(matrix(x, nrow = N_R*N_S_io, ncol = N_R*N_S_io, byrow = FALSE))
              rhs1 <- as.vector(t(X_R_hat))
              f1 <- lhs1-rhs1 #!!!!!!!!!!!!
              g1 <- kronecker(diag(N_R*N_S_io), matrix(1, nrow=1,ncol=N_R*N_S_io))
              return(list("constraints" = f1, "jacobian" = g1))
            }
            sum(f_constr_ineq(x0)$constraints>0)
            f_constr_eq <- function(x) {
              temp <- matrix(x, nrow = N_R*N_S_io, byrow = FALSE)
              temp <- aggregate(temp, by = list(sectors_SR$sector_code), FUN = sum)[,-1]
              temp <- t(temp)
              temp <- aggregate(temp, by = list(sectors_SR$sector_code), FUN = sum)[,-1]
              temp <- t(temp)
              lhs2 <- as.vector(temp)
              rhs2 <- as.vector(x_Sio_times_Sio)
              f2 <- lhs2-rhs2
              bb_1 <- kronecker(matrix(1, nrow=1, ncol=N_R), diag(N_S_io))
              bb_2 <- kronecker(diag(N_S_io), bb_1)
              g2 <- kronecker(matrix(1, nrow=1, ncol=N_R), bb_2)
              return(list("constraints" = f2, "jacobian" = g2))
            }
            sum(f_constr_eq(x0)$constraints < -0.00001 | f_constr_eq(x0)$constraints > 0.00001)
            test <- f_constr_eq(x0)$jacobian
            local_opts <- list("print_level" = 2)
            res <- nloptr( x0=x0,
                            eval_f=f_opt,
                            lb = lb,
                            eval_g_ineq = f_constr_ineq,
                            eval_g_eq = f_constr_eq,
                            opts = list("algorithm"="NLOPT_LD_MMA", "check_derivatives"=TRUE, "print_level"=2, "local_opts" = local_opts))
            x_R_iriolq_final <- matrix(res$solution, nrow = N_R*N_S_io, ncol = N_R*N_S_io, byrow = FALSE)
            save(x_R_iriolq_final, file = "iriolq_final_version.RData")
  
  ###
            
  ###
            
  load("iriolq_final_version.RData")

        sum(abs(x_R_iriolq_final - x_R_iriolq)>1)
        cor(as.vector(x_R_iriolq_final), as.vector(x_R_iriolq))
        cond1_LHS_v2 <- matrix(colSums(x_R_iriolq_final), nrow = N_R, byrow=TRUE)
        sum(x_Sio_times_Sio)
        sum(cond1_LHS_v2)
        sum(cond1_LHS_v2 <= X_R_hat)
        sum(cond1_LHS_v2 > X_R_hat)
        test2 <- cond1_LHS_v2 - X_R_hat
        test2[test2>0]

  x_R_iriolq_i2 <- x_R_iriolq_final   #region as slow dimension
  x_R_iriolq_i1 <- x_R_iriolq_final[indexing1,indexing1] #region as fast dimension
  heatmap(x_SR_i1, Rowv=NA, Colv=NA, revC=TRUE)
  heatmap(x_R_iriolq_i1, Rowv=NA, Colv=NA, revC=TRUE)
  heatmap(x_SR_i2, Rowv=NA, Colv=NA, revC=TRUE)
  heatmap(x_R_iriolq_i2, Rowv=NA, Colv=NA, revC=TRUE)
  
  errors_i2 <- x_SR_i2 - x_R_iriolq_i2
  cname <- paste0("comparison_", tblsrcv, deltav, gravv)
  eval(parse(text = paste0(cname, " <- list()")))
  
  eval(parse(text = paste0(cname, "$MSE_overall <- sum(errors_i2^2)/((N_R*N_S_io)^2)")))
  eval(parse(text = paste0(cname, "$RMSE_overall <- sqrt(",cname,"$MSE_overall)")))
  eval(parse(text = paste0(cname, "$MSE_diag <- sum(errors_i2^2 * d.block)/(N_R*(N_S_io^2))")))
  eval(parse(text = paste0(cname, "$RMSE_diag <- sqrt(",cname,"$MSE_diag)")))
  eval(parse(text = paste0(cname, "$MSE_nondiag <- sum(errors_i2^2 * nd.block)/((N_R*N_S_io)^2-N_R*(N_S_io^2))")))
  eval(parse(text = paste0(cname, "$RMSE_nondiag <- sqrt(",cname,"$MSE_nondiag)")))
  fname <- paste0("output/RMSE_", cname, ".RData")
  eval(parse(text = paste0("save(", cname, ", file = '", fname, "')")))
  
  ### when all options are explored, run...
  # cm <- ""
  # for (ee in ls(pattern = "comparison*")) {
  #   cm <- paste0(cm, ", ", ee)
  # }
  # cm <- substr(cm, 3, nchar(cm))
  # eval(parse(text = paste0("save(",cm,",file='comparison.RData')")))
  load("comparison.RData")
  results_table <- data.frame(scenario = rep(NA, 9),
                              RMSE_total = rep(NA, 9),
                              RMSEt_pc_spatial = rep(NA, 9),
                              RMSE_blockdiagonal = rep(NA, 9),
                              RMSEbd_pc_spatial = rep(NA, 9))
  table_seq <- c("spatial", "JpDhomGhom", "JpDhomGhet",
                            "IcioEurDhomGhom", "IcioEurDhomGhet",
                            "IcioDhomGhom", "IcioDhomGhet",
                            "JpDhetGhom", "JpDhetGhet")
  results_table$scenario <- table_seq
  for (rr in 1:9) {
    eval(parse(text = paste0("results_table$RMSE_total[rr] <- comparison_", table_seq[rr], "$RMSE_overall")))
    eval(parse(text = paste0("results_table$RMSE_blockdiagonal[rr] <- comparison_", table_seq[rr], "$RMSE_diag")))
    eval(parse(text = paste0("results_table$RMSEt_pc_spatial[rr] <- comparison_", table_seq[rr], "$RMSE_overall / comparison_spatial$RMSE_overall")))
    eval(parse(text = paste0("results_table$RMSEbd_pc_spatial[rr] <- comparison_", table_seq[rr], "$RMSE_diag / comparison_spatial$RMSE_diag")))
  }
  save(results_table, file = "comparison_results_table.RData")
}

if(country == "PL") {
  A_R <- A_R_ind[1:(N_S_io*N_R), 1:(N_S_io*N_R)]
  
  #Simulation baseline
  Sim_region <- 27 #PL514
  Sim_sector <- 45
  Sim_sales <- 200
  Sim_cost <- 100
  
  Sim_Acol <- A_R[,(Sim_sector-1)*N_R+Sim_region]
  Sim_Acol <- Sim_Acol / sum(Sim_Acol) * Sim_cost / Sim_sales
  A_R_sim <- rbind(cbind(A_R, as.matrix(Sim_Acol)),
                   matrix(0, nrow = 1, ncol = ncol(A_R)+1))
  dY <- matrix(c(rep(0, N_R*N_S_io), Sim_sales), ncol=1)
  dX <- solve(diag(nrow(A_R)+1)-A_R_sim, dY)
  dX_tab <- matrix(dX[1:(N_S_io*N_R)], nrow = N_R, ncol = N_S_io, byrow = FALSE)
  
  Sim_Acol_ind <- A_R_ind[,(Sim_sector-1)*N_R+Sim_region]
  A_R_ind_sim <- rbind(cbind(A_R_ind, as.matrix(Sim_Acol_ind)),
                       matrix(0, nrow = 1, ncol = ncol(A_R_ind)+1))
  dY_ind <- matrix(c(rep(0, N_R*(N_S_io+1)), Sim_sales), ncol=1)
  dX_ind <- solve(diag(nrow(A_R_ind)+1)-A_R_ind_sim, dY_ind)
  dX_ind_tab <- matrix(dX_ind[1:(N_S_io*N_R)], nrow = N_R, ncol = N_S_io, byrow = FALSE)
  
  dX_R_direct <- rep(0,N_R); dX_R_direct[Sim_region] <- Sim_sales
  dX_R_indirect <- rowSums(dX_tab)
  dX_R_induced <- rowSums(dX_ind_tab) - dX_R_indirect
  dX_R_total <- dX_R_direct + dX_R_indirect + dX_R_induced
  
  #Maps
  load("map_sf.RData")
  map_sf <- cbind(map_sf, dX_R_direct, dX_R_indirect, dX_R_induced, dX_R_total)
  
  png(filename = paste0("output/map_direct_", mode, ".png"), height = 700, width= 900)
  ggplot() +
    geom_sf(data = map_sf, mapping = aes(fill=dX_R_direct)) +
    scale_fill_distiller(palette = "Greys", direction = 1, name = "[m EUR]")+
    annotation_scale() +
    annotation_north_arrow(location = "br", which_north = "true") +
    theme(axis.title.y=element_blank(),axis.title.x=element_blank())
  dev.off()
  
  png(filename = paste0("output/map_indirect_", mode,".png"), height = 700, width= 900)
  ggplot() +
    geom_sf(data = map_sf, mapping = aes(fill=dX_R_indirect)) +
    scale_fill_distiller(palette = "Greys", direction = 1, trans = "log10", name = "[m EUR]")+
    annotate("text", x = 17, y = 49.5, label = paste0("M. Wrocław: ", round(dX_R_indirect[27], 2))) +
    annotate("segment", x = 17, y = 49.6, xend = 17, yend = 51.07, lwd = 1,
             arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
    annotate("text", x = 15, y = 50, label = paste0("Wrocławski: ", round(dX_R_indirect[31], 2))) +
    annotate("segment", x = 15, y = 50, xend = 16.6, yend = 51.1, lwd = 1,
             arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
    annotation_scale() +
    annotation_north_arrow(location = "br", which_north = "true") +
    theme(axis.title.y=element_blank(),axis.title.x=element_blank())
  dev.off()
  
  if(useDatarino == FALSE) {fname <- paste0("output/map_induced_noDR_", mode, ".png")}
  if(useDatarino == TRUE) {fname <- paste0("output/map_induced_DR_", mode, ".png")}
  png(filename = fname, height = 700, width= 900)
  ggplot() +
    geom_sf(data = map_sf, mapping = aes(fill=dX_R_induced)) +
    scale_fill_distiller(palette = "Greys", direction = 1, trans = "log10", name = "[m EUR]")+
    annotate("text", x = 17, y = 49.5, label = paste0("M. Wrocław: ", round(dX_R_induced[27], 2))) +
    annotate("segment", x = 17, y = 49.6, xend = 17, yend = 51.07, lwd = 1,
             arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
    annotate("text", x = 15, y = 50, label = paste0("Wrocławski: ", round(dX_R_induced[31], 2))) +
    annotate("segment", x = 15, y = 50, xend = 16.6, yend = 51.1, lwd = 1,
             arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
    annotation_scale() +
    annotation_north_arrow(location = "br", which_north = "true") +
    theme(axis.title.y=element_blank(),axis.title.x=element_blank())
  dev.off()
  
  if(useDatarino == FALSE) {fname <- paste0("output/map_total_noDR_", mode, ".png")}
  if(useDatarino == TRUE) {fname <- paste0("output/map_total_DR_", mode, ".png")}
  png(filename = fname, height = 700, width= 900)
  ggplot() +
    geom_sf(data = map_sf, mapping = aes(fill=dX_R_total)) +
    scale_fill_distiller(palette = "Greys", direction = 1, trans = "log10", name = "[m EUR]")+
    annotate("text", x = 17, y = 49.5, label = paste0("M. Wrocław: ", round(dX_R_total[27], 2))) +
    annotate("segment", x = 17, y = 49.6, xend = 17, yend = 51.07, lwd = 1,
             arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
    annotate("text", x = 15, y = 50, label = paste0("Wrocławski: ", round(dX_R_total[31], 2))) +
    annotate("segment", x = 15, y = 50, xend = 16.6, yend = 51.1, lwd = 1,
             arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
    annotation_scale() +
    annotation_north_arrow(location = "br", which_north = "true") +
    theme(axis.title.y=element_blank(),axis.title.x=element_blank())
  dev.off()
  
  #Simulation - confidence interval for the fraction of home indirect effect (PL, no drivetime)
  if (drivetime == FALSE) {
    N <- 1000 #4000 #how many replicas of the simulation
    NN <- nrow(MCMC_chain)
    NNN <- round(NN / N, 0)
    n_ind <- seq(from = 1, to = NN, by = NNN)
    dX_tab_NNN <- list()
    
    ####################
    for (rr in 1:length(n_ind)) {
      shape_vec_rr <- exp(MCMC_chain[n_ind[rr], 1:N_S_reg])
      scale_vec_rr <- exp(MCMC_chain[n_ind[rr], (N_S_reg+1):(2*N_S_reg)])
      theta_rr <- c(shape_vec_rr, scale_vec_rr)
      print(paste0("######## Iteration ", rr, " ##########"))
      A_R_ind_rr <- MRIO(theta_rr)
      A_R_rr <- A_R_ind_rr[1:(N_S_io*N_R), 1:(N_S_io*N_R)]
      A_R_sim_rr <- rbind(cbind(A_R_rr, as.matrix(Sim_Acol)),
                          matrix(0, nrow = 1, ncol = ncol(A_R)+1))
      rm(A_R_rr, A_R_ind_rr)
      dX_rr <- solve(diag(nrow(A_R)+1)-A_R_sim_rr, dY)
      dX_tab_NNN[[rr]] <- matrix(dX_rr[1:(N_S_io*N_R)], nrow = N_R, ncol = N_S_io, byrow = FALSE)
    }
    save(dX_tab_NNN, file = "output/dX_tab_NNN.RData")
    ##############
    
    load("output/dX_tab_NNN.RData")
    
    home_share_indirect <- rep(NA, length(dX_tab_NNN))
    for (rr in 1:length(dX_tab_NNN)) {
      dX_R_indirect_rr <- rowSums(dX_tab_NNN[[rr]])
      home_share_indirect[rr] <- dX_R_indirect_rr[Sim_region] / sum(dX_R_indirect_rr)
    }
    p_est <- dX_R_indirect[Sim_region] / sum(dX_R_indirect)
    png(filename = "output/histogram.png", height = 500, width = 500)
    hist(home_share_indirect, freq = FALSE, breaks = 10, xlab = "share of indirect effects on output in the impulse region",
         main = NULL, col = "darkred", ylab = "density")
    abline(v=p_est)
    text(x=p_est, y = 6.5, labels = "posterior means")
    q0025 = quantile(home_share_indirect, probs = 0.025)
    q0975 = quantile(home_share_indirect, probs = 0.975)
    q005 = quantile(home_share_indirect, probs = 0.05)
    q095 = quantile(home_share_indirect, probs = 0.95)
    dev.off()
  }
  
  #Bayes factors (drivetime / no drivetime)
  data2 <- list(data = data, hyper = hyper)
  log_posterior_wrap(theta0, data2)
  lb <- rep(-Inf, length(theta0))
  ub <- rep(Inf, length(theta0))
  names(lb) <- paste0("theta", 1:length(theta0))
  names(ub) <- paste0("theta", 1:length(theta0))
  colnames(MCMC_chainA) <- paste0("theta", 1:length(theta0))
  colnames(MCMC_chainB) <- paste0("theta", 1:length(theta0))
  combined.chains <- mcmc.list(mcmc(MCMC_chainA), 
                               mcmc(MCMC_chainB))
  Bridge <- bridge_sampler(samples = combined.chains,
                           data = data2,
                           log_posterior = log_posterior_wrap,
                           lb = lb,
                           ub = ub)
  if(drivetime == FALSE) {
    NoTimeBridge <- Bridge
    save(NoTimeBridge, file = "output/notimebridge.RData")
  } else {
    TimeBridge <- Bridge
    save(TimeBridge, file = "output/timebridge.RData")
  }
  load("output/notimebridge.RData")
  load("output/timebridge.RData")
  BF <- bayes_factor(TimeBridge, NoTimeBridge)$bf
  save(BF, file = "output/BF.RData")
}

