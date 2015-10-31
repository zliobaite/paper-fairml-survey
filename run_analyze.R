# 2015 07 21 I.Zliobaite
# discrimination testing

param_cik <- 100
param_n <- 10000
#param_pF <- 0.2
#param_pp <- 0.7
file_name_start <- 'results/out_'
param_table_cap <- 9

generate_data <- function(n,pF,pp,disc)
{
  nF <- round(pF*n)
  nM <- n - nF 
  
  npos <- round(pp*n)
    
  if (disc < 0)
  {
    do_reverse <- 1  
    disc <- -disc
  }else{
    do_reverse <- 0
  }
  
  data <- runif(n,0,1)
  s <- c(rep('F',nF),rep('M',nM)) 
  
  nswapF <- round(nF*disc)
  nswapM <- round(nM*disc)
  
  if ((nswapF + nswapM)>0)
  {
    ind_pick <- c(1:nswapF,(nF+1):(nF+nswapM))
    data_pick <- data[ind_pick]
    ind <- order(data_pick)
    ind_pick_sorted <- ind_pick[ind]  
    
    if (do_reverse)
    {
      s[ind_pick_sorted] <- c(rep('M',nswapM),rep('F',nswapF))
    }else{
      s[ind_pick_sorted] <- c(rep('F',nswapF),rep('M',nswapM))
    }
  }
  
  ind <- order(data,decreasing = TRUE)
  data <- data[ind]
  s <- s[ind]
  
  c <- rep(0,n)
  c[1:npos] <- 1
  
  data <- cbind(data,c,s)
  colnames(data) <- c('y','c','s')
  return(data)
}

compute_ent <-function(xx) #already a table
{
  ind <- which(xx>0)
  xx <- xx[ind]
  H <- -sum(log2(xx)*xx)
  return(H)
}

measure_disc <- function(data)
{
  n <- dim(data)[1]
  ps <- table(data[,'s'])/n
  pc <- table(data[,'c'])/n
  pjoint <- table(data[,'c'],data[,'s'])/n
  
  ddif <- pjoint['1','M']/ps['M'] - pjoint['1','F']/ps['F']
  
  if (ddif>0)
  {
    m1 <- pc['1']/ps['M']
    m2 <- pc['0']/ps['F']
  }else{
    m1 <- pc['1']/ps['F']
    m2 <- pc['0']/ps['M']
  }
  dmax <- min(m1,m2)
  ddnorm <- ddif/dmax
  
  dratio <- (pjoint['1','M']/ps['M'])/(pjoint['1','F']/ps['F'])
  delift <- (pjoint['1','M']/ps['M'])/pc['1']
  dolift <- (pjoint['1','M']/pjoint['1','F'])/(pjoint['0','M']/pjoint['0','F'])
  
  Hc <- compute_ent(pc)
  Hs <- compute_ent(ps)
  MI <- Hc + Hs - compute_ent(pjoint)
  MInorm <- MI / sqrt(Hc*Hs)
  
  dd <- cbind(ddif,ddnorm,dratio,delift,dolift,MInorm)
  return(dd)
}

for (pF in seq(0.1,0.9,0.2))
{
  for (pp in seq(0.1,0.9,0.2))
  {
    ddall <- c()
    for (disc in seq(-1,1,0.1))
    {
      dd <- c()
      for (sk2 in 1:param_cik)
      {
        data1 <- generate_data(param_n,pF,pp,disc)
        dd <- rbind(dd,measure_disc(data1))
      }
      dd <- round(apply(dd,2,mean),digits = 3)
      ddall <- rbind(ddall,c(disc,dd))
    }
    table_cap_neg <- -param_table_cap
    colnames(ddall)[1] <- 'disc'
    ind <- which(ddall==Inf)
    ddall[ind] <- param_table_cap
    ind <- which(ddall==-Inf)
    ddall[ind] <- table_cap_neg
    ind <- which(ddall>param_table_cap)
    ddall[ind] <- param_table_cap
    ind <- which(ddall<table_cap_neg)
    ddall[ind] <- table_cap_neg
    file_name_now <- paste(file_name_start,'pF',as.character(pF*100),'_pp',as.character(pp*100),'.dat',sep='')
    print(file_name_now)
    write.table(ddall, file = file_name_now, row.names = FALSE, col.names = TRUE, sep = ' ', quote = FALSE)
  }
}

