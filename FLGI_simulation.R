#Function: FLGIsimulation (and FLGIallocationprobability)
#Purpose: to simulate a clinical trial using the FLGI design 

#Inputs:
#K = number of treatment arms (including control) with
#J = number of patient blocks
#b = block size (J*b=N, the total number of patients in trial)
#fisherco = cutoff point used for the Fisher exact tests
#p = vector of success probabilities for each treatment
#sims = number of Monte Carlo simulations
#prior = matrix (K x 2) containing the prior parameters for each treatment.
#        use 1s for an uninformative prior
#mc = number of Monte Carlo simulations used to calculate the FLGI probabilities

#Outputs:
#fwer = family wise type 1 error rate
#power = power of the trial
#ENS = expected number of successes
#ENSsd = standard deviation of trial successess
#EAP = expected allocation proportion
#EAPsd = standard deviation of the allocation to the best treatment
#avgss = vector of average sample sizes for each treatment arm
#Succ = vector of number of successes in each trial
#(option to add further outputs if neccersary)

#Requires: 
Gittins_mat <-read.table("Gittins_mat.txt")

#Prior should be a K x 2 matrix of 1's. This represents uninformative priors for each treatment
FLGIsimulation <- function(K,J,b,fisherco,p,sims,prior,mc){
      
      if (K != dim(prior)[1]){
            stop("Prior size does not match number of treatments")
      }
      
      if (length(p) != K){
            stop("Probability vector does not match treatments")
      }
      
      treatments <- c(0:(K-1)) #Treatments index
      bcutoff <- fisherco/(K-1) #Adjusted cut-off point (Bonferroni)
      N <- J*b #Total sample size
      
      #Create matrix to store trial outcomes
      results <- matrix(, nrow = sims, ncol = K-1) 
      
      #Create vectors to store total successess and failures in each trial
      Succ <- vector()
      Fail <- vector()
      
      #Create matrix to store allocation numbers to calculate EAP
      ss <- matrix(, nrow = sims, ncol = K)
      
      #Create matrix to store all p values
      pvalues <- matrix(,nrow = sims,ncol = (K-1))
      
      #Create list to store all data
      bigdata <- list()
      
      #Allocation to control
      allocontrol <- vector()
      
      for (i in 1:sims){
            
            s0 <- prior[,1] #should equal 1 with an uninformative prior
            f0 <- prior[,2] #should equal 1 with an uninformative prior
            
            n <- matrix(0, nrow = J, ncol = K) #To store the sample sizes on each arm for each block
            
            s <- matrix(0, nrow = J, ncol = K) #Successes
            f <- matrix(0, nrow = J, ncol = K) #Failures
            
            #Assume that patient delay is the same length as a block period.
            #Hence first 2 blocks use equal allocation probabilities
            
            for (j in 1:2){ #First 2 blocks where no data is available
                  
                  allocation <- sample(treatments, b, replace = TRUE)
                  
                  for (k in 1:K){
                        n[j,k] <- sum(allocation==(k-1))
                  }
                  
                  for (k in 1:K){
                        s[j,k] <- rbinom(1,n[j,k],p[k])
                        f[j,k] <- n[j,k] - s[j,k]
                  }
                  
            }
            
            for (j in 3:3){
                  
                  allocation <- vector()
                  
                  #make posterior by adding inital value and observed successes/failures
                  #Only have access to block-2 and older data due to response delay
                  posterior <- matrix(cbind(s0+s[1,], f0 + f[1,]),
                                      byrow = FALSE, ncol = 2)
                  
                  #calculate probabilities for rest of block
                  prob <- FLGIallocationprobability(posterior, b, mc)
                  
                  #allocate remaining patients in block
                  allocation <- sample(treatments, b, replace = TRUE, prob = prob)
                  for (k in 1:K){
                        n[j,k] <- sum(allocation==(k-1))
                  }
                  
                  for (k in 1:K){
                        s[j,k] <- rbinom(1,n[j,k],p[k])
                        f[j,k] <- n[j,k] - s[j,k]
                  }
            }
            
            if (J>3){
            for (j in 4:J){
                  
                  allocation <- vector()
                  
                  #make posterior by adding inital value and observed successes/failures
                  #Only have access to block-2 and older data due to response delay
                  posterior <- matrix(cbind(s0+colSums(s[1:(j-2),]), f0 + colSums(f[1:(j-2),])),
                                      byrow = FALSE, ncol = 2)
                  
                  #calculate probabilities for rest of block
                  prob <- FLGIallocationprobability(posterior, b, mc)

                  #allocate remaining patients in block
                  allocation <- sample(treatments, b, replace = TRUE, prob = prob)
                  for (k in 1:K){
                        n[j,k] <- sum(allocation==(k-1))
                  }
                  
                  for (k in 1:K){
                        s[j,k] <- rbinom(1,n[j,k],p[k])
                        f[j,k] <- n[j,k] - s[j,k]
                  }
            }
            }
            
            allocontrol[i] <- colSums(n)[1]/N
            
            ss[i,] <- colSums(n)
            
            data <- rbind(colSums(s),colSums(f)) 
            
            bigdata[[i]] <- data

            Succ[i] <- sum(data[1,]) #Store total successes
            Fail[i] <- sum(data[2,]) #Store total failures
            
            #Perform Fisher's exact test with adjusted cut-off
            reject <- vector()
            for (k in 1:(K-1)){
                  pval <- fisher.test(data[,c(k+1,1)], alternative="greater")$p.value
                  if (pval < bcutoff){
                        reject[k] <- 1 #Denote 1 as rejecting the null
                  }
                  else{
                        reject[k] <- 0 #Fail to reject null
                  }
                  pvalues[i,k] <- pval
            }
            results[i,] <- reject
            
      }
      
      #Summarise weather there has been a type I error in each simulation
      temp <- vector()
      for (i in 1:sims){
            if (sum(results[i,])==0){
                  temp[i] <- 0 
            }
            else {
                  temp[i] <- 1
            }
      }
      results <- cbind(results, temp)
      
      #Calculate the family wise error rate
      fwer <- sum(results[,K])/sims
      
      #Find the best treatment and calculate power from there (does not
      #include multiple best treatments - only for unique)
      best <- which.max(p[c(1:K)])
      power <- sum(results[,(best-1)])/sims
      
      #Calculate ENS and ENF
      ENS <- mean(Succ)
      ENSsd <- sd(Succ)
      ENF <- mean(Fail)
      ENFsd <- sd(Fail)
      
      #Calculate EFP
      EFP <- mean(Fail/N)
      EFPsd <- sd(Fail/N)
      
      #Calculate EAP
      EAP <- mean(ss[,best]/N)
      EAPsd <- sd(ss[,best]/N)
      
      #Average ss on each arm
      avgss <- colMeans(ss)
      
      return <- list(fwer=fwer,
                     power=power,
                     ENS=ENS,
                     ENSsd=ENSsd,
                     EAP=EAP,
                     EAPsd=EAPsd,
                     avgss=avgss,
                     Succ=Succ)
      
}

#This function is modified from MATLAB code from supplementary materials of the following paper:
# Villar, S. S., Wason, J., and Bowden, J. (2015). Response-adaptive randomization for multiarm
# clinical trials using the forward looking gittins index rule. Biometrics, 71(4):969-978.

FLGIallocationprobability <- function(I0, block, noRuns){
      
      
      K = dim(I0)[1] # Set the number of treatments
      # set dimension of index vector and other vectors
      index <- matrix(0, K, 1)
      selected <- matrix(0, noRuns, block)
      prob <- matrix(0, K, block)
      
      for (j in 1:noRuns){
            n = rowSums(I0)
            s= I0[,1]
            f= I0[,2]
            for (t in 0:(block-1)){
                  for (k in 1:K){
                        index[k] = Gittins_mat[s[k],f[k]]
                  }
                  max_index <- max(index)
                  kmax <- which.max(index)
                  idx <- which(index==max(index))
                  
                  #Randomise if more than one treatment attains the maximum GI
                  
                  if (length(idx) > 1){
                        posi <- sample(1:length(idx),1)
                        kmax <- idx[posi]
                        max_index <- index[kmax]
                  }
                  
                  selected[j,t+1]= kmax
                  
                  snext=s
                  fnext=f
                  nnext=n #passive state dynamics
                  
                  probSuc_kmax = s[kmax]/(s[kmax]+f[kmax]) #prob. of a success on current patient with treatment nmax
                  
                  Pos = (runif(1) <= probSuc_kmax) #simulate a success event
                  
                  nnext[kmax]= n[kmax] + 1
                  
                  if (Pos==1){
                        snext[kmax] = s[kmax] + 1
                  } else{
                        fnext[kmax] = f[kmax]+1
                  }
                  
                  
                  s = snext
                  f = fnext
                  n = nnext
            }
      }
      
      for (i in 1:block){
            for (k in 1:K){
                  prob[k,i]= sum(selected[,i]==k)/noRuns
            }
      }
      allocation_probabilities = rowMeans(prob)
      allocation_probabilities
}  
