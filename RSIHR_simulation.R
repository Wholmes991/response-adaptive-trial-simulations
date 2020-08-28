#Function: RSIHRsimulation
#Purpose: to simulate a clinical trial using the RSIHR design 
#         (only for K=2 treatments)

#Inputs:
#K = number of treatment arms (including control) with
#J = number of patient blocks
#b = block size (J*b=N, the total number of patients in trial)
#fisherco = cutoff point used for the Fisher exact tests
#p = vector of success probabilities for each treatment
#sims = number of Monte Carlo simulations

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

RSIHRsimulation <- function(K,J,b,fisherco,p,sims){

      
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
      
      #Create matrix to store p values from fishers exact test
      pvalues <- matrix(,nrow = sims,ncol = (K-1))
      
      #Create list to store all data
      bigdata <- list()
      
      for (i in 1:sims){
            
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
            
            for (j in 3:J){
                  
                  allocation <- vector()
                  
                  #estimate ps
                  phat <- vector()
                  for (k in 1:K){
                        phat[k] <- (1+sum(s[1:(j-2),k]))/(2+sum(n[1:(j-2),k]))
                  }
                  
                  prob <- sqrt(phat[1])/(sqrt(phat[1])+sqrt(phat[K]))
                  
                  prob <- c(prob,1-prob)
                  
                  
                  #allocate patients in block
                  allocation <- sample(treatments, b, replace = TRUE, prob = prob)
                  for (k in 1:K){
                        n[j,k] <- sum(allocation==(k-1))
                  }
                  
                  for (k in 1:K){
                        s[j,k] <- rbinom(1,n[j,k],p[k])
                        f[j,k] <- n[j,k] - s[j,k]
                  }
            }
            
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
            
            print(i)      
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
