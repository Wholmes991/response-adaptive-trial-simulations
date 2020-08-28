#Function: EAsimulation
#Purpose: to simulate a clinical trial using Equal Allocation (EA)

#Inputs:
#K = number of treatment arms (including control) with
#N = number of patients in the trial
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

EAsimulation <- function(K,N,fisherco,p,sims){
      treatments <- c(0:(K-1)) #Treatments
      bcutoff <- fisherco/(K-1) #Adjusted cut-off point (Bonferroni correction)
      
      #Create matrix to store trial outcomes
      results <- matrix(, nrow = sims, ncol = K-1) 
      
      #Create vectors to store total successess and failures in each trial
      Succ <- vector()
      Fail <- vector()
      
      #Create matrix to store allocation numbers to calculate EAP
      ss <- matrix(, nrow = sims, ncol = K)
      
      #Create a pvalue matrix to store all pvalues
      pvalues <- matrix(,nrow = sims,ncol = (K-1))
      
      #Create list to store all data
      bigdata <- list()
      
      #Simulate trial sims times
      for (i in 1:sims){
            
            #Allocate the N patients to each treatment with equal probability
            allocation <- sample(treatments, N, replace = TRUE)
            
            #Denote the number of patients on each arm in the vector n
            n <- vector()
            for (k in 1:K){
                  n[k] <- sum(allocation==(k-1))
            }
            ss[i,] <- n #Store sample sizes for each trial
            
            #Simulate the trial and record the successes and failure on each arm
            s <- vector()
            f <- vector()
            for (k in 1:K){
                  s[k] <- rbinom(1,n[k],p[k])
                  f[k] <- n[k] - s[k]
            }
            
            data <- rbind(s,f)
            Succ[i] <- sum(data[1,]) #Store total successes
            Fail[i] <- sum(data[2,]) #Store total failures
            
            bigdata[[i]] <- data
            
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

