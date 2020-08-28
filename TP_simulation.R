#Function: TPsimulation
#Purpose: to simulate a clinical trial using the TP adapted for binary responses 

#Inputs:
#K = number of treatment arms (including control) with
#J = number of patient blocks
#b = block size (J*b=N, the total number of patients in trial)
#fisherco = cutoff point used for the Fisher exact tests
#p = vector of success probabilities for each treatment
#sims = number of Monte Carlo simulations
#a1 and a2 = vectors of length K representing prior parameters - use 1s for uninformative
#mc = number of Monte Carlo simulations used to calculate the TP probabilities

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


#a1 and a2 are the prior parameters and can be centered around 
#estimate of p0 (a1=p0, a2=1-p0)
TPsimulation <- function(K,J,b,fisherco,p,sims,a1,a2,mc){
      
      if ((K != length(a1)) | (K!=length(a2))){
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
            
            for (j in 3:3){ #Split j=3 up due to colSums restrictions
                  
                  allocation <- vector()
                  
                  #make posterior by adding inital value and observed successes/failures
                  #Only have access to block-2 and older data due to response delay
                  posterior <- matrix(cbind(a1 + s[1,], a2 + f[1,]),
                                      byrow = FALSE, ncol = 2)
                  
                  
                  #calculate probabilities of pk > p0 through monte carlo simulations
                  montecarlo <- matrix(, nrow = mc, ncol=K)
                  prob <- vector()
                  for (k in 1:K){
                        montecarlo[,k] <- rbeta(mc,posterior[k,1],posterior[k,2])
                  }
                  
                  for (k in 1:(K-1)){
                        prob[k] <- mean(montecarlo[,k+1]>montecarlo[,1])^(10*(j/J)^0.75)
                  }

                  if (sum(prob)!=0){
                     prob <- prob/sum(prob)
                  }
                                    
                  prob0 <- (1/(K-1))*exp(max(n[1,2:K]) - n[1,1])^(0.25*(j/J)) 
                  
                  prob <- c(prob0,prob)
                  
                  prob <- prob/sum(prob)

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
            
            if (J>3){
            for (j in 4:J){
                  
                  allocation <- vector()
                  
                  #make posterior by adding inital value and observed successes/failures
                  #Only have access to block-2 and older data due to response delay
                  posterior <- matrix(cbind(a1 + colSums(s[1:(j-2),]), a2 + colSums(f[1:(j-2),])),
                                      byrow = FALSE, ncol = 2)
                  
                  
                  #calculate probabilities of pk > p0 through monte carlo simulations
                  montecarlo <- matrix(, nrow = mc, ncol=K)
                  prob <- vector()
                  for (k in 1:K){
                        montecarlo[,k] <- rbeta(mc,posterior[k,1],posterior[k,2])
                  }
                  
                  for (k in 1:(K-1)){
                        prob[k] <- mean(montecarlo[,k+1]>montecarlo[,1])^(10*(j/J)^0.75)
                  }

                  if (sum(prob)!=0){
                     prob <- prob/sum(prob)
                  }
                  
                  prob0 <- (1/(K-1))*exp(max(colSums(n[1:(j-2),2:K])) - sum(n[1:(j-2),1]))^(0.25*(j/J))  
                  
                  prob <- c(prob0,prob)
                  
                  prob <- prob/sum(prob)

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
