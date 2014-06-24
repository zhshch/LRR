###########################################################################################
### Code for "Discord Search in Multivariate Time-series Data Using Local Recurrence Rates" 
### Code Strcture (written in function name):                                              
###     discord (main function) 
###             ---- PAA
###             ---- Best.wsize
###             ---- Segment
###                     ---- MakeRP
###                             ---- phase_space_more (phase_space)
###                             ---- dist_m (B_dist)
###             ---- LREC
###             ---- Outlier
###                     ---- MakeRP1
###                             ---- phase_space_more (phase_space)
###                             ---- dist_m1 (B_dist)
### Usage:
###     Use discord() function to get the results
###     Input: time series data, embedding dimension, commpression ratio
###            window size, slack margin, top discords number
###            window size optimiazation(T/F)...
###     Output: change points, discords, number of distance queries, window size
###             *plots show the final results
###
### Recommended Experiments Environment:
### R version 3.1.0 (2014-04-10)
### Platform: i686-pc-linux-gnu (32-bit)
### attached base packages:
### [1] stats  graphics  grDevices utils  datasets  methods   base   
###
### Time: 2014/5/24 UTC+8 
###########################################################################################

discord <- function(dat,m=1,c=0.1,w=2,s=30,p=1,o=F){
        ## Main function
        ## x: multivarite time series data; m: embedding dimension
        ## c: commpression ratio; w: window size
        ## s: slack margin; p: top discords number
        ## o: window size optimiazation(T/F)
        #####################################################################
        x <- PAA(dat,c)
        if (o==T){
                w <- best.w(x,m,s)
        }
        sn <- 0
        N <- length(x[,1])
        cp_all <- Segment(x,m,w,s,makeplot=F,sn)
        o_all <- outlier(dat,N,c,p,sn,cp_all,makeplot=T)
        discords <- o_all[[1]]
        discords <- as.data.frame(discords)
        colnames(discords) <- c("Length","Mean","Begin","End")
        sn <- o_all[[2]]
        return(list(discords,sn,w))
}

PAA <- function(x,c=0.1){
        ## Piecewise Aggregate Approximation (PAA)
        ## input: 
        ## x: original time series data
        ## c: compression ratio
        ## output:
        ## x.new: time series data after compression by PAA
        s=1/c
        n <- nrow(x)
        col <- ncol(x)
        x.new <- array(0,dim=c(floor(n/s),col))
        for(i in 1:floor(n/s)){
                x.new[i,]<-c(as.matrix(apply(x[((i-1)*s):(i*s),],2,mean)))
        }
        return(x.new)
}

best.w <- function(x,m,s,begin=2,end=10){
        ## Function for window size optimiazation
        ## Input: 
        ##      x: datasets after PAA
        ##      m: embedding dimension
        ##      s: slack margin;
        ## Output:
        ##      best window size
        ###############################################
        N <- length(x[,1])
        res <- array(0)
        for(w in begin:end){
                cp_all <- Segment(x,m,w,s,makeplot=F,sn=0)
                if (length(cp_all[[1]])>2){
                        cp <- cp_all[[1]]
                        cp_entr_max_origin <- cp_all[[4]]
                        cp <- sort(cp)
                        cp_temp <- array(0)
                        cp_temp_entr <- array(0)
                        n <- length(cp)
                        j <- 1
                        for (i in 1:(n-1)){
                                if (cp[i+1]-cp[i]>10){
                                        cp_temp[j] <- cp[i]
                                        cp_temp_entr[j] <- cp_entr_max_origin[i]
                                        j <- j+1
                                }
                        }
                        cp <- cp_temp
                        cp_entr_max_origin <- cp_temp_entr
                        n <- length(cp)
                        cp_ratio_1 <- array(0)
                        cp_ratio_2 <- array(0)
                        cp_n <- array(0)
                        for (i in 1:(n-1)){
                                cp_n[i] <- cp[i+1]-cp[i]+1
                        }
                        res[w] <- sd(cp_n)
                }
                else{
                        res[w] <- NA
                }
                
        }
        res <- res[-1]
        plot(res,type="l")
        return(which(res==min(res,na.rm=T)[1])+1)
}

Segment <- function(x,m,w,s,makeplot=F,sn){
        ## Function for Segment Phase
        ## Input:
        ##      x: datasets after PAA
        ##      m: embedding dimension
        ##      s: slack margin;
        ##      sn: number of distance queries
        ## Output:
        ##      entr.max: change points after adjustment
        ##      D.01: local recurrence matrix
        ##      W.entr: Local Recurrence Rate (LREC)
        ##      entr.max.origin: change points before adjustment
        ##      sn: number of distance queries
        ##########################################################
        D <- MakeRP(x,t=0.75,step=1,m=1,sn,w=w+s)
        D.01 <- D[[1]]
        sn <- D[[2]]
        entr <- LREC(D.01,w,makeplot)
        entr.max <- entr[[1]]
        w.entr <- entr[[2]]
        entr.max.origin <- entr[[3]]
        return(list(entr.max,D.01,w.entr,entr.max.origin,sn))
}
LREC <- function(X,w,makeplot=T){
        ## Function for Local Recurrence Rate (LREC) extraction
        ## Input:
        ##      X: local recurrence matrix; w: window size
        ## Output:
        ##      entr.max: change points after adjustment
        ##      W.entr: Local Recurrence Rate (LREC)
        ##      entr.max.origin: change points before adjustment
        ##########################################################
        X.n <- nrow(X)
        W.entr <- array(0)
        ret <- array(1)
        W.m <- phase_space(data=c(1:X.n),N=X.n,m=w,tau=1)
        W.all <- w*w
        for (i in 1:ncol(W.m)){
                W.entr[i] <- sum(X[W.m[,i],W.m[,i]])/W.all
        }
        W.entr <- W.entr/(max(W.entr)-min(W.entr))
        ret <- ret*(W.entr>0)
        entr.w <- phase_space(data=ret,N=length(ret),m=2,tau=1)
        entr.new <- apply(entr.w,2,sd)
        entr.max <- which(entr.new==max(entr.new))+2
        entr.max.origin <- entr.max
        for (i in 1:length(entr.max)){
                if(ret[entr.max[i]+1]!=0){
                        entr.max[i] <- entr.max[i]+w
                }
        }
        
        if (makeplot==T){
                image(1:X.n,1:X.n,1-X)
        }
        return(list(entr.max,W.entr,entr.max.origin))
}

outlier <- function(dat,N,c,p,sn,cp_all,makeplot=T){
        ## Function for Outlier Dectection Phase
        ## Input:
        ##      dat: original datasets
        ##      N: number of records after PAA
        ##      m: embedding dimension
        ##      c: compression ratio
        ##      p: top discords number
        ##      sn: number of distance queries
        ## Output:
        ##      discord_show: discord vectors
        ##      sn: number of distance queries
        ##########################################################
        cp <- cp_all[[1]]
        cp <- sort(cp)
        cp_m <- cp_all[[2]]
        cp_entr <- cp_all[[3]]
        cp_entr_max_origin <- cp_all[[4]]
        sn <- cp_all[[5]]
        
        ## filtration of change points
        cp_temp <- array(0)
        cp_temp_entr <- array(0)
        n <- length(cp)
        j <- 1
        for (i in 1:(n-1)){
                if (cp[i+1]-cp[i]>=10){
                        cp_temp[j] <- cp[i]
                        cp_temp_entr[j] <- cp_entr_max_origin[i]
                        j <- j+1
                }
        }
        cp <- cp_temp
        cp_entr_max_origin <- cp_temp_entr
        
        
        n <- length(cp)
        cp_n <- array(0)
        cp_begin <- array(0)
        cp_end <- array(0)
        for (i in 1:(n-1)){
                cp_n[i] <- cp[i+1]-cp[i]+1
                cp_begin[i] <- cp[i]*(1/c)
                cp_end[i] <- cp[i+1]*(1/c)
        }
        
        n <- length(cp_entr_max_origin)
        cp <- cp_entr_max_origin
        cp_mean <- array(0)
        for (i in 1:(n-1)){
                cp_mean[i] <- mean(cp_entr[cp[i]:cp[i+1]])
        }
        cp_k <- cbind(cp_n,cp_mean,cp_begin,cp_end)
        rp1 <- MakeRP1(cp_k[,1:2],sn=sn)
        cp_rp <- rp1[[1]]
        sn <- rp1[[2]]
        deg <- sapply(1:nrow(cp_rp), function(i) {
                return(sum(cp_rp[i,]))
        })
        deg.percent <- deg/(n-1)
        image(cp_rp)
        plot(deg.percent)
        discord_factor <- array(1,dim=N*(1/c))
        if (max(deg.percent)>=0.4){
                dn <- order(deg.percent,decreasing=T)[1:p]
        }
        discord_show <- array(0,dim=c(length(dn),ncol(cp_k)))
        for (i in 1:length(dn)){
                discord_show[i,] <- cp_k[dn[i],]
        }    
        discord_factor[1:cp_begin[1]] <- 1
        for (i in 1:length(dn)){
                j <- dn[i]
                discord_factor[cp_begin[j]:cp_end[j]] <- 2   
        } 
        ## make plots
        if (makeplot){
                col <- ncol(dat)
                op <- par(mfrow = c(col+1, 1), mar = c(2, 4, 2, 2))
                on.exit(par(op))
                plot(cp_entr,type="l",xlab="Timestamps",ylab="LREC",lwd=2)
                
                for (k in 1:col){
                        plot(dat[,k],type="l",lwd=2,xlab="",ylab=paste("Variable ",as.character(k)))
                        for (i in 1:length(dn)){
                                j <- dn[i]           
                                lines(c(cp_begin[j]:cp_end[j]),dat[cp_begin[j]:cp_end[j],k],col="red",lwd=3)
                        }
                }
        }

        return(list(discord_show,sn))
}

MakeRP <- function(x,t=0.75,step=1,m=1,sn,w){
        ## Function for generating local recurrence matrix
        tau <- step
        x <- t(phase_space_more(data=x,m,tau))
        x.sc<-apply(x,2,function(x){
                y=(x-min(x))/(max(x)-min(x))+1
        })
        p<-ncol(x.sc)
        n<-nrow(x.sc)
        
        D <- dist_m(x.sc,w,sn)
        D1<- D[[1]]
        D2<- D[[2]]
        sn<- D[[5]]
        
        D1.01<-array(1,dim=dim(D1))
        D2.01<-array(1,dim=dim(D2))
        D1.theta<-quantile(D[[3]],probs=t)
        D2.theta<-quantile(D[[4]],probs=t)
        
        D1.01 <- D1.01*((D1-D1.theta)>0)
        D2.01 <- D2.01*((D2-D2.theta)>0)
        D.01 <- D1.01*D2.01
        return(list(D.01,sn))
}
dist_m <- function(x,w,sn){
        ## Function for distance computing
        ## d1: distance matrix by Euclidean distance
        ## d2: distance matrix by Bhattacharyya distance
        n=nrow(x)
        d1=d2=array(0,dim=c(n,n))
        d1_dia=d2_dia=array(0)
        dia_i=1
        for(i in 1:(n-w)){
                for (j in (i+1):(i+w)){
                        d1[i,j]=d1[j,i]=B_dist(x[i,],x[j,])
                        d2[i,j]=d2[j,i]=sqrt(sum((x[i,]-x[j,])^2))
                        d1_dia[dia_i]=d1[i,j]
                        d2_dia[dia_i]=d2[i,j]
                        sn = sn+1
                }
        }
        return(list(d1,d2,c(d1_dia,d1_dia),c(d2_dia,d2_dia),sn))
}

MakeRP1 <- function(x,t=0.75,step=1,m=1,sn){
        ## Function for generating full recurrence matrix
        tau <- step
        x <- t(phase_space_more(data=x,m=m,tau))
        x.sc<-apply(x,2,function(x){
                y=(x-min(x))/(max(x)-min(x))+1
        })
        p<-ncol(x.sc)
        n<-nrow(x.sc)
        
        D <- dist_m1(x.sc,sn)
        D1<-D[[1]]
        D2<-D[[2]]
        sn<-D[[3]]
        
        D1.01<-array(1,dim=dim(D1))
        D2.01<-array(1,dim=dim(D2))
        D1.theta<-quantile(D1,probs=t)
        D2.theta<-quantile(D2,probs=t)
        
        D1.01 <- D1.01*((D1-D1.theta)>0)
        D2.01 <- D2.01*((D2-D2.theta)>0)
        
        D.01=1-(1-D1.01)*(1-D2.01)
        return(list(D.01,sn))
}
dist_m1 <- function(x,sn){
        ## Function for distance computing
        ## d1: distance matrix by Euclidean distance
        ## d2: distance matrix by Bhattacharyya distance
        n=nrow(x)
        d1=d2=array(0,dim=c(n,n))
        for(i in 1:n){
                for (j in i:n){
                        d1[i,j]=d1[j,i]=B_dist(x[i,],x[j,])
                        d2[i,j]=d2[j,i]=sqrt(sum((x[i,]-x[j,])^2))
                        sn = sn+1
                }
        }
        return(list(d1,d2,sn))
}

B_dist<-function(x1,x2){
        ## Bhattacharyya distance
        sum1=sum(x1)
        sum2=sum(x2)
        sumup=sqrt(x1*x2)
        sumdown=sqrt(sum1*sum2)
        sumup=sum(sumup)
        sumdown=round(sumdown,3)
        sumup=round(sumup,3)
        if(sumup/sumdown==1) 
                return(0)
        else 
        {dist=sqrt(1-sumup/sumdown)
         return(dist)}
}

phase_space <- function(data,N,m,tau){
        ## Phase space reconstruction for 1-dimensional dataset
        M <- N-(m-1)*tau
        res <- array(0,dim=c(m,M))
        for (j in 1:M){
                for (i in 1:m){
                        res[i,j] <- data[(i-1)*tau+j]
                }
        }
        return(res)
}
phase_space_more <- function(data,m,tau){
        ## Phase space reconstruction for multi-dimensional dataset
        c <- ncol(data)
        N <- nrow(data)
        M <- N-(m-1)*tau
        res <- list()
        for (i in 1:c){
                res[[i]] <- phase_space(data[,i],N,m,tau)
        }
        ret <- array(0)
        for (i in 1:length(res)){
                ret <- rbind(ret,res[[i]])
        }
        return(ret[-1,])
}


