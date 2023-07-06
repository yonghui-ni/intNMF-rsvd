library(mclust)
library(NMF)
library(MASS)
library(iCluster)
library(InterSIM)
library(IntNMF)
library(rsvd)

#simulation settings;
#prop <- c(0.65,0.35)                  # true cluster number is 2
#prop <- c(0.30,0.40,0.30)            # true cluster number is 3
#prop <- c(0.20,0.30,0.27,0.23)       # true cluster number is 4
#prop <- c(0.23,0.17,0.25,0.15,0.2)   # true cluster number is 5
#size<-0.5
#size<-1
#size<-2
#size<-3
#size<-4

sim_intnmf_rsvd<-function(prop=c(0.65,0.35),size=0.5,iter=100){
        
        #clustering performance when k=2,3,4,5
        ARI.2<-c();ARI.3<-c();ARI.4<-c();ARI.5<-c()
        Purity.2<-c();Purity.3<-c();Purity.4<-c();Purity.5<-c()
        Entropy.2<-c();Entropy.3<-c();Entropy.4<-c();Entropy.5<-c()
        SI.2<-c();SI.3<-c();SI.4<-c();SI.5<-c()
        
        for (i in 1:iter){
                set.seed(i)
                prop <- prop               # sample proportion in clusters
                effect <- size             # cluster effect size
                sim.D <- InterSIM(n.sample=100,cluster.sample.prop=prop,delta.methyl=effect,delta.expr=effect,delta.protein=effect,p.DMP=0.25,p.DEG=NULL,p.DEP=NULL,
                                  do.plot=F, sample.cluster=T, feature.cluster=T)
                dat1 <- sim.D$dat.methyl
                dat2 <- sim.D$dat.expr
                dat3 <- sim.D$dat.protein
                true.cluster.assignment <- sim.D$clustering.assignment
                
                true.cluster <- true.cluster.assignment$cluster.id
                names(true.cluster) <- true.cluster.assignment$subjects
                
                d1.new <- dat1
                tmp1 <- rsvd(d1.new)
                d1.new <- d1.new %*% tmp1$v
                
                d2.new <- dat2
                tmp2 <- rsvd(d2.new)
                d2.new <- d2.new %*% tmp2$v
                
                d3.new <- dat3
                tmp3 <- rsvd(d3.new)
                d3.new <- d3.new %*% tmp3$v
                
                # shift the data to make all entries positive
                if (!all(d1.new>=0)) d1.new <- pmax(d1.new + abs(min(d1.new)), .Machine$double.eps) 
                if (!all(d2.new>=0)) d2.new <- pmax(d2.new + abs(min(d2.new)), .Machine$double.eps) 
                if (!all(d3.new>=0)) d3.new <- pmax(d3.new + abs(min(d3.new)), .Machine$double.eps) 
                
                # Globally rescale data to avoid potential overflow/underflow
                d1.new <- d1.new/max(d1.new)   
                d2.new <- d2.new/max(d2.new)
                d3.new <- d3.new/max(d3.new)
                dat.new <- list(d1.new,d2.new,d3.new)
                
                # integrative NMF
                fit.1 <- nmf.mnnals(dat=dat.new,k=2,maxiter=500,st.count=10,n.ini=10,ini.nndsvd=TRUE,seed=TRUE)
                fit.2 <- nmf.mnnals(dat=dat.new,k=3,maxiter=500,st.count=10,n.ini=10,ini.nndsvd=TRUE,seed=TRUE)
                fit.3 <- nmf.mnnals(dat=dat.new,k=4,maxiter=500,st.count=10,n.ini=10,ini.nndsvd=TRUE,seed=TRUE)
                fit.4 <- nmf.mnnals(dat=dat.new,k=5,maxiter=500,st.count=10,n.ini=10,ini.nndsvd=TRUE,seed=TRUE)
                
                
                ## Cluster Purity
                Purity.2[i]<-ClusterPurity(ComputedClusters=fit.1$clusters, TrueClasses=true.cluster)
                Purity.3[i]<-ClusterPurity(ComputedClusters=fit.2$clusters, TrueClasses=true.cluster)
                Purity.4[i]<-ClusterPurity(ComputedClusters=fit.3$clusters, TrueClasses=true.cluster)
                Purity.5[i]<-ClusterPurity(ComputedClusters=fit.4$clusters, TrueClasses=true.cluster)
                
                ## Cluster Entropy
                Entropy.2[i]<-ClusterEntropy(ComputedClusters=fit.1$clusters, TrueClasses=true.cluster)
                Entropy.3[i]<-ClusterEntropy(ComputedClusters=fit.2$clusters, TrueClasses=true.cluster)
                Entropy.4[i]<-ClusterEntropy(ComputedClusters=fit.3$clusters, TrueClasses=true.cluster)
                Entropy.5[i]<-ClusterEntropy(ComputedClusters=fit.4$clusters, TrueClasses=true.cluster)
                
                
                #Silhouette width
                SI.2[i]<-mean((silhouette(fit.1$clusters,dmatrix=1-fit.1$consensus))[,3])
                SI.3[i]<-mean((silhouette(fit.2$clusters,dmatrix=1-fit.2$consensus))[,3])
                SI.4[i]<-mean((silhouette(fit.3$clusters,dmatrix=1-fit.3$consensus))[,3])
                SI.5[i]<-mean((silhouette(fit.4$clusters,dmatrix=1-fit.4$consensus))[,3])
                
                
                ## Adjusted Rand Index
                ARI.2[i] <-adjustedRandIndex(fit.1$clusters,true.cluster)
                ARI.3[i]<-adjustedRandIndex(fit.2$clusters,true.cluster)
                ARI.4[i]<-adjustedRandIndex(fit.3$clusters,true.cluster)
                ARI.5[i]<-adjustedRandIndex(fit.4$clusters,true.cluster)
                
}
        ARI<-cbind(ARI.2,ARI.3,ARI.4,ARI.5)
        SI<-cbind(SI.2,SI.3,SI.4,SI.5)
        Purity<-cbind(Purity.2,Purity.3,Purity.4,Purity.5)
        Entropy<-cbind(Entropy.2,Entropy.3,Entropy.4,Entropy.5)
        
        return(cbind(ARI,SI,Purity,Entropy))
 
}


sim_intnmf<-function(prop=c(0.65,0.35),size=0.5,iter=100){
        
        #clustering performance when k=2,3,4,5
        ARI.2<-c();ARI.3<-c();ARI.4<-c();ARI.5<-c()
        Purity.2<-c();Purity.3<-c();Purity.4<-c();Purity.5<-c()
        Entropy.2<-c();Entropy.3<-c();Entropy.4<-c();Entropy.5<-c()
        SI.2<-c();SI.3<-c();SI.4<-c();SI.5<-c()
        
        for (i in 1:iter){
                set.seed(i)
                prop <- prop               # sample proportion in clusters
                effect <- size             # cluster effect size
                sim.D <- InterSIM(n.sample=100,cluster.sample.prop=prop,delta.methyl=effect,delta.expr=effect,delta.protein=effect,p.DMP=0.25,p.DEG=NULL,p.DEP=NULL,
                                  do.plot=F, sample.cluster=T, feature.cluster=T)
                dat1 <- sim.D$dat.methyl
                dat2 <- sim.D$dat.expr
                dat3 <- sim.D$dat.protein
                true.cluster.assignment <- sim.D$clustering.assignment
                
                true.cluster <- true.cluster.assignment$cluster.id
                names(true.cluster) <- true.cluster.assignment$subjects
                
        
                ## Make all data positive by shifting to positive direction.
                ## Also rescale the datasets so that they are comparable. 
                if (!all(dat1>=0)) dat1 <- pmax(dat1 + abs(min(dat1)), .Machine$double.eps) 
                dat1 <- dat1/max(dat1)   
                if (!all(dat2>=0)) dat2 <- pmax(dat2 + abs(min(dat2)), .Machine$double.eps) 
                dat2 <- dat2/max(dat2)
                if (!all(dat3>=0)) dat3 <- pmax(dat3 + abs(min(dat3)), .Machine$double.eps) 
                dat3 <- dat3/max(dat3)
                
                dat <- list(dat1,dat2,dat3)
                
                # integrative NMF
                fit.1 <- nmf.mnnals(dat=dat,k=2,maxiter=500,st.count=10,n.ini=10,ini.nndsvd=TRUE,seed=TRUE)
                fit.2 <- nmf.mnnals(dat=dat,k=3,maxiter=500,st.count=10,n.ini=10,ini.nndsvd=TRUE,seed=TRUE)
                fit.3 <- nmf.mnnals(dat=dat,k=4,maxiter=500,st.count=10,n.ini=10,ini.nndsvd=TRUE,seed=TRUE)
                fit.4 <- nmf.mnnals(dat=dat,k=5,maxiter=500,st.count=10,n.ini=10,ini.nndsvd=TRUE,seed=TRUE)
                
                
                ## Cluster Purity
                Purity.2[i]<-ClusterPurity(ComputedClusters=fit.1$clusters, TrueClasses=true.cluster)
                Purity.3[i]<-ClusterPurity(ComputedClusters=fit.2$clusters, TrueClasses=true.cluster)
                Purity.4[i]<-ClusterPurity(ComputedClusters=fit.3$clusters, TrueClasses=true.cluster)
                Purity.5[i]<-ClusterPurity(ComputedClusters=fit.4$clusters, TrueClasses=true.cluster)
                
                ## Cluster Entropy
                Entropy.2[i]<-ClusterEntropy(ComputedClusters=fit.1$clusters, TrueClasses=true.cluster)
                Entropy.3[i]<-ClusterEntropy(ComputedClusters=fit.2$clusters, TrueClasses=true.cluster)
                Entropy.4[i]<-ClusterEntropy(ComputedClusters=fit.3$clusters, TrueClasses=true.cluster)
                Entropy.5[i]<-ClusterEntropy(ComputedClusters=fit.4$clusters, TrueClasses=true.cluster)
                
                
                #Silhouette width
                SI.2[i]<-mean((silhouette(fit.1$clusters,dmatrix=1-fit.1$consensus))[,3])
                SI.3[i]<-mean((silhouette(fit.2$clusters,dmatrix=1-fit.2$consensus))[,3])
                SI.4[i]<-mean((silhouette(fit.3$clusters,dmatrix=1-fit.3$consensus))[,3])
                SI.5[i]<-mean((silhouette(fit.4$clusters,dmatrix=1-fit.4$consensus))[,3])
                
                
                ## Adjusted Rand Index
                ARI.2[i] <-adjustedRandIndex(fit.1$clusters,true.cluster)
                ARI.3[i]<-adjustedRandIndex(fit.2$clusters,true.cluster)
                ARI.4[i]<-adjustedRandIndex(fit.3$clusters,true.cluster)
                ARI.5[i]<-adjustedRandIndex(fit.4$clusters,true.cluster)
                
        }
        ARI<-cbind(ARI.2,ARI.3,ARI.4,ARI.5)
        SI<-cbind(SI.2,SI.3,SI.4,SI.5)
        Purity<-cbind(Purity.2,Purity.3,Purity.4,Purity.5)
        Entropy<-cbind(Entropy.2,Entropy.3,Entropy.4,Entropy.5)
        
        return(cbind(ARI,SI,Purity,Entropy))
        
}

#simulation result example
result_intnmf_rsvd<-sim_intnmf_rsvd(prop=c(0.30,0.40,0.30),size=2,iter=5)
result_intnmf<-sim_intnmf(prop=c(0.30,0.40,0.30),size=2,iter=5)
