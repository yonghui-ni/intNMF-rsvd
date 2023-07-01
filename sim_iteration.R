### new proposed method
#setwd("~/Documents/paper_clustering")
rm(list=ls())
library(mclust)
#library(RColorBrewer)
library(NMF)
library(MASS)
library(iCluster)
library(InterSIM)
library(IntNMF)
#library(SNFtool)
#require(graphics)
library(rsvd)


#Purity<-c()
#Entropy<-c()
#SI<-c()
ARI<-c()
#prop<-c(0.65,0.35)
#prop <- c(0.30,0.40,0.30)
prop <- c(0.20,0.30,0.27,0.23)
size<-3

system.time(for (i in 1:100){
        set.seed(i)
        prop <- prop
        effect <- size
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
        #k<-choose.rank(tmp1$d)
        d1.new <- d1.new %*% tmp1$v
        
        d2.new <- dat2
        tmp2 <- rsvd(d2.new)
        #k<-choose.rank(tmp2$d)
        d2.new <- d2.new %*% tmp2$v
        
        d3.new <- dat3
        tmp3 <- rsvd(d3.new)
        #k<-choose.rank(tmp3$d)
        d3.new <- d3.new %*% tmp3$v
        
        if (!all(d1.new>=0)) d1.new <- pmax(d1.new + abs(min(d1.new)), .Machine$double.eps) # shift the data to make all entries positive
        if (!all(d2.new>=0)) d2.new <- pmax(d2.new + abs(min(d2.new)), .Machine$double.eps) 
        if (!all(d3.new>=0)) d3.new <- pmax(d3.new + abs(min(d3.new)), .Machine$double.eps) 
        d1.new <- d1.new/max(d1.new)   # Globally rescale data to avoid potential overflow/underflow
        d2.new <- d2.new/max(d2.new)
        d3.new <- d3.new/max(d3.new)
        dat.new <- list(d1.new,d2.new,d3.new)
        
        ### Use optimum k = 4 and find cluster membership 
        #fit.1.new<- nmf.mnnals(dat=d1.new,k=length(prop),maxiter=200,st.count=10,n.ini=10,ini.nndsvd=TRUE,seed=TRUE)
        #fit.2.new<- nmf.mnnals(dat=d2.new,k=length(prop),maxiter=200,st.count=10,n.ini=10,ini.nndsvd=TRUE,seed=TRUE)
        #fit.3.new<- nmf.mnnals(dat=d3.new,k=length(prop),maxiter=200,st.count=10,n.ini=10,ini.nndsvd=TRUE,seed=TRUE)
        fit.new <- nmf.mnnals(dat=dat.new,k=length(prop),maxiter=400,st.count=10,n.ini=10,ini.nndsvd=TRUE,seed=TRUE)
        
        #cluster.mem.nmf.1.new <- fit.1.new$clusters
        #cluster.mem.nmf.2.new <- fit.2.new$clusters
        #cluster.mem.nmf.3.new <- fit.3.new$clusters
        cluster.mem.nmf.new <- fit.new$clusters
        
        
        ## Cluster Purity
        #Purity[i]<-ClusterPurity(ComputedClusters=cluster.mem.nmf.new, TrueClasses=true.cluster) 
        ## Cluster Entropy
        #Entropy[i]<-ClusterEntropy(ComputedClusters=cluster.mem.nmf.new, TrueClasses=true.cluster)
        ## Consensus matrix plot
        #ConsensusMatPlot(fit.new,rowLab=TRUE,colLab=TRUE)
        ## Silhouette plot
        #SilhouettePlot(fit.new,cluster.col=NULL)
        #SI[i]<-mean((silhouette(cluster.mem.nmf.new,dmatrix=1-fit.new$consensus))[,3])
        #SI[i]<-mean(CancerSubtypes::silhouette_SimilarityMatrix(cluster.mem.nmf.new,fit.new$consensus)[,3])
        #adjustedRandIndex
        ARI[i] <-adjustedRandIndex(cluster.mem.nmf.new,true.cluster)
})

#save(Purity,Entropy,SI,ARI,file='sim_new_iteration_ce1.RData')
#save(Purity,Entropy,SI,ARI,file='sim_new_iteration_ce0.5.RData')
#save(Purity,Entropy,SI,ARI,file='sim_new_iteration_ce2.5.RData')
#save(Purity,Entropy,SI,ARI,file='sim_new_iteration_ce2.RData')

#save(ARI,file='sim100_new_iteration_ce2.RData')
#save(ARI,file='sim100_new_iteration_ce1.RData')
#save(ARI,file='sim100_new_iteration_ce3.RData')
#save(ARI,file='sim100_new_iteration_ce0.5.RData')

save(ARI,file='sim100_new_iteration_ce2_1.RData')
save(ARI,file='sim100_new_iteration_ce1_1.RData')
save(ARI,file='sim100_new_iteration_ce3_1.RData')
save(ARI,file='sim100_new_iteration_ce0.5_1.RData')


#save(ARI,file='sim100k2_new_iteration_ce0.5.RData')
#save(ARI,file='sim100k2_new_iteration_ce1.RData')
#save(ARI,file='sim100k2_new_iteration_ce2.RData')
#save(ARI,file='sim100k2_new_iteration_ce3.RData')

#save(ARI,file='sim100k3_new_iteration_ce0.5.RData')
#save(ARI,file='sim100k3_new_iteration_ce1.RData')
#save(ARI,file='sim100k3_new_iteration_ce2.RData')
#save(ARI,file='sim100k3_new_iteration_ce3.RData')
