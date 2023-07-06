## Take AML data for methods application example of optimal cluster number and cluster assignment
## load('AML_1000_norm.RData')
## dat1 expression: RNAseqv2  level 3 RSEM genes normalized; log2 transformed 
## dat2 methylation: Illumina-450k level 3
## dat3 miRNA: Illumina mirnaseq level 3 miR gene expressio; log2 transformed 

#########################################################
#                      IntNMF-RSVD                      #
#########################################################
library(IntNMF)
library(rsvd)
load('AML_1000_norm.RData')
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

# obtain optimal cluster number
opt.k.new <- nmf.opt.k(dat=dat.new)      
k.fast<- which.max(rowMeans(opt.k.new))+1 

fit.new <- nmf.mnnals(dat=dat.new,k=k.fast)
cls.fast<-fit.new$clusters


#########################################################
#                        IntNMF                         #
#########################################################
load('AML_1000_norm.RData')
# shift the data to make all entries positive
if (!all(dat1>=0)) dat1 <- pmax(dat1 + abs(min(dat1)), .Machine$double.eps) 
dat1 <- dat1/max(dat1)   
if (!all(dat2>=0)) dat2 <- pmax(dat2 + abs(min(dat2)), .Machine$double.eps) 
dat2 <- dat2/max(dat2)
if (!all(dat3>=0)) dat3 <- pmax(dat3 + abs(min(dat3)), .Machine$double.eps) 
dat3 <- dat3/max(dat3)

dat.intnmf <- list(dat1,dat3,dat2)

# obtain optimal cluster number
opt.k <- nmf.opt.k(dat=dat.intnmf)
k.intnmf<-which.max(rowMeans(opt.k))+1 

fit.intnmf <- nmf.mnnals(dat=dat.intnmf, k=k.intnmf,ini.nndsvd=TRUE, seed=TRUE)
cls.intnmf<-fit.intnmf$clusters

#########################################################
#                        iCluster                       #
#########################################################
library(iCluster)
load('AML_1000_norm.RData')
dat <- list(dat1,dat2,dat3)

# obtain optimal cluster number
fit<-list()
pod<-list()
for (k in 2:8) {
        fit[[k]]<-iCluster(dat, k=k, lambda=c(0.1,0.1,0.1))
        pod[[k]]<-compute.pod(fit[[k]])
}
pod<-unlist(pod)
k.icluster<-which.min(pod)+1

fit.icluster<-iCluster(dat, k=k.icluster, lambda=c(0.1,0.1,0.1,0.1))
cls.icluster<-fit.icluster$clusters


#########################################################
#                        moCluster                      #
#########################################################
library(mogsa)
load('AML_1000_norm.RData')
dat_t<-list(t(dat1),t(dat2),t(dat3))

moa<-mbpca(dat_t, ncomp=10, method='globalScore', k = "all", center = TRUE,
           scale = FALSE, option = "uniform", maxiter = 1000,
           moa = TRUE, verbose = TRUE, svd.solver ="fast.svd",
           k.obs = "all", w = NA, w.obs = NA,
           unit.p = FALSE, unit.obs = FALSE, pos = FALSE)      
system(bootMbpca(moa, mc.cores = 1, B = 20, replace = TRUE,
                 resample = c("sample", "gene", "total"), log = "y", ncomp = NULL, method = NULL,
                 maxiter = 1000, svd.solver = c("svd", "fast.svd", "propack"), plot = TRUE))
# black dot significantly higher than the permutation eigenvalues (red box) represents the concordant structures across the data sets

get.elbow <- function(values, is.max) {
        second.derivatives = c()
        for (i in 2:(length(values) - 1)) {
                second.derivative = values[i + 1] + values[i - 1] - 2 * values[i]
                second.derivatives = c(second.derivatives, second.derivative)
        }
        print(second.derivatives)
        if (is.max) {
                return(which.max(second.derivatives) + 1)
        } else {
                return(which.min(second.derivatives) + 1)
        }
}
ncomp = get.elbow(moa@eig, is.max=T) 

moa<-mbpca(dat_t, ncomp=ncomp, method='globalScore', k = "all", center = TRUE,
           scale = FALSE, option = "uniform", maxiter = 1000,
           moa = TRUE, verbose = TRUE, svd.solver ="fast.svd",
           k.obs = "all", w = NA, w.obs = NA,
           unit.p = FALSE, unit.obs = FALSE, pos = FALSE)
  
scr <- moaScore(moa)
# obtain optimal cluster number
gap <- moGap(moa, K.max = 8, cluster = "kmeans")
k1<-gap$nClust[1];k2<-gap$nClust[3];if (k1>1) {k.mo<-k1} else {k.mo<-k2} 

kmean.mocluster <- kmeans(scr,k.mo)
cls.mocluster <- kmean.mocluster$cluster

#########################################################
#                        PINSPlus                       #
#########################################################
library(PINSPlus)
load('AML_1000_norm.RData')
dat <- list(dat1,dat2,dat3)
result.pins <- SubtypingOmicsData(dataList = dat,clusteringMethod = "pam")

k.pins<-max(unique(result.pins$cluster2))
cls.pins<-result.pins$cluster2


#########################################################
#                        SNF                            #
#########################################################
load('AML_1000_norm.RData')
library(SNFtool)
dat1<-standardNormalization(dat1)
dat2<-standardNormalization(dat2)
dat3<-standardNormalization(dat3)
K = round(ncol(dat1)/10)    # number of neighbors, usually (10~30)
alpha = 0.5  	            # hyperparameter, usually (0.3~0.8)
T = 20 	
dat <- list(dat1,dat2,dat3)

dist <- lapply(dat, function(x) (dist2(x, x)^1/2))          #distance matrix
A<-lapply(dist, function(x) affinityMatrix(x, K, alpha))    #similarity matrix
W<-SNF(A, K, T)                                             #fused similarity matrix

# obtain optimal cluster number
estimationResult = estimateNumberOfClustersGivenGraph(W, 2:8) 
k.snf<-estimationResult$`Eigen-gap best`
cls.snf <- spectralClustering(W,k.snf)



#########################################################
#                        NEMO                           #
#########################################################
load('AML_1000_norm.RData')
library(SNFtool)                   # standardNormalization()
dat1<-standardNormalization(dat1)
dat2<-standardNormalization(dat2)
dat3<-standardNormalization(dat3)

dat_t<-list(t(dat1),t(dat2),t(dat3))
NUM.NEIGHBORS.RATIO = 6 
nemo.affinity.graph <- function(raw.data, k=NA) {
        if (is.na(k)) {
                k = as.numeric(lapply(1:length(raw.data), function(i) round(ncol(raw.data[[i]]) / NUM.NEIGHBORS.RATIO)))
        } else if (length(k) == 1) {
                k = rep(k, length(raw.data))
        }
        sim.data = lapply(1:length(raw.data), function(i) {affinityMatrix(dist2(as.matrix(t(raw.data[[i]])),
                                                                                as.matrix(t(raw.data[[i]]))), k[i], 0.5)})
        affinity.per.omic = lapply(1:length(raw.data), function(i) {
                sim.datum = sim.data[[i]]
                non.sym.knn = apply(sim.datum, 1, function(sim.row) {
                        returned.row = sim.row
                        threshold = sort(sim.row, decreasing = TRUE)[k[i]]
                        returned.row[sim.row < threshold] = 0
                        row.sum = sum(returned.row)
                        returned.row[sim.row >= threshold] = returned.row[sim.row >= threshold] / row.sum
                        return(returned.row)
                })
                sym.knn = non.sym.knn + t(non.sym.knn)
                return(sym.knn)
        })
        patient.names = Reduce(union, lapply(raw.data, colnames))
        num.patients = length(patient.names)
        returned.affinity.matrix = matrix(0, ncol = num.patients, nrow=num.patients)
        rownames(returned.affinity.matrix) = patient.names
        colnames(returned.affinity.matrix) = patient.names
        
        shared.omic.count = matrix(0, ncol = num.patients, nrow=num.patients)
        rownames(shared.omic.count) = patient.names
        colnames(shared.omic.count) = patient.names
        
        for (j in 1:length(raw.data)) {
                curr.omic.patients = colnames(raw.data[[j]])
                returned.affinity.matrix[curr.omic.patients, curr.omic.patients] = returned.affinity.matrix[curr.omic.patients, curr.omic.patients] + affinity.per.omic[[j]][curr.omic.patients, curr.omic.patients]
                shared.omic.count[curr.omic.patients, curr.omic.patients] = shared.omic.count[curr.omic.patients, curr.omic.patients] + 1
        }
        
        final.ret = returned.affinity.matrix / shared.omic.count
        lower.tri.ret = final.ret[lower.tri(final.ret)]
        final.ret[shared.omic.count == 0] = mean(lower.tri.ret[!is.na(lower.tri.ret)])
        
        return(final.ret)
}
W_nemo<- nemo.affinity.graph(dat_t)

# obtain optimal cluster number
estimationResult = estimateNumberOfClustersGivenGraph(W_nemo, 2:8)
k.nemo<-estimationResult$`Eigen-gap best`
cls.nemo = spectralClustering(W_nemo,k.nemo)


