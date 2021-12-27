getwd()
setwd('C:/Users/asus_pinar/desktop/files')
getwd()

#library(genlasso)
library(data.table)
library(rpart)
library(repr)
require(scatterplot3d)
require(FNN)
library(ggplot2)
require(randomForest)
require(TunePareto)
options(scipen=999)
options(repr.plot.width=15, repr.plot.height=8)

#function definitions

nn_classify_cv=function(dist_matrix,train_class,test_indices,k){
  
  test_distances_to_train=as.matrix(dist_matrix[test_indices,])
  test_distances_to_train=test_distances_to_train[,-test_indices]
  trainclass=train_class[-test_indices]
  #print(str(test_distances_to_train))
  ordered_indices=apply(test_distances_to_train,1,order)
  if(k==1){
    nearest_class=(trainclass[as.numeric(ordered_indices[1,])])
    nearest_class=data.table(id=test_indices,nearest_class)
  } else {
    nearest_class=apply(ordered_indices[1:k,],2,function(x) {trainclass[x]})
    nearest_class=as.data.table(nearest_class)
    nearest_class=data.table(id=test_indices,t(nearest_class))
  }
  
  long_nn_class=melt(nearest_class,'id')
  
  class_counts=long_nn_class[,.N,list(id,value)]
  class_counts[,predicted_prob:=N/k]
  wide_class_prob_predictions=dcast(class_counts,id~value,value.var='predicted_prob')
  wide_class_prob_predictions[is.na(wide_class_prob_predictions)]=0
  class_predictions=class_counts[,list(predicted=value[which.max(N)]),by=list(id)]
  
  
  return(list(prediction=class_predictions,prob_estimates=wide_class_prob_predictions))
  
}


RF_classify_cv=function(data_df, test_indices, ntr){
                        
  data_df = long_repr          
  test_df=data_df[test_indices,]
  train_df=data_df[-test_indices,]
  
  test_features = test_df[,-c("BagID", "BagLabel")]
  test_target = test_df[,"BagLabel"]
  test_df = test_df[,-("BagID")]
  
  train_features = train_df[,-c("BagID", "BagLabel")]
  train_target = train_df[,"BagLabel"]
  train_df = train_df[,-("BagID")]
  
  
  fit_to_repr_rf=randomForest(train_df[,-c('BagLabel'),with=F],
                         as.factor(train_df$BagLabel),
                         importance=TRUE,
                         ntree=ntr) 

  predicted=predict(fit_to_repr_rf,test_df)

  return(predicted)
  
}

musk_data = read.csv("Musk1.csv")
musk_data = as.data.table(musk_data)

colnames(musk_data)[1]=c("BagLabel")
colnames(musk_data)[2]=c("BagID")

head(musk_data)
tail(musk_data)
summary(musk_data)

c(unique(musk_data[,"BagID"]))

print("1 Labeled Rows Count ")
nrow(musk_data[BagLabel==1,])

print("0 Labeled Rows Count ")
nrow(musk_data[BagLabel==0,])

print("1 Labeled Bags Count ")
sum(musk_data[, list(Labs=mean(BagLabel)), BagID]$Labs==1)

print("0 Labeled Rows Count ")
sum(musk_data[, list(Labs=mean(BagLabel)), BagID]$Labs==0)

print("NA variables")
sum(is.na(musk_data))


#create bagged df
bagged_df = data.table(unique(musk_data[,"BagID"])) #bagID
bagged_df[,nof_instances:=as.numeric(table(musk_data$BagID))]   #number of instances in each bag
bagged_df[,BagLabel:=musk_data[, list(Labs=mean(BagLabel)), BagID]$Labs]  #labels of each bag
summary(bagged_df)

#The first column is a binary column that shows if the bag is musk or non-musk. Instances are already clustered into bags. All instances in a bag assumed to have the same label as the bag. Second column shows the bag ID's. There are 92 bags of which 47 are labeled as 1 and 45 as 0. Number of observations in each bag differs with a minimum of 2 and a maximum of 40.

#Bags represent 92 molecules which are classified as musk or non-musk by human experts. Due to the rotation of bonds a molecule could have different shapes, conformations. 475 observations shows different conformations of molecules represented with the features in remaining 166 columns. All features are numeric with different ranges

#The goal is to learn to predict whether a new molecule is a musk or a non-musk.


## EXPLORING FEATURES

#Renaming features X1 to X166
feature_list=c()
for( i in 1:166){
  feature_list = append(feature_list,(paste("X", i, sep="")))
}
colnames(musk_data)[3:168] = feature_list
#head(musk_data,5)

musk_data[,index:=1:.N]

## PLOTTING FEATURES
musk_melted=(melt(musk_data[,-c(2)] ,  id.vars = c('index','BagLabel'), variable.name = 'feature'))
musk_melted=musk_melted[order(musk_melted$index)]

for( ft in colnames(musk_data)){
  print(ggplot(data = musk_melted[feature==ft])+geom_point(aes(y=value,
                                                               x=index,
                                                               col=as.factor(musk_melted[feature==ft]$BagLabel))))
}
########################

#DISTANCE VISUALIZATION WITH MDS
#for instances (euclidean)

normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}

musk_feature_df=musk_data[,-c(1,2,169)]
musk_normalized=as.data.frame(apply(musk_feature_df, FUN=normalize,MARGIN = 2))


#get euclidean distances (instance level)
instance_distances = dist(x = musk_normalized,method = "euclidean",upper = TRUE)
instance_distances=as.matrix(instance_distances)

#apply mds (to 2d)
mds_musk=cmdscale(instance_distances, k=2, eig=T)
mds_musk_df = as.data.frame(mds_musk$points)
mds_musk_df = cbind(mds_musk_df, index=musk_data$index)
mds_musk_df = cbind(mds_musk_df, labels=musk_data$BagLabel)
#mds_musk$eig


(ggplot(mapping = aes(x=V1 , y=V2),data=(mds_musk_df))+geom_point(aes(col=as.factor(labels)) )
  +labs(title="MDS 2D Musk Data")
  +theme_minimal())


#apply mds (to 3d)
mds_musk_3d=cmdscale(instance_distances, k=3, eig=T)
mds_musk3d_df = as.data.frame(mds_musk_3d$points)
mds_musk3d_df = cbind(mds_musk3d_df, index=musk_data$index)
mds_musk3d_df = cbind(mds_musk3d_df, labels=musk_data$BagLabel)

scatterplot3d(x=mds_musk3d_df$V1, y=mds_musk3d_df$V2, z=mds_musk3d_df$V3,main="MDS 3D Musk Data",
              xlab = "X", ylab= "Y", zlab="Z", color=as.factor((mds_musk3d_df$labels)+3))


#pca
pca_musk = prcomp(musk_normalized)
summary(pca_musk)

#idea 1 = referenced from paper k-means clustering and represent with instances in each cluster but normalize
#idea 2 = represent with mean distance of instances to each cluster

#classification approach 1 : NN classification
#classification approach 2 : random forest


musk_normalized=as.data.frame(apply(musk_feature_df, FUN=normalize,MARGIN = 2))

############## KMEANS

#iter.max is set to 25 in the original application and therefore it is also set to 25 in this work
#nstarts is set to 5 in the original application tuned with a list of (4,5,6,7) in this work



######## 

#idea 1 = referenced from paper k-means clustering and represent with instances in each cluster but normalize

n_fold=10
nof_rep=1

nofcenters_set = c(16,20,24,30,34,40) #even numbers
nstarts_set = c(5,10)

k_set=c(1,3,5,7,10)
ntree_set=c(500,700,800)

result_nn=vector('list',length(nstarts_set)*length(nofcenters_set)*nof_rep*n_fold*length(k_set))
iter_nn=1

result_rf=vector('list',length(nstarts_set)*length(nofcenters_set)*nof_rep*n_fold*length(ntree_set))
iter_rf=1

for(center in nofcenters_set){
  for (nst in nstarts_set){
    
    #center=15
    #nst=5
    
    #kmeans
    set.seed(1234)
    km=kmeans(musk_normalized,centers=center, nstart = nst, iter.max = 25)
    
    #find nearest clusters for each instance
    nearestCluster=get.knnx(km$centers,musk_normalized,k=1)$nn.index
    
    #obtain new representation for bags: number of instances in each cluster
    instancelevelClusterAssignments=data.table(BagID=musk_data$BagID,Cluster=array(km$cluster))
    clusterrepr=table(instancelevelClusterAssignments)
    
    long_repr=data.table(clusterrepr)
    #long_repr=long_repr[order(long_repr[,"BagID"])]
    
    long_repr=dcast(long_repr,BagID~as.numeric(Cluster), fun.aggregate=sum, value.var ='N')
    long_repr$BagID = as.numeric(long_repr$BagID)
    long_repr=long_repr[order(long_repr$BagID)]
    
    k=long_repr[,-c("BagID")] 
    k= k/bagged_df$nof_instances
    k[,BagID := long_repr$BagID]
    long_repr=k
    
    long_repr = cbind(long_repr,BagLabel=bagged_df$BagLabel)

    #obtain cv folds 10-fold cv stratified
    labs = long_repr$BagLabel
    n_fold=10
    nof_rep=1
    set.seed(2301)
    cv_indices=generateCVRuns(labs, 
                              ntimes =nof_rep, 
                              nfold = n_fold,
                              leaveOneOut = FALSE, 
                              stratified = TRUE)
    
    #approach 1: nearest neighbor classification 
    
    #get distance matrix with new representation
    features_dummy = long_repr[,-c("BagID", "BagLabel") ]
    repr_dist = dist(features_dummy,method = "euclidean",upper = TRUE)
    repr_dist = as.matrix(repr_dist)
    diag(repr_dist)=999999

    
    for(i in 1:nof_rep){
        this_fold=cv_indices[[i]]
        for(j in 1:n_fold){
          test_indices=this_fold[[j]]
          for(k in 1:length(k_set)){
            
            current_k=k_set[k]
            current_fold=nn_classify_cv(repr_dist,labs,test_indices,k=current_k)
            
            accuracy=sum(labs[test_indices]==current_fold$prediction$predicted)/length(test_indices)
            tmp=data.table(kmeanscenters=center, startsets=nst ,repid=i,foldid=j,
                           k=current_k,acc=accuracy)
            result_nn[[iter_nn]]=tmp
            iter_nn=iter_nn+1
            
        }
      }
    }
    
    ######
    #approach 2 -Random Forest for Classification

    for(i in 1:nof_rep){
      this_fold=cv_indices[[i]]
      for(j in 1:n_fold){
        test_indices=this_fold[[j]]
        for(noftree in 1:length(ntree_set)){
          
          current_ntree = ntree_set[noftree]

          current_fold= RF_classify_cv(data_df = long_rep, 
                                       test_indices = test_indices,
                                       ntr = current_ntree)
          
          accuracy=sum(labs[test_indices]==current_fold)/length(test_indices)
          tmp=data.table(kmeanscenters=center, startsets=nst ,repid=i,foldid=j,
                         ntree=current_ntree,acc=accuracy)
          result_rf[[iter_rf]]=tmp
          iter_rf=iter_rf+1

        }
      }
    }

    ################

  }
}

NN_results_1 = rbindlist(result_nn)
NN_summary_1 = NN_results_1[,list(avg_acc=mean(acc),sdev_acc=sd(acc),result_count=.N),by=list(startsets,kmeanscenters,k)]

#order by avg accuracy
NN_summary_1=NN_summary_1[order(NN_summary_1$avg_acc, decreasing = TRUE)]
head(NN_summary_1,25)

#Best average accuracy 0.87 is obtained with 24-means centers and 1-NN classifier. Number of starting sets tried did not affect the performance. Trying 5 sets in start is preferred for simplicity.
#The standard deviations are high, considering the complexity and application time another set of parameters could be preferred. 

#During the trials it is observed that initial cluster centers (since not specified determined according to the seed) strongly affects the performance of the method.


RF_results_1 = rbindlist(result_rf)
RF_summary_1 = RF_results_1[,list(avg_acc=mean(acc),sdev_acc=sd(acc),result_count=.N),by=list(startsets,kmeanscenters,ntree)]

#order by avg accuracy
RF_summary_1=RF_summary_1[order(RF_summary_1$avg_acc, decreasing = TRUE)]
head(RF_summary_1,10)

#Best average accuracy 0.88 is obtained with 30-means centers and 500 trees. Although, average accuracy values with 500, 700 and 800 trees are same 500 trees are preferred to have a simpler model. Trying 5 sets in start is preferred for simplicity.


#########################

#idea 2 = represent with mean distance of instances to each cluster

#sets same as the first representation
nofcenters_set = c(16,20,24,30,34,40) #even numbers
nstarts_set = c(5,10)

k_set=c(1,3,5,7,10)
ntree_set=c(500,700,800)

result_nn=vector('list',length(nstarts_set)*length(nofcenters_set)*nof_rep*n_fold*length(k_set))
iter_nn=1

result_rf=vector('list',length(nstarts_set)*length(nofcenters_set)*nof_rep*n_fold*length(ntree_set))
iter_rf=1

for(center in nofcenters_set){
  for (nst in nstarts_set){
    
    center=15
    nst=5
    
    #kmeans
    set.seed(1234)
    km=kmeans(musk_normalized,centers=center, nstart = nst, iter.max = 25)
    
    #find nearest clusters for each instance
    nearestCluster=get.knnx(km$centers,musk_normalized,k=1)$nn.index
    
    #obtain new representation for bags:
    #for each instance calculate average distance to each cluster: mean of distances to observations in that cluster
    nearestCluster = data.table(nearestCluster)
    colnames(nearestCluster) = "clusterID"
    nearestCluster[, obs_index := 1:.N]

    newdist = dist(musk_feature_df, method = "euclidean",upper = TRUE)
    newdist = as.matrix(newdist)
    
    newReprinstance = data.table()

    
    for(c in 1:center){
      
      clusterindices = nearestCluster[clusterID == c]$obs_index
      
      #get distance with all observations in that cluster
      dummy_dist = newdist[,clusterindices]
      
      #get row means
      newReprinstance = cbind(newReprinstance,rowMeans(dummy_dist))
      
    }
    
    #renaming columns 
    cnames_list = ((function (x) paste("cluster",x,sep="")) (1:center))
    colnames(newReprinstance) = cnames_list
    
    newReprinstance[,obs_index := musk_data$index]
    newReprinstance[,BagID:=musk_data$BagID]
    newReprinstance[,BagLabel:=musk_data$BagLabel]
    
    #get repr of bags: mean distances of instances to clusters

    newReprbag = data.table()
    
    for(b in 1:max(unique(musk_data[,"BagID"]))){
      dmmy_df = rbind(colMeans(newReprinstance[BagID==b,-c("obs_index")]))
      newReprbag = rbind(newReprbag,dmmy_df)
    }
 
    #obtain cv folds 10-fold cv stratified
    labs = newReprbag$BagLabel
    n_fold=10
    nof_rep=1
    set.seed(2301)
    cv_indices=generateCVRuns(labs, 
                              ntimes =nof_rep, 
                              nfold = n_fold,
                              leaveOneOut = FALSE, 
                              stratified = TRUE)
    
    #approach 1: nearest neighbor classification 
    
    #get distance matrix with new representation
    features_dummy = newReprbag[,-c("BagID", "BagLabel")]
    repr_dist = dist(features_dummy,method = "euclidean",upper = TRUE)
    repr_dist = as.matrix(repr_dist)
    diag(repr_dist)=999999
    
    
    for(i in 1:nof_rep){
      this_fold=cv_indices[[i]]
      for(j in 1:n_fold){
        test_indices=this_fold[[j]]
        for(k in 1:length(k_set)){
          
          current_k=k_set[k]
          current_fold=nn_classify_cv(repr_dist,labs,test_indices,k=current_k)
          
          accuracy=sum(labs[test_indices]==current_fold$prediction$predicted)/length(test_indices)
          tmp=data.table(kmeanscenters=center, startsets=nst ,repid=i,foldid=j,
                         k=current_k,acc=accuracy)
          result_nn[[iter_nn]]=tmp
          iter_nn=iter_nn+1
          
        }
      }
    }
    
    ######
    #approach 2 -Random Forest for Classification
    
    for(i in 1:nof_rep){
      this_fold=cv_indices[[i]]
      for(j in 1:n_fold){
        test_indices=this_fold[[j]]
        for(noftree in 1:length(ntree_set)){
          
          current_ntree = ntree_set[noftree]
          
          current_fold= RF_classify_cv(data_df = newReprbag, 
                                       test_indices = test_indices,
                                       ntr = current_ntree)
          
          accuracy=sum(labs[test_indices]==current_fold)/length(test_indices)
          tmp=data.table(kmeanscenters=center, startsets=nst ,repid=i,foldid=j,
                         ntree=current_ntree,acc=accuracy)
          result_rf[[iter_rf]]=tmp
          iter_rf=iter_rf+1
          
        }
      }
    }
    
    ################
    
  }
}


NN_results_2 = rbindlist(result_nn)
NN_summary_2 = NN_results_2[,list(avg_acc=mean(acc),sdev_acc=sd(acc),result_count=.N),by=list(startsets,kmeanscenters,k)]

#order by avg accuracy
NN_summary_2=NN_summary_2[order(NN_summary_2$avg_acc, decreasing = TRUE)]
head(NN_summary_2)

#Best average accuracy 0.88 is obtained with 30-means centers and 1-NN classifier. 
#Considering the high standard deviations algorithm complexity could be taken into account therefore 16-means with 1-NN classifier could be used.

RF_results_2 = rbindlist(result_rf)
RF_summary_2 = RF_results_2[,list(avg_acc=mean(acc),sdev_acc=sd(acc),result_count=.N),by=list(startsets,kmeanscenters,ntree)]

#order by avg accuracy
RF_summary_2=RF_summary_2[order(RF_summary_2$avg_acc, decreasing = TRUE)]
head(RF_summary_2,10)

#The algorithm worked with 0.75 accuracy independent of the parameters. 

#ADD REFERENCES

