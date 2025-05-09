# GCCMSingle(xEmbedings,yPred,lib_size,pred,totalRow,totalCol,b,cores)
GCCMSingle<-function(xEmbedings,yPred,lib_size,pred,totalRow,totalCol,b,winStepRatio,cores=NULL,dir,validRatio)
{
  x_xmap_y <- data.frame()
  
  
  if(is.null(cores))
  {
    
    for(r in seq(1,(totalRow-lib_size+1),round(1+winStepRatio* lib_size)))
    {
      
      x_xmap_y<-rbind(x_xmap_y,GCCMSingleInner(xEmbedings,yPred,lib_size,pred,totalRow,totalCol,b,r,winStepRatio,dir,validRatio))
    }
    
  }else
  {
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    clusterExport(cl,deparse(substitute(GCCMSingleInner)))
    clusterExport(cl,deparse(substitute(locate)))
    clusterExport(cl,deparse(substitute(projection)))
    clusterExport(cl,deparse(substitute(distance_Com)))
    clusterExport(cl,deparse(substitute(compute_stats)))
    
    x_xmap_y <- foreach(r=seq(1,(totalRow-lib_size+1),round(1+winStepRatio* lib_size)), .combine='rbind') %dopar% GCCMSingleInner(xEmbedings,yPred,lib_size,pred,totalRow,totalCol,b,r,winStepRatio,dir,validRatio)
    stopCluster(cl)
  }

  
  return(x_xmap_y)
}

# GCCMSingleInner(xEmbedings,yPred,lib_size,pred,totalRow,totalCol,b,r)
GCCMSingleInner<-function(xEmbedings,yPred,lib_size,pred,totalRow,totalCol,b,r,winStepRatio,dir,validRatio)
{
  x_xmap_y <- data.frame()
  
  for(c in seq(1,(totalCol-lib_size+1),round(1+winStepRatio* lib_size)))
  {
    
    
    pred_indices <- rep.int(FALSE, times = totalRow*totalCol)
    lib_indices<- rep.int(FALSE, times = totalRow*totalCol)
    
    pred_indices[locate(pred[,1],pred[,2],totalRow,totalCol) ]<-TRUE
    
    pred_indices[which(is.na(yPred)) ]<-FALSE
    
    
    lib_rows<-seq(r,(r+lib_size-1))
    lib_cols<-seq(c,(c+lib_size-1))
    
    lib_ids<-merge(lib_rows,lib_cols)
    
    
    lib_indices[locate(lib_ids[,1],lib_ids[,2],totalRow,totalCol)]<-TRUE
    
    
    if(length(which(is.na(yPred[which(lib_indices)]))) > ((lib_size*lib_size)/2))
    {
      next
    }
    
    # run cross map and store results
    results <-  projection(xEmbedings,yPred,lib_indices ,pred_indices,b,dir,validRatio)
    
    
    x_xmap_y <- rbind(x_xmap_y, data.frame(L = lib_size, rho = results$stats$rho)) 
  }
  
  return(x_xmap_y)
}


# GCCM(xMatrixM, yMatrixM, lib_sizes, lib, pred, E,cores=32)
GCCM<-function(xMatrix, yMatrix, lib_sizes, lib, pred, E, tau = 1, b = E+2, winStepRatio=0,cores=NULL,dir=0,validRatio=0)
{
  
  ### xMatrix: Two dimension matrix of X variable
  ### yMatrix: Two dimension matrix of Y variable
  ### lib_sizes: the sizes of libs to be calculated, they will appears at the horizontal axis
  ### lib: is not used in this version. It is kept for users to define a certain part of the matrix to be used in the cross mapping prediction 
  #         In this version, we use all spatial units of xMatrix and yMatrix by default.
  ### pred: indices of spatial units to be predicted. Due to spatial correlation, the prediction precision of nearby units are  
  #         almost the same, so we can skip some spatial units to  speed up the algorithm
  ### E: is the number of spatial lags used to construct the embeding
  ### b: number of nearest neighbors to use for prediction. We set it default to E+2 according to 
  #               suggestions from Katharina M Bracher，Stuart King and Ian Maccormick of University of Edinburgh, Edinburgh.
  #               Here we acknowledge their contributions 
  ### tau:is the step of spatial lags. If tau=1, the first ring around the focus unit is used as the when E=1; while if tau=2,the 
  #       second ring around the focus unit is used as the when E=1  
  ### winStepRatio: is a speedup parameter. In each prediction, a sliding window with the lib_size is used to confine the number of points in the state space.
  #                 If the matrix is very large, it is time consuming. We can increase the sliding step with winStepRatio. The winStepRatio will be multiplied with 
  #                 the width/height to set the sliding steps  
  #                 
  ### cores: the number of cores can be used for parallel computing
  ### dir: direction parameter for anisotropy. It used to select spatial units from  spatial lags to be used to reconstruct the state sapce. That means only spatial lags 
  ###      in the direction defined by dir take into the reconstruction of the state space. 
  ###      dir=0, all directions are considered. dir=1, Northeast; dir=2, North; dir=3, Northwest; dir=4, West; dir=5, Southwest; dir=6, South; dir=7, Southeast; dir=8, East;  
  ###      dir can also be a vector with more than one directions. For example, dir=c(1,2), then both Northeast and North will be used; dir=c(1,5),Northeast and Southwest
  ### validRatio: is the parameters to handle NA values. When the study area have too many NA (or Nodata) values, we would get NA results. To handle the NA values，we will 
  #               neglect the NA values of target variables. But it will also causes NA results or unstable predictions. So the validRatio is used to expand the neighbors 
  #               to farther sate points: maxDistacne+meanDistance*validRatio，where maxDistacne is the largest distance of the closest b neighbors, the meanDistance is the 
  #               average distance f the closest b neighbors. If validRatio=0.01, then the candidates set were expanded to 0.01*meanDistance. Then after removing NA values of the 
  #               closest b neighbors, farther sate points in candidates set could be used to replenish
  
  
  imageSize<-dim(xMatrix)
  totalRow<-imageSize[1]
  totalCol<-imageSize[2]
  
  yPred<- as.vector(t(yMatrix))
  
  xEmbedings<-list()
  xEmbedings[[1]]<- as.vector(t(xMatrix))
  
  for(i in 1:E)
  {
    xEmbedings[[i+1]]<-laggedVariableAs2Dim(xMatrix, i*tau)  #### row first
  
  }
  
  x_xmap_y <- data.frame()

  for(lib_size in lib_sizes)
  {
   
     x_xmap_y<-rbind(x_xmap_y,GCCMSingle(xEmbedings,yPred,lib_size,pred,totalRow,totalCol,b,winStepRatio,cores,dir,validRatio))
      
  }
  
  return (x_xmap_y)
}

locate<-function(curRow,curCOl,totalRow,totalCol)
{
  return ((curRow-1)*totalCol+curCOl) 
}

# projection(xEmbedings,yPred,lib_indices ,pred_indices,b)
projection<-function(embedings,target,lib_indices, pred_indices,num_neighbors,dir,validRatio=0)
{
  makepred <- rep.int(NaN, times = length(target))
  
  for(p in which (pred_indices))
  {
    temp_lib <- lib_indices[p]
    lib_indices[p] <- FALSE
    
    libs <- which(lib_indices)
    
    distances<-distance_Com(embedings,libs,p,dir)
    #distances<-colMeans(distances)
    distances_sorted<-sort(distances,na.last = TRUE)
    nearest_distances<-distances_sorted[1:num_neighbors]
    meandis<-mean(nearest_distances,na.rm=TRUE)
    maxdis<-max(nearest_distances,na.rm = TRUE)
    distances_Ex<-maxdis+validRatio*meandis
    distances[distances>distances_Ex]<-NA
    
    
    
    ###########deal with Na###########
    target4Lib<-target[libs]
    distances[which(is.na(target4Lib))]<-NA
    ####################
    # find nearest neighbors
    # neighbors <- order(distances)[1:num_neighbors]
    neighbors <- order(distances,na.last = TRUE)[1:num_neighbors]
    min_distance <- distances[neighbors[1]]
    if(is.na(min_distance))
    {
      
      lib_indices[p] <- temp_lib 
      
      next
    }
    # compute weights
    if(min_distance == 0) # perfect match
    {
      weights <- rep.int(0.000001, times = num_neighbors)
      weights[distances[neighbors] == 0] <- 1
      flagstore<-c()
      for (flag in 1:num_neighbors)
      {
        if(is.na(target[libs[neighbors]][flag]))# if(is.na(test[flag]))
        {
          next
        }
        flagstore<-c(flagstore,flag)
        
      }
      numerator=0
      denominator=0
      if(!(is.null(flagstore)))
      {
        for (calcuflag in flagstore)
        {
          numerator<-numerator+(weights[calcuflag]*target[libs[neighbors]][calcuflag])
          
          denominator<-denominator+weights[calcuflag]
          
        }
      }
    }else{
     
      if(length(which(is.na(distances[neighbors])))!=0)
      {
        
        weights <- exp(-distances[neighbors]/min_distance)
        weights[weights < 0.000001] <- 0.000001
        
        Truenieghbors<-neighbors[which(!is.na(weights))]
      
       
        flagstore<-c()
        for (flag in 1:length(Truenieghbors))
        {
          
          if(is.na(target[libs[Truenieghbors]][flag]))# if(is.na(test[flag]))
          {
            next
          }
          flagstore<-c(flagstore,flag)
          
        }
        
        numerator=0
        denominator=0
       
        if(!(is.null(flagstore)))
        {
          for (calcuflag in flagstore) 
          {
            numerator<-numerator+(weights[calcuflag]*target[libs[neighbors]][calcuflag])
            
            denominator<-denominator+weights[calcuflag]
            
           
            
          }
        }
      }else{
       
        weights <- exp(-distances[neighbors]/min_distance)
        weights[weights < 0.000001] <- 0.000001
        flagstore<-c()
       
        for (flag in 1:num_neighbors)
        {
          if(is.na(target[libs[neighbors]][flag]))# if(is.na(test[flag]))
          {
            next
          }
          flagstore<-c(flagstore,flag)
          
        }
        numerator=0
        denominator=0
        if(!(is.null(flagstore)))
        {
          for (calcuflag in flagstore) 
          {
            numerator<-numerator+(weights[calcuflag]*target[libs[neighbors]][calcuflag])
           
            denominator<-denominator+weights[calcuflag]
           
          } 
        }
      }
      
    }
    
    
    
    makepred[p] <-numerator/denominator
    
    lib_indices[p] <- temp_lib 
    numerator=NaN
    denominator=NaN
  }
  
  # return output & stats
  return(list(pred = makepred, stats = compute_stats(target[pred_indices], makepred[pred_indices])))
  
  
}


# distance_Com(embedings,libs,p)
distance_Com<-function(embeddings,libs,p,dir=0)
{
  distances<-c()
  
  emd<-embeddings[[1]]
  
  distances<-cbind(distances,abs(emd[libs]-emd[p]))
  
  
  
  for(e in 2:length(embeddings))
  {
    emd<-embeddings[[e]]
    
    dirIndices<-seq(1:dim(emd)[2])
    
    if(dir[1]!=0)
    {
      allDirIndices<-seq(1,by=(e-1),length.out=8)
      dirIndices<-allDirIndices[dir]
    }
    
    
    #q <- matrix(rep(emd[p], length(libs)), nrow = length(libs), byrow = T)
    p_matrix<-matrix(emd[p,dirIndices], nrow = length(libs), ncol = length(dirIndices), byrow = TRUE )
    distances<-cbind(distances,abs(emd[libs,dirIndices]-p_matrix))
    
    # distances<-cbind(distances,abs(emd[libs,]-emd[p,]))
    
  }
  
  return (rowMeans(distances,na.rm=TRUE))
  
}

