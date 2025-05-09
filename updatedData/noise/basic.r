
confidence<- function(r,n,level=0.05)
{
  z<-1/2*log((1+r)/(1-r))
  ztheta<-1/sqrt(n-3)
  
  qZ<-qnorm(1-level/2)
  
  upper<-z+qZ*ztheta
  
  lower<-z-qZ*ztheta
  
  r_upper<-(exp(2*upper)-1)/(exp(2*upper)+1)

  r_lower<-(exp(2*lower)-1)/(exp(2*lower)+1)
  
  
  return (cbind( r_upper, r_lower))
  
}



significance<- function(r,n)
{
  t<-r*sqrt((n-2)/(1-r*r))
 
  return (1-pt(t,n-2))*2
  
}

expandMatrix<- function(dataMatrix,lagNum)
{
  
  if(lagNum<0)
  {
    return (dataMatrix)
  }
  
  lagNum<-round(lagNum)
  
  if(lagNum>1)
  {
    dataMatrix<-expandMatrix(dataMatrix,lagNum-1)
  }
  
  ColNum<-ncol(dataMatrix)
  RowNum<-nrow(dataMatrix)
  
  dataMatrix<-rbind(rep(NA,ColNum),dataMatrix)
  dataMatrix<-rbind(dataMatrix,rep(NA,ColNum))  
  dataMatrix<-cbind(rep(NA,ColNum)+2,dataMatrix) 
  dataMatrix<-cbind(dataMatrix,rep(NA,ColNum)+2) 
  
  return(dataMatrix)
  
}


laggedVariable<-function(dataMatrix,lagNum)  
{
  ColNum<-ncol(dataMatrix)
  RowNum<-nrow(dataMatrix)
  exDataMatrix<-expandMatrix(dataMatrix,lagNum)
  
  laggedVar<-array(rep(NA,ColNum*RowNum*8*lagNum),dim=c(RowNum*ColNum,8*lagNum))
  
  for(r in 1:RowNum)
  {
    for(c in 1:ColNum)
    {
      item<-1
      exr<-r+lagNum
      exc<-c+lagNum
      #############Start from Notheast, fisr North################
      for(la in seq(lagNum,(1-lagNum)))
      {
        #####Row is  row-lagNum
        laggedVar[r,c,item]<- exDataMatrix[exr-lagNum,exc+la]
        item<-item+1
      }
      
      #############Then West################
      for(ra in seq(-lagNum,(lagNum-1)))
      {
        #####Row is  row-lagNum
        laggedVar[r,c,item]<- exDataMatrix[exr+ra,exc-lagNum]
        item<-item+1
      }
      
      #############Then South################
      
      for(la in seq(-lagNum,(lagNum-1)))
      {
        #####Row is  row-lagNum
        laggedVar[r,c,item]<- exDataMatrix[exr+lagNum,exc+la]
        item<-item+1
      }
      
      
      #############Then East################
      for(ra in seq(lagNum,(1-lagNum)))
      {
        #####Row is  row-lagNum
        laggedVar[r,c,item]<- exDataMatrix[exr+ra,exc+lagNum]
        item<-item+1
      }
      
    }
  }
  return (laggedVar)
}


laggedVariableAs2Dim<-function(dataMatrix,lagNum)  
{
  ColNum<-ncol(dataMatrix)
  RowNum<-nrow(dataMatrix)
  exDataMatrix<-expandMatrix(dataMatrix,lagNum)
  
  laggedVar<-array(rep(NA,ColNum*RowNum*8*lagNum),dim=c(RowNum*ColNum,8*lagNum))
  
  ###By row; from top to down
  for(r in 1:RowNum)
  {
    for(c in 1:ColNum)
    {
      item<-1
      exr<-r+lagNum
      exc<-c+lagNum
      #############Start from Notheast, fisr North################
      for(la in seq(lagNum,(1-lagNum)))
      {
        #####Row is  row-lagNum
        laggedVar[(r-1)*ColNum+c,item]<- exDataMatrix[exr-lagNum,exc+la]
        item<-item+1
      }
      
      #############Then West################
      for(ra in seq(-lagNum,(lagNum-1)))
      {
        #####Row is  row-lagNum
        laggedVar[(r-1)*ColNum+c,item]<- exDataMatrix[exr+ra,exc-lagNum]
        item<-item+1
      }
      
      #############Then South################
      
      for(la in seq(-lagNum,(lagNum-1)))
      {
        #####Row is  row-lagNum
        laggedVar[(r-1)*ColNum+c,item]<- exDataMatrix[exr+lagNum,exc+la]
        item<-item+1
      }
      
      
      #############Then East################
      for(ra in seq(lagNum,(1-lagNum)))
      {
        #####Row is  row-lagNum
        laggedVar[(r-1)*ColNum+c,item]<- exDataMatrix[exr+ra,exc+lagNum]
        item<-item+1
      }
      
    }
  }
  return (laggedVar)
}



distance<-function(emd1,emd2)
{
  if(is.factor(emd1))
  {
    return (distanceCat(emd1,emd2))
  }else if(is.numeric(emd1))
  {
    return (distanceCon(emd1,emd2))
  }else
  {
    return (NA)
  }
  
}

 
  
distanceCon<-function(emd1,emd2)
{
  naIn1<-which(is.na(emd1))
  naIn2<-which(is.na(emd2))
  
  emd2[naIn1]<-NA
  emd1[naIn2]<-NA
  return(sum(abs(emd1-emd2),na.rm = TRUE)) 
  
}


distanceCat<-function(emd1,emd2)
{
  
  naIn1<-which(is.na(emd1))
  naIn2<-which(is.na(emd2))
  
  emd2[naIn1]<-NA
  emd1[naIn2]<-NA
  
  table1<-table(emd1)
  table2<-table(emd2)
  
  fractions1<-as.dataframe(table1/sum(table1))
  
  fractions2<-as.dataframe(table2/sum(table2))
  
  merged<-merge(fractions1,fractions2,by.x=1,by.y=1,all=TRUE)
  merged[is.na(merged)]<-0
  
  return (sum(abs(merged[,2]-merged[,3]),na.rm = TRUE))
  
}





distanceCross<-function(trainData,embeddings,neigborNum,targetVar,xName="x",yName="y")
{
  sampleNum<-length(targetVar);
  
  values<-c()
  geoDis<-c()
  
  reOrders<-order(trainData[,targetVar])
  trainData<-trainData[reOrders,]
  
  N<-nrow(trainData) 
  
  
  for(i in  1:(N-1))
  {
    for(j in (i+1):(min(i+neigborNum, N)))
    {
      values<-c(values,abs(trainData[i,targetVar]-trainData[j,targetVar]))
      geoDis<-c(geoDis,sqrt((trainData[i,xName]-trainData[j,xName])^2+(trainData[i,yName]-trainData[j,yName])^2))
    }
      
  }
  
  resulsDistance<-data.frame(values)
  
  for(e in 1:length(embeddings) )
  {
    emd<-embeddings[[e]]
    distances<-c()
    if(is.null(dim(emd)))
    {
      emd<-emd[reOrders]
      
      N<-length(emd)
      
      for(i in  1:(N-1))
      {
        for(j in (i+1): (min(i+neigborNum, N)))
        {
          distances<-c(distances,abs(emd[i]- emd[j]))
        }
      }
      
    }else
    {
      emd<-emd[reOrders,]
      N<-nrow(emd)
      
      
      for(i in  1:(N-1))
      {
        for(j in (i+1):(min(i+neigborNum, N)))
        {
          
          distances<-c(distances,distance(emd[i,],emd[j,]))
        }
      }
    }
    
    resulsDistance<-cbind(resulsDistance,distances)
    
  }
  
  resulsDistance<-cbind(resulsDistance,geoDis)
  return (resulsDistance)
}

distanceCrossGeo<-function(sample,xName="x",yName="y")
{
  sampleNum<-nrow(sample)
  
  
  
  values<-c()
  
  for(i in  1:(sampleNum-1))
  {
    for(j in (i+1):sampleNum)
    {
      
      values<-c(values,sqrt((sample[i,xName]-sample[j,xName])^2+(sample[i,yName]-sample[j,yName])^2))
    }
  }
  
  return (values)
}

distanceBtwGeo<-function(test,train, xName="x",yName="y")
{
  
  values<-c()
  
  for(i in 1:nrow(train))
  {
    for(j in 1:nrow(test))
    {
      values<-c(values,sqrt((train[i,xName]-test[j,xName])^2+(train[i,yName]-test[j,yName])^2))
    }
  }
  
  return (values)
}


distanceBtw<-function(predEmd,embeddings)
{
  resulsDistance<-c()
  for(e in 1:length(embeddings) )
  {
    emd<-embeddings[[e]]
    distances<-c()
    if(is.null(dim(emd)))
    {
      for(i in  1:(length(emd)))
      {
        distances<-c(distances,abs(emd[i]- predEmd[[e]]))
      }
      
    }else
    {
      for(i in  1:(nrow(emd)))
      {
        distances<-c(distances, distance(emd[i,],predEmd[[e]]))
      }
      
    }
    resulsDistance<-cbind(resulsDistance,distances)
    
  }
  return (resulsDistance)
}


getGridEmbedings<-function(data,nrow,ncol,byrow,atts,lags)
{
  embeddings<-list()
  
  num<-0
  for(a in 1:length(atts))
  {
    for(l in lags)  
    {
      
      num<-num+1
      if(l==0)
      {
        embeddings[[num]]<-data[,atts[a]]
      }else
      {
        embeddings[[num]]<-laggedVariableAs2Dim(matrix(data[,atts[a]],nrow=nrow,ncol=ncol,byrow=byrow),l)
      
      }
      
    }
  }
  return (embeddings)
}


getSampleEmbedings<-function(ids,embeddings)
{
  sampleEmbedings<-list()
  for(e in 1:length(embeddings))
  {
    tmpEmbedings<-embeddings[[e]]
    
    if(is.null(dim(tmpEmbedings)))
    {
      sampleEmbedings[[e]]<-tmpEmbedings[ids]
    }else
    {
      sampleEmbedings[[e]]<-tmpEmbedings[ids,]
    }
  }
  
  return (sampleEmbedings)
}

normalize <- function(x) 
{
  return((x - min(x)) / (max(x) - min(x)))
}



dominate<-function(solutions)
{
  paretos<-solutions[1,]
  paretoIndex<-c(1)
  
  for(i in 2:nrow(solutions))
  {
    removes<-c()
    added<-TRUE
    s<-solutions[i,]
    for(j in 1:nrow(paretos))
    {
      p<-paretos[j,]
      if(all((p-s)>0))                          #all()   any()
      {
        removes<-c(removes,j)
        
      }else if(any((p-s)>0))
      {
        next
      }else
      {
        added<-FALSE
        break
      }
    }
    if(length(removes)>0)
    {
      paretos<-paretos[-removes,]
      paretoIndex<-paretoIndex[-removes]
    }
    if(added)
    {
      paretos<-rbind(paretos, s)
      paretoIndex<-c(paretoIndex,i)
    }
  }
  return (paretoIndex)
}

distanceCrossSign<-function(data,targetVars)
{
  results<-c()
  n<-nrow(data)

  for(targetV in targetVars)
  {
    if(is.numeric(data[,targetV]))
    {
      values<-c()
      for(i in  1:n)
      {
        for(j in 1:n)
        {
          values<-c(values, (data[i,targetV]- data[j,targetV]))
        }
      }
      results<-cbind(results,values)
    }else
    {
      values1<-c()
      values2<-c()
      for(i in  1:n)
      {
        for(j in 1:n)
        {
          values1<-c(values1,data[i,targetV])
          values2<-c(values2,data[j,targetV])
        }
      }
      results<-cbind(results,values1,values2)
    }
    
  }
  
  results<-as.data.frame(results)
  
  names(results)<-targetVars
  
  return (results)
  
  #  distance<-function(emd1,emd2)
  #  {
  #    if(is.factor(emd1))
  #    {
  #      return (distanceCat(emd1,emd2))
  #    }else if(is.numeric(emd1))
  #    {
  #      return (distanceCon(emd1,emd2))
  #    }else
  #    {
  #      return (NA)
  #    }
  
  #  }
  
}


distanceSign<-function(former,latter,targetVars)
{
  results<-c()
  
 for(f in 1:nrow(former) )
 {
   for(l in 1:nrow(latter) )
   {
     values<-c()
     for(targetV in targetVars)
     {
       if(is.numeric(former[f,targetV]))
       {
         values<-c(values,former[f,targetV]-latter[l,targetV])
       }else
       {
         values<-c(values,former[f,targetV],latter[l,targetV])
       }
     }
     
     results<-rbind(results,values)
   }
 }
  
  results<-as.data.frame(results)
  
  names(results)<-targetVars
  
  return (results)
  
}


compute_stats <- function(obs, pred)
{
  # computes performance metrics for how well predictions match observations
  # obs = vector of observations
  # pred = vector of prediction
  
  N = sum(is.finite(obs) & is.finite(pred))
  rho = cor(obs, pred, use = "pairwise.complete.obs")
  mae = mean(abs(obs-pred), na.rm = TRUE)
  rmse = sqrt(mean((obs-pred)^2, na.rm = TRUE))
  return(data.frame(N = N, rho = rho, mae = mae, rmse = rmse))
}

