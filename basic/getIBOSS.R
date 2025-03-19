if(TRUE){
  samp.data=function(data,beta,num.covariate,num.datalines,num.samp,delta)
  {
    dim = num.covariate
    x.matrix=data$x.matrix
    y=data$y
    samp.y=c()
    samp.x=data.frame()
    dim.samp=ceiling(num.samp/(2*(dim-1)))
    op.result=optimize(min.fun,num.covariate=num.covariate,interval=c(-10,10))
    c=op.result$minimum
    
    
    c_index=as.matrix(x.matrix)%*%beta
    distance1=abs(c_index-c)
    distance2=abs(c_index+c)
    index=(distance1>distance2)
    distance1[index]=0
    distance2[!index]=0
    distance=distance1+distance2
    
    index.update=(distance<delta)
    distance=distance[index.update]
    c_index=c_index[index.update];
    
    samp.x=x.matrix[index.update,]
    
    samp.y=y[index.update]
    x.matrix=samp.x
    y=samp.y
    
    u=length(y)
    
    samp.y=c()
    samp.x=c()
    
    use_column <- x.matrix[,2]
    x.matrix_temp <- x.matrix
    
    for(i in 2:(dim))
    {
      index <- order(use_column)[c(1:dim.samp,(length(use_column)-dim.samp+1):length(use_column))]
      
      samp.x<-rbind(samp.x,x.matrix[index,])
      samp.y<-c(samp.y, y[index])
      y<- y[-index]
      if (i < dim) {
        use_column <- x.matrix[-index,i+1]
      }
      x.matrix <- x.matrix[-index,]
    }
    
    return(list(y=samp.y,x.matrix=samp.x))
  }
}else{
  src <- '
#include <algorithm>
#include <iostream>
Rcpp::NumericVector x(Rx);
Rcpp::NumericVector f(Rf);
Rcpp::NumericVector r(Rr);//subsample size
Rcpp::NumericVector max(Rmax);
int n = x.size(), nf = f.size();//so essentially x and f are same thing for different use
int k = r[0], j=0, mx = max[0], loc[mx];
double y[n];
for (int i = 0; i < n; i++) y[i] = x[i];
//nth_element put smaller front, larger back
std::nth_element(y, y + k - 1, y + sizeof(y)/sizeof(*y));//then y become partial sorted ver of x
double  kl = y[k-1];//the kth smallest number, used as threshold
for (int i = 0; i < n; i++) y[i] = -x[i];
std::nth_element(y, y + k - 1, y + sizeof(y)/sizeof(*y));//then y become partial sorted ver of -x


double  ku = -y[k-1];//the kth largest number, used as threshold
for (int i = 0; i < nf; i++) {
if (f[i] <= kl || f[i] >= ku)
loc[j++] = i + 1;//keep those smaller than lower thresh and those larger than upper thresh
}
Rcpp::NumericVector rt(j);
for (int l = 0; l < j; l++) {
rt[l] = loc[l];//turn int array to numbericVector?
}
return rt;
'
getIdx <- cxxfunction(signature(Rx="numeric", Rf="numeric",
                                Rr="numeric", Rmax="numeric"),
                      src, plugin = "Rcpp")

  samp.data=function(data,beta,num.covariate,num.datalines,num.samp,delta)
  {
    dim=num.covariate
    x.matrix=data$x.matrix
    y=data$y
    samp.y=c()
    samp.x=data.frame()
    dim.samp=ceiling(num.samp/(2*(dim-1)))
    op.result=optimize(min.fun,num.covariate=num.covariate,interval=c(-10,10))
    c=op.result$minimum
    #print(c)
    # min.v=op.result$objective
    
    c_index=as.matrix(x.matrix)%*%beta
    distance1=abs(c_index-c)
    distance2=abs(c_index+c)
    index=(distance1>distance2)
    distance1[index]=0
    distance2[!index]=0
    distance=distance1+distance2
    
    index.update=(distance<delta)
    distance=distance[index.update]
    c_index=c_index[index.update];
    # print(c(max(c_index),min(c_index)))
    # ratio=min.v/min.fun(c_index[index.update[num]],num.covariate)
    samp.x=x.matrix[index.update,]
    # x.matrix=x.matrix[-index.update,]
    
    samp.y=y[index.update]
    x.matrix=samp.x
    y=samp.y
    
    u=length(y)
    print(u)
    gc()
    #print(u)
    samp.y=c()
    samp.x=c()
    #plot(x.matrix[,1],x.matrix[,2]);
    #for(i in 1:(dim-1))
    # {
    
    # index.update=order(x.matrix[,i])[c(1:dim.samp,(u-dim.samp+1):u)]
    #u=u-2*dim.samp
    # samp.y=c(samp.y,y[index.update])
    #samp.x=rbind(samp.x,x.matrix[index.update,])
    #y=y[-index.update]
    
    # x.matrix=x.matrix[-index.update,]
    #x.matrix=x.matrix[-]
    
    # }
    index <- c()
    use_column <- x.matrix[,1]
    # x.matrix_temp <- x.matrix
    index_raw<-c(1:nrow(x.matrix))
    for(i in 1:(dim-1))
    {
      index_pick <- getIdx(use_column, use_column, dim.samp, num.samp)
      index<- c(index,index_raw[index_pick])
      index_raw<- index_raw[-index_pick]
      use_column <- x.matrix[-index,i+1]
      #x.matrix <- x.matrix[-index_pick,]
    }
    
    #print(c(length(index),length(unique(index)),length(index)!=length(unique(index))))
    samp.y<-y[index]
    samp.x <- x.matrix[index,]
    
    print(length(samp.y))
    
    #plot(samp.x[,1],samp.x[,2])
    
    return(list(y=samp.y,x.matrix=samp.x))
  }
}


