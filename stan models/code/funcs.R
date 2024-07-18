make_design_matrix=function(x,grp){
  x2=matrix(nrow=length(x),ncol=length(unique(grp)))
  for(i in 1:length(unique(grp))){
    x2[,i]=ifelse(grp==levels(factor(grp))[i],1,0)*x
  }
  return(x2)
}

rag_n=function(x){
  rle=rle(x) #get run lengths of strings
  start_n=c(1,1+cumsum(rle$lengths)) #starting points
  start_n=start_n[-length(start_n)] #drop last value
  end_n=c(start_n[-1]-1) #end points
  end_n[length(start_n)]=length(x) #final obs.
  
  return(cbind(start_n,end_n))
}
