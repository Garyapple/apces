myDiscretized <- function(data){
  res = apply(data,2,function(x){
    tmp=array(1,length(x))
    tmp[x <= -0.5]=0
    tmp[x >= 0.5]=2
    x=tmp
  })
  return(res)
}