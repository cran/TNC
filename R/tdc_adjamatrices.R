.tdc_adjamatrices <- function(x, startsnapshot,
                endsnapshot,
                directed,
                normalize,
                centrality_evolution){
  if(!(as.integer(startsnapshot)==startsnapshot) |
     !(as.integer(endsnapshot)==endsnapshot)){
    stop("startsnapshot and endsnapshot have to be whole numbers")
  }
  if(endsnapshot<=startsnapshot | endsnapshot<2){
    stop("endsnapshot has to be bigger than startsnapshot and bigger than 1")
  }
  if(!is.list(x) | length(x)<2){
    stop("x has to be a list of adjacency matrices with minimum length of 2")
  }
  if(!all(c(sapply(x, min)) %in% c(0,1)) |
     !all(c(sapply(x, max)) %in% c(0,1))){
    stop("x must contain only 1 and 0")
  }
  if(!zero_range(c(sapply(x, dim)))){
    stop("Stop! Different number of nodes a different snapshots")
  }

  nV <- ncol(x[[1]])
  lapply(x,function(z){
    if(normalize==TRUE){
      apply(z, 1, function(y) {
        sum(y)/((nV-1)*(endsnapshot-startsnapshot))
      })
    }else{
      apply(z, 1, sum)
    }

  }
  )

  if(centrality_evolution==TRUE){
    return(
      list(
        TDC=Reduce("+",lapply(x,function(z){
          if(normalize==TRUE){
            apply(z, 1, function(y) {
              sum(y)/((nV-1)*(endsnapshot-startsnapshot))
            })
          }else{
            apply(z, 1, sum)
          }
        }
        )),
        CentEvo=t(do.call(rbind,lapply(x,function(z){
          if(normalize==TRUE){
            apply(z, 1, function(y) {
              sum(y)/((nV-1)*(endsnapshot-startsnapshot))
            })
          }else{
            apply(z, 1, sum)
          }
        }
        )))
      )
    )
  }else{
    return(
      Reduce("+",lapply(x,function(z){
        if(normalize==TRUE){
          apply(z, 1, function(y) {sum(y)/
              ((nV-1)*(endsnapshot-startsnapshot))})
        }else{
          apply(z, 1, sum)
        }
      }
      ))
    )
  }

}
