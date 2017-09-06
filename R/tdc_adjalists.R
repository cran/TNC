.tdc_adjalists <- function(x, startsnapshot,
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
  if(!zero_range(sapply(x, length))){
    stop("Stop! Different number of nodes a different snapshots")
  }


  nV <- length(x[[1]])
  lapply(x,function(z){
    if(normalize==TRUE){
      sapply(z, function(y) {
        length(y)/((nV-1)*(endsnapshot-startsnapshot))
      })
    }else{
      sapply(z, length)
    }

  }
  )

  if(centrality_evolution==TRUE){
    return(
      list(
        TDC=Reduce("+",lapply(x,function(z){
          if(normalize==TRUE){
            sapply(z, function(y) {
              length(y)/((nV-1)*(endsnapshot-startsnapshot))
            })
          }else{
            sapply(z, length)
          }
        }
        )),
        CentEvo=t(do.call(rbind,lapply(x,function(z){
          if(normalize==TRUE){
            sapply(z, function(y) {
              length(y)/((nV-1)*(endsnapshot-startsnapshot))
            })
          }else{
            sapply(z, length)
          }
        }
        )))
      )
    )
  }else{
    return(
      Reduce("+",lapply(x,function(z){
        if(normalize==TRUE){
          sapply(z, function(y) {
            length(y)/((nV-1)*(endsnapshot-startsnapshot))})
        }else{
          sapply(z, length)
        }
      }
      ))
    )
  }

}
