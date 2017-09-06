.tcc_adjalists <- function(x, startsnapshot,
                endsnapshot,
                vertexindices,
                directed,
                normalize,
                centrality_evolution){
  if(!(as.integer(startsnapshot)==startsnapshot) |
     !(as.integer(endsnapshot)==endsnapshot)){
    stop("startsnapshot and endsnapshot must be integers")
  }
  if(endsnapshot<=startsnapshot | endsnapshot<2){
    stop("endsnapshot needs to be larger than 1 and larger than startsnapshot")
  }
  if(!is.list(x) | length(x)<2){
    stop("x needs to be a list of adjacency matrices with minimum length of 2")
  }
  if(!zero_range(sapply(x, length))){
    stop("Stop! Different number of nodes a different snapshots")
  }

  nV <- length(x[[1]])

  if(centrality_evolution==TRUE){
    CentEvo <- matrix(0,nrow=nV, ncol=endsnapshot-startsnapshot+1)
  }
  TCC <- rep(0, nV)
  if(is.null(vertexindices)){
    indices <- 1:nV
  }else{
    indices <- vertexindices
  }
  for(n in indices){
    endknot <- n

    TCC_t <- rep(0,nV)

    reached <- rep(Inf, nV)

    for(j in (endsnapshot):startsnapshot){

      direct_connected <- as.vector(x[[j]][[endknot]])
      connected_before <- which(reached!=Inf)
      for(i in connected_before){
        indirect_connected <- as.vector(x[[j]][[i]])
        indirect_connected <-
          indirect_connected[which(indirect_connected != endknot &
                                     !(indirect_connected %in% direct_connected))]
        for(l in indirect_connected){
          reached[l] <-
            min(reached[as.vector(x[[j]][[l]])], reached[l])
        }
      }
      reached[direct_connected] <- j
      TCC_t <- TCC_t + 1/(reached-j+1)
      if(centrality_evolution==TRUE){

        CentEvo[,j] <-
          CentEvo[,j]+(1/(reached-j+1))/((nV-1)*(endsnapshot-startsnapshot))
      }
    }
    TCC <- TCC+TCC_t
  }
  if(normalize==TRUE){
    TCC <- TCC/((nV-1)*(endsnapshot-startsnapshot))
  }
  if(centrality_evolution==TRUE){
    out <- list(TCC=TCC, CentEvo=CentEvo)
    return(out)
  }
  else{
    return(TCC)
  }
}
