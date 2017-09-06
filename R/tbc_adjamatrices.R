.tbc_adjamatrices <- function(x, startsnapshot,
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
  if(!all(c(sapply(x, min)) %in% c(0,1)) |
     !all(c(sapply(x, max)) %in% c(0,1))){
    stop("x must contain only 1 and 0")
  }
  if(!zero_range(c(sapply(x, dim)))){
    stop("Stop! Different number of nodes a different snapshots")
  }

  nV <- ncol(x[[1]])

  TBC <- rep(0, nV)

  if(centrality_evolution==TRUE){
    CentEvo <- matrix(0,nrow=nV, ncol=endsnapshot-startsnapshot+1)
  }

  if(is.null(vertexindices)){
    indices <- 1:nV
  }else{
    indices <- vertexindices
  }


  for(n in indices){

    endknot <- n

    TBCs_tminus1 <- matrix(rep(0,nV^2), ncol=nV)

    TBC_n <- vector(mode="list", length=endsnapshot-startsnapshot+1)

    TBC_t <- rep(0,nV)


    reached <- rep(Inf, nV)

    before_matrix <- matrix(rep(0,nV^2), ncol=nV)
    before_matrix[endknot,endknot] <- 1


    connected <- NULL
    for(j in (endsnapshot):startsnapshot){

      if(directed==TRUE){
        TMP <- t(x[[j]])
      }else{
        TMP <- x[[j]]
      }

      diag(TMP) <- 0

      TBC_t <- rep(0,nV)

      reached_old <- reached

      before_matrix_old <- before_matrix

      ones <- (1:nV)[TMP[endknot,]==1]

      before_matrix[,ones] <- rep(0,nV)
      before_matrix[endknot ,ones] <- 1
      diag(before_matrix)[ones] <- 1
      reached[ones] <- j

      TBCs_tminus1[ones,] <- 0


      if(length(connected)>0){

        if(length(ones)>(nV-3)){
          connecter <- (1:nV)[-c(ones, endknot)][sum(
            TMP[c(endknot, connected),-c(endknot, ones)])>0]
        }else{
          connecter <-
            (1:nV)[-c(endknot, ones)][apply(
              TMP[c(endknot, connected),-c(endknot, ones)],
              2,sum)>0]
        }




        for(i in connecter){

          i_connections <- (1:nV)[TMP[, i]==1]

          times_connections <- reached_old[i_connections]

          shortest_arrival <- min(times_connections)
          shortest_i_connections <-
            i_connections[times_connections==shortest_arrival]



          if((reached_old[i]==Inf) |
             (shortest_arrival < reached_old[i])){

            if(length(shortest_i_connections) == 1){
              before_matrix[,i] <-
                before_matrix_old[,shortest_i_connections]
            }else{
              before_matrix[,i] <- apply(
                before_matrix_old[,shortest_i_connections],
                1, sum)
            }


            Sum_SPs_t <- before_matrix[endknot,i]

            TBCs_tminus1[i,-c(endknot,i)] <-
              before_matrix[-c(endknot,i),i]/Sum_SPs_t



            before_matrix[i,i] <- Sum_SPs_t


            reached[i] <- shortest_arrival


          }

          if(shortest_arrival == reached_old[i]){

            if(length(shortest_i_connections)==1){
              before_matrix[,i] <- before_matrix[,i] +
                before_matrix_old[,shortest_i_connections]
            }else{
              before_matrix[,i] <- before_matrix[,i] +
                apply(before_matrix_old[,shortest_i_connections],
                      1, sum)
            }


            Sum_SPs_t <- before_matrix[endknot,i]


            TBCs_tminus1[i,-c(endknot,i)] <-
              before_matrix[-c(endknot,i),i]/Sum_SPs_t


            before_matrix[i,i] <- Sum_SPs_t

          }


        }

      }


      TBC_t <- apply(TBCs_tminus1, 2, sum)
      if(centrality_evolution==TRUE){
        CentEvo[,j] <- CentEvo[,j]+TBC_t
      }
      TBC_n[[j]] <- TBC_t

      connected <- (1:nV)[reached!=Inf]

    }

    TBC_n <- do.call(rbind, TBC_n)
    TBC_n <- apply(TBC_n, 2, sum)

    TBC <- TBC + TBC_n

  }

  if(normalize==TRUE){
    TBC <- TBC*(1/((nV-1)*(nV-2)*(endsnapshot-startsnapshot)))
  }

  if(centrality_evolution==TRUE){
    if(normalize==TRUE){
      out <-
        list(TBC=TBC,
             CentEvo=CentEvo*
               (1/((nV-1)*(nV-2)*(endsnapshot-startsnapshot))))
    }else{
      out <- list(TBC=TBC, CentEvo=CentEvo)
    }
    return(out)
  }
  else{
    return(TBC)
  }
}
