#' Temporal closeness centrality
#'
#' \code{tcc} returns the temporal closeness centrality for each node in a
#' dynamic network  (sequence of graph snapshots).
#'
#' @param x A list of adjacency matrices or a list of adjacency lists.
#' @param type Data format of \code{x}. Possible formats are \code{"M"} for a
#'   list of adjacency matrices (containing only 1s and 0s) and \code{"L"} for a
#'   list of adjacency lists (adjacency lists of the igraph package are
#'   supported). Default is \code{NULL}.
#' @param startsnapshot Numeric. Entry of \code{x} to start the calculation of
#'   \code{tcc}. Default is 1.
#' @param endsnapshot Numeric. Entry of \code{x} to end the calculation of
#'   \code{tcc}. Default is the last element of \code{x}.
#' @param vertexindices Numeric. A vector of nodes. Only shortest temporal paths
#'   ending at nodes in \code{vertexindices} are considered for calculating
#'   \code{tcc}. Can be used to parallel the calculation of \code{tcc} (see
#'   section Examples). Default is \code{NULL}.
#' @param directed Logical. Set \code{TRUE} if the dynamic network is a directed
#'   network. Default is \code{FALSE}.
#' @param normalize Logical. Set \code{TRUE} if centrality values should be
#'   normalized with \eqn{1/((|V|-1)*m)} where \eqn{|V|} is the number of nodes
#'   and \eqn{m =} \code{endsnapshot} \eqn{-} \code{startsnapshot}. Default is
#'   \code{TRUE}.
#' @param centrality_evolution Logical. Set \code{TRUE} if an additional matrix
#'   should be returned containing the centrality values at each snapshot. Rows
#'   correspondent to nodes, columns correspondent to snapshots. Default is
#'   \code{FALSE}.
#' @details \code{tcc} calculates the temporal closeness centrality (Kim and
#'   Anderson, 2012). To keep the computational effort linear in the number of
#'   snapshots the Reversed Evolution Network algorithm (REN; Hanke and Foraita,
#'   2017) is used to find all shortest temporal paths.
#' @section Warning: Using adjacency matrices as input exponentially increases
#'   the required memory. Use adjacency lists to save memory.
#' @return The (normalized) temporal betweenness centrality values of all nodes
#'   (\code{TCC}). If \code{centrality_evolution} is \code{TRUE}, an additional
#'   \eqn{(|V| x T)} matrix is returned (\code{CentEvo}), containing the
#'   temporal centrality value at each snapshot between \code{startsnapshot} and
#'   \code{endsnapshot}.
#' @seealso \code{\link{tbc},\link{tdc}}
#' @examples
#' # Create a list of adjacency matrices, plot the corresponding graphs
#' # (using the igraph package) and calculate tcc
#'
#' A1 <- matrix(c(0,1,0,0,0,0,
#'                1,0,1,0,0,0,
#'                0,1,0,0,0,0,
#'                0,0,0,0,0,0,
#'                0,0,0,0,0,0,
#'                0,0,0,0,0,0), ncol=6)
#'
#' A2 <- matrix(c(0,0,0,0,0,0,
#'                0,0,1,0,0,0,
#'                0,1,0,1,1,0,
#'                0,0,1,0,0,0,
#'                0,0,1,0,0,0,
#'                0,0,0,0,0,0), ncol=6)
#'
#' A3 <- matrix(c(0,0,0,0,0,0,
#'                0,0,0,0,0,0,
#'                0,0,0,0,0,0,
#'                0,0,0,0,0,0,
#'                0,0,0,0,0,0,
#'                0,0,0,0,0,0), ncol=6)
#'
#' A4 <- matrix(c(0,1,0,0,0,0,
#'                1,0,0,1,0,0,
#'                0,0,0,0,0,0,
#'                0,1,0,0,0,0,
#'                0,0,0,0,0,0,
#'                0,0,0,0,0,0), ncol=6)
#'
#' library(igraph)
#' par(mfrow=c(2,2))
#'
#' Layout <-
#'  layout_in_circle(graph_from_adjacency_matrix(A1, mode = "undirected"))
#'
#' plot(graph_from_adjacency_matrix(A1, "undirected"), layout=Layout)
#' plot(graph_from_adjacency_matrix(A2, "undirected"), layout=Layout)
#' plot(graph_from_adjacency_matrix(A3, "undirected"), layout=Layout)
#' plot(graph_from_adjacency_matrix(A4, "undirected"), layout=Layout)
#'
#' As <- list(A1,A2,A3,A4)
#'
#' tcc(As, "M", centrality_evolution=TRUE)
#'
#' ### Create list of adjacency lists
#' Ls <- lapply(seq_along(As), function(i){
#'   sapply(1:6, function(j){which(As[[i]][j,]==1)})
#' })
#'
#' tcc(Ls, "L", centrality_evolution=TRUE)
#'
#' ### Run tbc in parallel ###
#' library(parallel)
#' # Calculate the number of cores
#' cores_avail <- detectCores()-1
#' # Initiate cluster
#' cl <- makeCluster(2)
#' clusterExport(cl, c("As", "tcc"))
#'
#' TCC <- parLapply(cl, 1:6, function(x){
#'   tcc(As, "M", vertexindices = x)
#'  }
#' )
#'
#' stopCluster(cl)
#'
#' Reduce("+", TCC)
#' @references Kim, Hyoungshick and Anderson, Ross (2012). \emph{Temporal node
#'   centrality in complex networks}. Physical Review E, 85 (2).
#'
#'   Hanke, Moritz and Foraita, Ronja (2017). \emph{Clone temporal centrality
#'   measures for incomplete sequences of graph snapshots}. BMC Bioinformatics,
#'   18 (1).
#' @export

tcc <- function(x,
                type=NULL,
                startsnapshot = 1,
                endsnapshot = length(x),
                vertexindices = NULL,
                directed = FALSE,
                normalize = TRUE,
                centrality_evolution = FALSE){

  if(is.null(type)){
    stop("Please select a data format for the input x")
  }
  if(type=="L"){
    return(
      .tcc_adjalists(x, startsnapshot = startsnapshot,
                     endsnapshot = endsnapshot,
                     vertexindices = vertexindices,
                     directed = directed,
                     normalize = normalize,
                     centrality_evolution = centrality_evolution)
    )
  }
  if(type=="M"){
    return(
      .tcc_adjamatrices(x, startsnapshot = startsnapshot,
                        endsnapshot = endsnapshot,
                        vertexindices = vertexindices,
                        directed = directed,
                        normalize = normalize,
                        centrality_evolution = centrality_evolution)
    )
  }
}
