# File:         wheresCroc_rex.r 
# Description:  HMM 
#               for Project 2, Artificial Intelligence 2019, UU
# Author:       Rex Ruan

# Install the package
# install.packages("DeliveryMan_1.1.0.tar", repos = NULL, type="source")

# Load the library
library("WheresCroc")

normalize <- function(data) {
  t <- sum(data)
  return(data/t)
}

getEmissions <- function(readings,probs) {
  emissionsProb <- c()
  for(col in 1:40) {
      
      salinity <- dnorm(readings[1],probs[["salinity"]][col,1],probs[["salinity"]][col,2],F)
      phosphate <- dnorm(readings[2],probs[['phosphate']][col,1],probs[['phosphate']][col,2],F)
      nitrogen <- dnorm(readings[3],probs[['nitrogen']][col,1],probs[['nitrogen']][col,2],F)
      emissionsProb[col] <- prod(salinity,phosphate,nitrogen)
  }
  
  return(normalize(emissionsProb))
}

findNeighbors <- function(point,edges) {
  return (c(edges[which(edges[,1]==point),2],edges[which(edges[,2]==point),1]))
}

findPaths <- function(pointA,pointB,edges) {
  if(pointA==pointB){return(c(pointA))}
  visitedNode <- c(pointA)
  paths <- matrix(pointA)
  while(!pointB %in% visitedNode) {
    col <- dim(paths)[2]
    row <- dim(paths)[1]
    paths <- cbind(paths,0)
    while(row>0){
      newNeighbors <- setdiff(findNeighbors(paths[row,col],edges),visitedNode)
      if(length(newNeighbors)>0){
        visitedNode <- c(visitedNode,newNeighbors)
        paths[row,col+1] <- newNeighbors[1]
        for(num in tail(newNeighbors,-1)){
          paths <- rbind(paths,0)
          paths[dim(paths)[1],1:col] <- paths[row,1:col]
          paths[dim(paths)[1],col+1] <- num
        }
      }else{paths[row,col+1]<-0}
      row <- row-1
    }
    zeroRows <- which(paths[,col+1]==0)
    for(r in sort(zeroRows,decreasing=T)){
      paths <- paths[-r,]
    }
  }
  col <- dim(paths)[2]
  if(length(col)==0){return(paths)}
  return(c(paths[which(paths[,col]==pointB),]))
}

initialState <- function(backpackers) {
  if(backpackers[1] < 0 || backpackers[2] < 0){
    return(probDistrIfbpDeath(backpackers[which(backpackers<0)]))
  }
  initialProb <- replicate(40,0)
  initialProb <- initialProb + 1/38
  initialProb[backpackers[1]] <- 0
  initialProb[backpackers[2]] <- 0
  return(initialProb)
}

getTransProbs <- function(backpackers,prevTrans,edges){
  ##############compute transition probs###########
  transProbs <- replicate(40,0)
  for(i in 1:2){
    if(!is.na(backpackers[i])){
      if(backpackers[i]<0){
        prevTrans = transProbs
        prevTrans[backpackers[i]] = 1
        break
      } else {prevTrans[backpackers[i]] = 0} 
    }
  }
 
  for(i in 1:40){
    neighbors <- findNeighbors(i,edges)
    prob <- prevTrans[i] * 1/(length(neighbors)+1)
    for(neighbor in neighbors){
      num <- length(findNeighbors(neighbor,edges))+1
      prob <- prob + prevTrans[neighbor]/num
    }
    transProbs[i] <- prob
  }
  
  
  return(normalize(transProbs))
}

myFunction <- function(moveInfo, readings, positions, edges,probs){
  if(moveInfo$mem$status==0 || moveInfo$mem$status==1){
    moveInfo$mem$prevProbs <- initialState(positions[1:2])
    moveInfo$mem$status=2
    }
  
  emissions = getEmissions(readings,probs)
  transitions = getTransProbs(positions[1:2],moveInfo$mem$prevProbs,edges)
  hmm <- emissions*transitions
  
  index <- which.max(hmm)
  path <- findPaths(positions[3],index,edges)
 
  if(length(path)==1){
    moveInfo$moves <- c(0,0)
    hmm[index] <- 0
    moveInfo$mem$prevProbs <- hmm
    return(moveInfo)
  }
  
  if(length(path)==2){
    moveInfo$moves <- c(index,0)
    hmm[index] <- 0
    moveInfo$mem$prevProbs <- hmm
    return(moveInfo)
  }
  
  if(length(path)>2){
    moveInfo$moves <- path[2:3]
    moveInfo$mem$prevProbs <- hmm
    return(moveInfo)
  }
}