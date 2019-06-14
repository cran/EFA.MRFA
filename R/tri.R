tri <-function(R){
#R -> simetric matrix

n <- dim(R)[2]

A <- c()

for (i in 1:n){
  for (j in (i+1):n){
    if (i==n){
      break
    }
    A = rbind(A,R[j,i])
  }
}

return(A)

}
