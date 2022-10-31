library(foreach)
library(doParallel)

# this is a bigger matrix that should give you time to notice the R processes that start up behind the scenes later
ourBigMatrix <- matrix(rnorm(100000), ncol = 5)

# we create a cluster (a group of R processes): make 4 of them, their communication type is "SOCK" (compatible with both Windows and Linux, probably MacOS too - don't worry about it). Bonus: skip loading the `methods` package for some extra speed - that's what the third argument is for :)
parallelCluster <- makeCluster(4, 
                               type = "SOCK", 
                               methods = FALSE)

# set the number of threads to 1 on each of the child processes
clusterEvalQ(cl = cl, {
  setMKLthreads(1)
})
a_list=NULL
system.time(for (i in 3e2:6e4) {
  output= sqrt(i)
  a_list= sum(a_list,output)
})
a_list
a_list=NULL
library(foreach)
system.time(foreach (i=3e2:6e4) %do% {
  output= sqrt(i)
  a_list= sum(a_list,output)
})
a_list