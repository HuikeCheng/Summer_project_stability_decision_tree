numrep <- 1000
association <- "Linear"
n <- c(50, 500)
pk <- 100
ev_xy <- c(0.01, 0.05, 0.1, 0.2)
nu_xy <- c(0.01, 0.05, 0.1, 0.2)
nchunks <- 50

a <- expand.grid(numrep=numrep, association=association, n=n, pk=pk, ev_xy=ev_xy, 
                 nu_xy=nu_xy, nchunks=nchunks)

colnames(a)

# write.table(a, file = "../Parameters/params_linear.txt",
#            quote = FALSE, row.names = FALSE, col.names = FALSE)

################################
numrep <- 1000
association <- c("Exponential", "Sigmoidal", "Indicator", "Quadratic")
n <- c(50, 500)
pk <- 100
ev_xy <- c(0.01, 0.05, 0.1, 0.2)
nu_xy <- c(0.01, 0.05, 0.1, 0.2)
nu_nl <- c(0.2, 0.6, 1)
nchunks <- 50

a <- expand.grid(numrep=numrep, association=association, n=n, pk=pk, ev_xy=ev_xy, 
                 nu_xy=nu_xy, nu_nl=nu_nl, nchunks=nchunks)

colnames(a)

# write.table(a, file = "../Parameters/params_nonlinear.txt",
#            quote = FALSE, row.names = FALSE, col.names = FALSE)

################################
numrep <- 1000
height <- c(3, 4, 5)
n <- c(50, 500)
pk <- 100
ev_xy <- c(0.01, 0.05, 0.1, 0.2)
nchunks <- 50

a <- expand.grid(numrep=numrep, height=height, n=n, pk=pk, ev_xy=ev_xy, 
                 nchunks=nchunks)

colnames(a)

# write.table(a, file = "../Parameters/params_tree.txt",
#            quote = FALSE, row.names = FALSE, col.names = FALSE)