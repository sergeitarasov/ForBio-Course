
##############################################################################################
# Testing Correlated Character Evolution
#
# Aim: Test if simulated characters are (un)correlated using AIC
#
# 1. Simulate tree
# 2. Simulate character evolution under correlated model
# 3. Fit correlated and uncorrelated model
# 4. Compare models using AIC
#
#
##############################################################################################


library("corHMM")
library("phytools")

setwd("~/Documents/Courses/Bergen_2020/R")

# useful functions
source("SMM_functions.R")


#   ____________________________________________________________________________
#   Let's make a correlated character                                       ####


# initialize shape character
char.state<-c("T", "C")
rate.param<-c(0.1, 0.1)
TL<-init_char_matrix(char.state, rate.param, diag.as=0)
diag(TL)<-rate.param*-1
TL

# initialize color character
char.state<-c("r", "b")
rate.param<-c(0.1, 0.1)
COL<-init_char_matrix(char.state, rate.param, diag.as=0)
diag(COL)<-rate.param*-1
COL

Q.ind<-comb2matrices(TL, COL, controlling.state=2, name.sep="", diag.as="")
Q.ind

# make character correlated
Q.cor <- Q.ind
Q.cor[2,4] <- 0.3
Q.cor[3,4] <- 0.3
diag(Q.cor) <- 0
diag(Q.cor) <- -rowSums(Q.cor)
Q.cor

#   ____________________________________________________________________________
#   Simulate history                                                        ####


tree<-pbtree(n=100, scale=200)
sim.h<-sim.history(tree, Q.cor, nsim=1, anc=setNames(c(1, 0, 0, 0), colnames(Q.cor) ) )
plot(sim.h)


#   ____________________________________________________________________________
#   ML Inference                                                            ####


# matrix for inference
char.state<-c("T", "C")
rate.param<-c(1, 1)
TL<-init_char_matrix(char.state, rate.param, diag.as=0)

char.state<-c("r", "b")
rate.param<-c(1, 1)
COL<-init_char_matrix(char.state, rate.param, diag.as=0)

# uncorrelated matrix
M.equal<-comb2matrices(TL, COL,  name.sep="", diag.as="", non.rate.as=NA)
M.equal

# correlated matrix
M.corr <- M.equal
M.corr[2,4] <- 2
M.corr[3,4] <- 2
M.corr

#   ____________________________________________________________________________
#   reconstruct ancestral states using corHMM                               ####

# uncorrelated matrix
taxa <- cbind(hist$tip.label, sim.h$states)
out.equal <- rayDISC(sim.h, taxa, rate.mat=M.equal, node.states="marginal", model="ARD", root.p="maddfitz")

# correlated matrix
out.corr <- rayDISC(sim.h, taxa, rate.mat=M.corr, node.states="marginal", model="ARD", root.p="maddfitz")


#   ____________________________________________________________________________
#   Compare Results                                                         ####


out.equal
out.corr

out.equal$AIC-out.corr$AIC


