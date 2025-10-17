# auxiliar script for saving the data used in models 

X.MDQ <- list()
X.MHQ <- list()
X.MHO <- list()
X.MDO <- list()
for (station in estaciones){
  X.MDQ[[station]] <- MDQ[[station]]$X
  X.MDO[[station]] <- MDO[[station]]$X
  X.MHQ[[station]] <- MHQ[[station]]$X
  X.MHO[[station]] <- MHO[[station]]$X
}

qsave(X.MDQ, 'X.MDQ.qs')
qsave(X.MDO, 'X.MDO.qs')
qsave(X.MHQ, 'X.MHQ.qs')
qsave(X.MHO, 'X.MHO.qs')

