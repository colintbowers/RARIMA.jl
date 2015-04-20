
#Set fixed locations
testDataDir = "/home/colin/TestData/"
testDataFPData = paste(testDataDir, "x.csv", sep="")
testDataFPARIMAInput = paste(testDataDir, "aI.csv", sep="")
testDataFPRCoefOutput = paste(testDataDir, "ROutCoef.csv", sep="")
testDataFPRvcovOutput = paste(testDataDir, "ROutvcov.csv", sep="")
testDataFPROtherOutput = paste(testDataDir, "ROutOther.csv", sep="")
testDataFPRResidOutput = paste(testDataDir, "ROutResid.csv", sep="")


#Read in the data
x = read.csv(testDataFPData, header=FALSE)
x = ts(x$V1)

#Read in model specs
y = scan(testDataFPARIMAInput, what="", sep="\n")

#Get relevant variables from model specs

#order
o1 = c(as.numeric(y[1]), as.numeric(y[2]), as.numeric(y[3]))

#seasonal order
so1 = list(order=c(as.numeric(y[4]), as.numeric(y[5]), as.numeric(y[6])), period=as.numeric(y[7]))

#include.mean
if (identical("true", y[8])) {
    im=TRUE
} else {
    im=FALSE
}

#transform.pars
if (identical("true", y[9])) {
    tP=TRUE
} else {
    tP=FALSE
}

#fixed
if (identical(y[10], "[]")) {
    fxd=NULL
} else {
    fxdTemp = strsplit(y[10], ":", fixed=TRUE)
    fxdTemp[[1]][1] = substr(fxdTemp[[1]][1], 2, nchar(fxdTemp[[1]][1]))
    K = length(fxdTemp[[1]])
    if (K > 1) {
        fxdTemp[[1]][K] = substr(fxdTemp[[1]][K], 1, nchar(fxdTemp[[1]][K])-1)
    }
    fxd = vector(mode="double", length=K)
    for (k in 1:K) {
        fxd[k] = as.numeric(fxdTemp[[1]][k])
    }
}

#init
if (identical(y[11], "[]")) {
    init1=NULL
} else {
    init1Temp = strsplit(y[11], ":", fixed=TRUE)
    init1Temp[[1]][1] = substr(init1Temp[[1]][1], 2, nchar(init1Temp[[1]][1]))
    K = length(init1Temp[[1]])
    if (K > 1) {
        init1Temp[[1]][K] = substr(init1Temp[[1]][K], 1, nchar(init1Temp[[1]][K])-1)
    }
    init1 = vector(mode="double", length=K)
    for (k in 1:K) {
        init1[k] = as.numeric(init1Temp[[1]][k])
    }
}

#method
method1 = y[12]

#Estimate arima model
m1 = arima(x, order=o1, seasonal=so1, include.mean=im, transform.pars=tP, fixed=fxd[1:6], init=init1, method=method1)

#Export numbers for testing
write.table(as.numeric(m1$coef), file=testDataFPRCoefOutput, row.names=FALSE, col.names=FALSE, sep=",")
write.table(m1$var.coef, file=testDataFPRvcovOutput, row.names=FALSE, col.names=FALSE, sep=",")
write.table(m1$residuals, testDataFPRResidOutput, row.names=FALSE, col.names=FALSE, sep=",")
otherVec = c(m1$loglik, m1$aic)
write.table(otherVec, testDataFPROtherOutput, row.names=FALSE, col.names=FALSE, sep=",")