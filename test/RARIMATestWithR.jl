module RARIMATest
#-----------------------------------------------------------
#PURPOSE
#	Colin T. Bowers module for testing RARIMA.
#NOTES
#	This module needs to be combined with R-code placed on the R path, and so should only really be used by people who know exactly what they are doing
#	This module is not designed to test if RARIMA loads and the functions work (use runtests.jl for that). Rather it is used to compare output from ARIMA with output from manually porting data to R via csv file and calling the relevant functions, and making sure the two are "close"
#-----------------------------------------------------------


#Load any entire modules that are needed (use import ModuleName1, ModuleName2, etc)
using 	RARIMA


#Load any specific variables/functions that are needed (use import ModuleName1.FunctionName1, ModuleName2.FunctionName2, etc)

#Specify the variables/functions to export (use export FunctionName1, FunctionName2, etc)


#******************************************************************************

#----------------------------------------------------------
#SET CONSTANTS FOR MODULE
#----------------------------------------------------------
const testDataDir = "/home/" * string(ENV["USER"]) * "/.julia/v" * string(VERSION)[1:3] * "/RARIMA/test/"::ASCIIString
const testDataFPData = testDataDir * "x.csv"::ASCIIString
const testDataFPARIMAInput = testDataDir * "aI.csv"::ASCIIString
const testDataFPROutCoef = testDataDir * "ROutCoef.csv"::ASCIIString
const testDataFPROutOther = testDataDir * "ROutOther.csv"::ASCIIString
const testDataFPROutResid = testDataDir * "ROutResid.csv"::ASCIIString
const testDataFPROutvcov = testDataDir * "ROutvcov.csv"::ASCIIString



#Write an ARIMAInput to csv (should only ever be needed for testing purposes that I can forsee)
function writecsvARIMAInput(a::ARIMAInput)
	x = Array(ASCIIString, 12)
	x[1] = string(a.p)
	x[2] = string(a.d)
	x[3] = string(a.q)
	x[4] = string(a.P)
	x[5] = string(a.D)
	x[6] = string(a.Q)
	x[7] = string(a.seasonalPeriod)
	x[8] = string(a.includeIntercept)
	x[9] = string(a.transformParams)
	x[10] = replace(string(a.fixedParams), ",", ":")
	x[11] = replace(string(a.initParams), ",", ":")
	x[12] = string(a.method)
	writecsv(testDataFPARIMAInput, x)
end


#Simulate a random ARIMAInput
function simARIMAInput(bindMethod::Bool=false)
	p = rand(0:2)
	d = rand(0:2)
	q = rand(0:2)
	P = rand(0:2)
	D = rand(0:2)
	Q = rand(0:2)
	seasonalPeriod = rand(1:5)
	includeIntercept = convert(Bool, rand(0:1))
	transformParams = convert(Bool, rand(0:1))
	if includeIntercept == true
		if d > 0 || D > 0
			includeIntercept = false
		end
	end
	nC = p + q + P + Q
	nCinit = p + q + P + Q + convert(Int, includeIntercept)
	rand() > 0.99 ? fixedParams = 0.2 * (rand(nC) - 0.5) : fixedParams = Array(Float64, 0)
	rand() > 0.5 ? initParams = 0.2 * (rand(nC) - 0.5) : initParams = Array(Float64, 0)
	if bindMethod == false
		z = rand()
		if z < 0.4
			method = "CSS-ML"
		elseif z < 0.7
			method = "CSS"
		else
			method = "ML"
		end
	else
		method = "CSS-ML"
	end
	return(ARIMAInput(p, d, q, P, D, Q, seasonalPeriod, includeIntercept, transformParams, fixedParams, initParams, method))
end

#Test ARIMA estimate function by performing calculations independently in R
function testEstimateR(numIter::Int, numObs::Int)
	for kk = 1:numIter
		aI = simARIMAInput(true)
		z = rand()
		if z < 0.5
			x = randn(numObs)
		elseif z < 0.8
			x = 0.05 * rand(-100:100, numObs)
		else
			x = 0.05 * rand(1:1000, numObs)
		end
		writecsv(testDataFPData, x)
		writecsvARIMAInput(aI)
		aM = estimate(x, aI) #estimate model in Julia
		RSystemCallString = `Rscript /home/colin/TestData/Julia_ARIMATest_System_Call_For_R_1.R`
		run(RSystemCallString)
		RCoef = vec(readcsv(testDataFPROutCoef))
		ROther = vec(readcsv(testDataFPROutOther))
		RResid = vec(readcsv(testDataFPROutResid))
		Rvcov = readcsv(testDataFPROutvcov)
		sumabs(RCoef - coef(aM)) > 1e-12 && error("Coefficient estimates deviate by large amount")
		sumabs(ROther - [aM.logLik, aM.aic]) > 1e-10 && error("log-likelihood and aic estimates deviate by large amount")
		sumabs(RResid - aM.resid) > 1e-10 && error("Residuals deviate by large amount")
		sumabs(Rvcov - aM.covCoef) > 1e-8 && error("coefficient covariance matrix estimates deviate by large amount")
	end
	println("Test passed.")
end



function simARIMAModel(; d::Int=-1, includeIntercept::Bool=false, includeSeasonal::Bool=false)
	p = rand(0:2)
	d < 0 && (d = rand(0:2))
	q = rand(0:2)
	includeSeasonal ? (D = rand(0:1); P = rand(0:1); Q = rand(0:1)) : (D = 0; P = 0; Q = 0)
	arCoef = Array(Float64, 0)
	maCoef = Array(Float64, 0)
	sarCoef = Array(Float64, 0)
	smaCoef = Array(Float64, 0)
	for k = 1:p
		rand() > 0.5 ? (mult = 1) : (mult = -1)
		push!(arCoef, mult * 0.5 * (p / k) * rand() * 0.5)
	end
	for k = 1:q
		rand() > 0.5 ? (mult = 1) : (mult = -1)
		push!(maCoef, mult * 0.5 * (p / k) * rand() * 0.5)
	end
	rand() > 0.5 ? (mult = 1) : (mult = -1)
	P == 1 && push(sarCoef, mult * 0.5 * rand())
	Q == 1 && push(smaCoef, mult * 0.5 * rand())
	if includeIntercept == true
		d == 0 ? (intercept = randn()) : (intercept = NaN)
	else
		intercept = NaN
	end
	sigma2 = 3 * rand()
	numCoef = p + q + P + Q + convert(Int, !(isnan(intercept)))
	covCoef = zeros(Float64, numCoef, numCoef)
	covCoefName = Array(ASCIIString, numCoef)
	seasonalPeriod = 1
	logLik = 0.0
	aic = 0.0
	resid = Array(Float64, 0)
	optimCode = 0
	rModelName = ""
	return(ARIMAModel(arCoef, maCoef, sarCoef, smaCoef, intercept, sigma2, covCoef, covCoefName, d, D, seasonalPeriod, logLik, aic, resid, optimCode, rModelName))
end





function testSimulateARIMA(numIter::Int)
	for kk = 1:numIter
		numObs = rand(10:1000)
		numBurn = rand(10:100)
		a = simARIMAModel(d=0, includeIntercept=false, includeSeasonal=false)
		globalEnv[:numObs] = numObs
		globalEnv[:numBurn] = numBurn
		globalEnv[:arCoef] = a.arCoef
		globalEnv[:maCoef] = a.maCoef
		reval("eVec = rnorm(numObs + 100)")
		reval("eBurn = rnorm(numBurn + 100)")
		xR = rcopy("as.numeric(arima.sim(n=numObs, list(ar=arCoef, ma=maCoef), innov=eVec, n.start=numBurn, start.innov=eBurn))")
		eVec = rcopy("eVec")
		eBurn = rcopy("eBurn")
		x = simulate(numObs, a, numBurn)



	end
end






end #end Module ctbRARIMA
