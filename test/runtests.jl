using RARIMA

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


function testBasic(numIter::Int)
	for kk = 1:numIter
		aI = simARIMAInput(true)
		aMSim = simARIMAModel(includeIntercept=false, includeSeasonal=false)
		numObs = rand(50:1000)
		x = simulate(numObs, aMSim) #test simulate in ARIMAModel
		aM = estimate(x, aI) #test estimate on ARIMAInput
		aAI = AutoARIMAInput()
		aMAuto = estimate(x, aAI) #test estimate on AutoARIMAInput
		if aM.rModelName != ""
			(pointFore, boundFore) = forecast(aM, 5) #test forecast using R
		end
		if aMAuto.D == 0 && length(aMAuto.smaCoef) == 0 && length(aMAuto.sarCoef) == 0
			(pointFore, boundFore) = forecast(x, aMAuto, 5) #test forecast using Julia
		end
	end
	println("Test passed.")
end


#Perform the tests for 100 iterations
testBasic(1)
testBasic(100)

