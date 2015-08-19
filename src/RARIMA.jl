module RARIMA
#-----------------------------------------------------------
#PURPOSE
#	Colin T. Bowers module for ARIMA modelling
#NOTES
#	This module wraps ARIMA functionality in the stats package of the R programming language as well as the forecast library of Rob Hyndman
#	This module uses Douglas Bates RCall package to interact with R
#IMPLEMENTED IN R
#	all estimate methods call R
#	forecast method that provides confidence bounds currently calls R
#IMPLEMENTED IN JULIA
#	forecast method that does NOT provide confidence bounds
#	simulate method
#SUB-OPTIMAL ROUTINES
#	The arima_statdiff and arima_statdiffinv routines create a temporary array for every order of differencing/integration. There is probably a smarter way to do this.
#LICENSE
#	Since I wrap the forecast library of Rob Hyndman, an MIT license would be invalid, so instead this module is under the GNU GPL v2. See github repository for more detail: https://github.com/colintbowers/RARIMA.jl.git
-----------------------------------------------------------


#Load any entire modules that are needed (use import ModuleName1, ModuleName2, etc)
using 	RCall

#Load an libraries needed by R (use reval("library(libraryNameHere)") to do this)
reval("library(forecast)") #needed for the auto.arima function

#Load any specific variables/functions that are needed (use import ModuleName1.FunctionName1, ModuleName2.FunctionName2, etc)
import 	Base.show,
  		Base.copy,
		Base.deepcopy

#Specify the variables/functions to export (use export FunctionName1, FunctionName2, etc)
export	ARIMAInput, #Type that specifies input parameters for estimating an ARIMA model
		AutoARIMAInput, #Type that specifies input parameters for automatically selecting an ARIMA model then estimating it
		ARIMAModel, #Output of estimation (or auto-estimation) of an ARIMA model, ie complete specification of a model
		estimate, #Function for estimating a model given some data
		estimateARIMA, #Keyword wrapper on estimate function
		estimateAutoARIMA, #Keyword wrapper on estimate function
		forecast, #Function for forecasting using a model and data
		simulate, #Function for simulating a model
		simulateARIMA, #Keyword wrapper on simulate function
		coef, #Function for retrieving coefficients of ARIMAModel
		vcov, #Function for retrieving variance-covariance matrix of parameters of ARIMAModel
		coefVar #Function for retrieving coefficients of ARIMAModel along with their variance estimates



#******************************************************************************

#----------------------------------------------------------
#SET CONSTANTS FOR MODULE
#----------------------------------------------------------
const validMethodType = Union(ASCIIString, Symbol)::UnionType #type union nesting valid ways to express an optimization/estimation method



#----------------------------------------------------------
#SET ABSTRACT TYPES FOR MODULE
#----------------------------------------------------------
#None


#----------------------------------------------------------
#TYPE
#	ARIMAInput
#FIELDS
#	p::Int: Number of AR lags
#	d::Int: Order of integration
#	q::Int: Number of MA lags
#	P::Int: Number of seasonal AR lags
#	D::Int: Order of seasonal integration
#	Q::Int: Number of seasonal MA lags
#	seasonalPeriod::Int: Seasonal period of data (only relevant if P, D, or Q are > 0)
#	includeIntercept::Bool: Include intercept term in model. Not relevant for d > 0.
#	transformParams::Bool: AR parameters transformed to ensure they remain in region of stationarity
#	fixedParams::Vector{Float64}: Array of size p+q+P+Q containing fixed values for coefficients. Set value to NaN if you don't want to fix a particular parameter. Set to empty array to not fix any parameters during estimation.
#	initParams::Vector{Float64}: Array of size p+q+P+Q of initial values for estimation. Set value to NaN if you don't have a initial value for a particular parameter (0 will be used). Set to empty array to not supply any initial values.
#	method<:validMethodType: Variable that indicates method of estimating parameters. If ASCIIString, set to "CSS-ML", "CSS", or "ML". See R docs for arima function for more detail. Other types for this variable not coded up yet.
#PURPOSE
#	This type stores all of the input, excluding the data, needed to estimate the parameters of an ARIMA model.
#CONSTRUCTORS
#	ARIMAInput{T<:validMethodType}(p::Int, d::Int, q::Int, P::Int, D::Int, Q::Int, seasonalPeriod::Int, includeIntercept::Bool, transformParams::Bool, fixedParams::Vector{Float64}, initParams::Vector{Float64}, method::T): Constructor for every field
#	ARIMAInput(): White noise constructor
#KEYWORD CONSTRUCTOR
#	ARIMAInput(; p::Int=0, d::Int=0, q::Int=0, P::Int=0, D::Int=0, Q::Int=0, seasonalPeriod::Int=1, includeIntercept::Bool=true, transformParams::Bool=true, fixedParams::Vector{Float64}=Array(Float64, 0), initParams::Vector{Float64}=Array(Float64, 0), method::ASCIIString="CSS-ML")
#METHODS
#	copy(a::ARIMAInput): Standard copy
#	deepcopy(a::ARIMAInput): Standard deepcopy
#	show(io::IO, a::ARIMAInput): Standard show
#	show(a::ARIMAInput): Show defaulting to STDOUT
#NOTES
#	Any new estimation methods can be allowed for by using a custom type for the method field
#----------------------------------------------------------
#--------------------
#Type definition
#--------------------
type ARIMAInput{T<:validMethodType}
	p::Int
	d::Int
	q::Int
	P::Int
	D::Int
	Q::Int
	seasonalPeriod::Int
	includeIntercept::Bool
	transformParams::Bool
	fixedParams::Vector{Float64}
	initParams::Vector{Float64}
	method::T
	function ARIMAInput{T<:validMethodType}(p::Int, d::Int, q::Int, P::Int, D::Int, Q::Int, seasonalPeriod::Int, includeIntercept::Bool, transformParams::Bool, fixedParams::Vector{Float64}, initParams::Vector{Float64}, method::T)
		(p < 0 || d < 0 || q < 0 || P < 0 || D < 0 || Q < 0) && error("p, d, q, P, D, and Q must all be non-negative")
		seasonalPeriod < 0 && error("Seasonal period must be non-negative")
		if includeIntercept == true
			if d > 0 || D > 0
				println("WARNING: includeIntercept altered to false since either d > 0 or D > 0")
				includeIntercept = false
			end
		end
		numCoefTemp = p + q + P + Q
		if length(fixedParams) > 0
			length(fixedParams) != numCoefTemp && error("Fixed parameter vector is different length to number of model coefficients")
		end
		if length(initParams) > 0
			length(initParams) != numCoefTemp && error("Initialization parameter vector is different length to number of model coefficients")
		end
		if typeof(method) == ASCIIString
			(method != "CSS-ML" && method != "ML" && method != "CSS") && error("Invalid estimation method")
		elseif typeof(method) == Symbol
			error("method field of type Symbol has not been coded yet")
		else
			error("method field is an unexpected type")
		end
		new(p, d, q, P, D, Q, seasonalPeriod, includeIntercept, transformParams, fixedParams, initParams, method)
	end
end
#--------------------
#Constructors (Note, one constructor is defined after ARIMAModel type definition)
#--------------------
ARIMAInput{T<:validMethodType}(p::Int, d::Int, q::Int, P::Int, D::Int, Q::Int, seasonalPeriod::Int, includeIntercept::Bool, transformParams::Bool, fixedParams::Vector{Float64}, initParams::Vector{Float64}, method::T) = ARIMAInput{typeof(method)}(p, d, q, P, D, Q, seasonalPeriod, includeIntercept, transformParams, fixedParams, initParams, method)
ARIMAInput() = ARIMAInput(0, 0, 0, 0, 0, 0, 1, true, true, Array(Float64, 0), Array(Float64, 0), "CSS-ML")
#--------------------
#Keyword constructor
#--------------------
ARIMAInput(; p::Int=0, d::Int=0, q::Int=0, P::Int=0, D::Int=0, Q::Int=0, seasonalPeriod::Int=1, includeIntercept::Bool=true, transformParams::Bool=true, fixedParams::Vector{Float64}=Array(Float64, 0), initParams::Vector{Float64}=Array(Float64, 0), method::ASCIIString="CSS-ML") = ARIMAInput(p, d, q, P, D, Q, seasonalPeriod, includeIntercept, transformParams, fixedParams, initParams, method)
#--------------------
#Methods
#--------------------
#copy, deepcopy
copy(a::ARIMAInput) = ARIMAInput(copy(a.p), copy(a.d), copy(a.q), copy(a.P), copy(a.D), copy(a.Q), copy(a.seasonalPeriod), copy(a.includeIntercept), copy(a.transformParams), copy(a.fixedParams), copy(a.initParams), copy(a.method))
deepcopy(a::ARIMAInput) = ARIMAInput(deepcopy(a.p), deepcopy(a.d), deepcopy(a.q), deepcopy(a.P), deepcopy(a.D), deepcopy(a.Q), deepcopy(a.seasonalPeriod), deepcopy(a.includeIntercept), deepcopy(a.transformParams), deepcopy(a.fixedParams), deepcopy(a.initParams), deepcopy(a.method))
#show
function show(io::IO, a::ARIMAInput)
	println(io, "ARIMAInput contents:")
	println(io, "    p = " * string(a.p))
	println(io, "    d = " * string(a.d))
	println(io, "    q = " * string(a.q))
	println(io, "    P = " * string(a.P))
	println(io, "    D = " * string(a.D))
	println(io, "    Q = " * string(a.Q))
	println(io, "    seasonal period = " * string(a.seasonalPeriod))
	println(io, "    include intercept = " * string(a.includeIntercept))
	println(io, "    transform params for stationarity = " * string(a.transformParams))
	if length(a.fixedParams) == 0
		println(io, "    fixed parameters = none")
	else
		println(io, "    fixed parameters = some allocated")
	end
	if length(a.initParams) == 0
		println(io, "    initial parameters = none")
	else
		println(io, "    initial parameters = some allocated")
	end
	println(io, "    optimisation method = " * string(a.method))
end
show(a::ARIMAInput) = show(STDOUT, a)






#----------------------------------------------------------
#TYPE
#	AutoARIMAInput
#FIELDS
# 	d::Int: order of integration (use to fix a specific value). Set to negative value if d not fixed.
# 	D::Int: order of seasonal integration (use to fix a specific value). Set to negative value if D not fixed.
# 	pMax::Int: maximum number of AR lags allowed in search procedure
# 	qMax::Int: maximum number of MA lags allowed in search procedure
# 	PMax::Int: maximum number of seasonal AR lags allowed in search procedure
# 	QMax::Int: maximum number of seasonal MA lags allowed in search procedure
# 	numCoefMax::Int: maximum number of coefficients total allowed in search procedure (not relevant for step-wise search procedure)
# 	dMax::Int: maximum order of integration allowed in search procedure (over-ridden by d)
# 	DMax::Int: maximum seasonal order of integration allowed in search procedure (over-ridden by D)
# 	pStart::Int: Initial value for p in step-wise search procedure
# 	qStart::Int: Initial value for q in step-wise search procedure
# 	PStart::Int: Initial value for P in step-wise search procedure
# 	QStart::Int: Initial value for Q in step-wise search procedure
# 	stationary::Bool: Set to true to restrict search procedure to stationary models
# 	seasonal::Bool: Set to false to restrict search procedure to non-seasonal models
# 	ic::ASCIIString: Information criteria to use in search procedure. Use "aic", "aicc", or "bic" (see R auto-arima docs for more detail)
# 	stepwise::Bool: If true, use step-wise search procedure. If false, check every feasible model (slower).
# 	parallel::Bool: If true and stepwise=false, then do checks in parallel if possible.
# 	numCore::Int: Number of cores to use if doing checks in parallel.
#PURPOSE
#	This type stores all of the input, excluding the data, needed to automatically select an ARIMA model and estimate it
#CONSTRUCTORS
#	AutoARIMAInput(d::Int, D::Int, pMax::Int, qMax::Int, PMax::Int, QMax::Int, numCoefMax::Int, dMax::Int, DMax::Int, pStart::Int, qStart::Int, PStart::Int, QStart::Int, stationary::Bool, seasonal::Bool): Default methodology options
#	AutoARIMAInput(): No seasonality, all method options set to default, all order and step-wise start values set to default.
#KEYWORD CONSTRUCTOR
#	AutoARIMAInput(; d::Int=-1, D::Int=-1, pMax::Int=5, qMax::Int=5, PMax::Int=2, QMax::Int=2, numCoefMax::Int=5, dMax::Int=2, DMax::Int=1, pStart::Int=2, qStart::Int=2, PStart::Int=1, QStart::Int=5, stationary::Bool=false, seasonal::Bool=false, ic::ASCIIString="aic", stepwise::Bool=true, parallel::Bool=false, numCore::Int=2)
#METHODS
#	copy(a::AutoARIMAInput): Standard copy
#	deepcopy(a::AutoARIMAInput): Standard deepcopy
#	show(io::IO, a::AutoARIMAInput): Standard show
#	show(a::AutoARIMAInput): Show defaulting to STDOUT
#NOTES
#----------------------------------------------------------
#--------------------
#Type definition
#--------------------
type AutoARIMAInput
	d::Int
	D::Int
	pMax::Int
	qMax::Int
	PMax::Int
	QMax::Int
	numCoefMax::Int
	dMax::Int
	DMax::Int
	pStart::Int
	qStart::Int
	PStart::Int
	QStart::Int
	stationary::Bool
	seasonal::Bool
	ic::ASCIIString
	stepwise::Bool
	parallel::Bool
	numCore::Int
	function AutoARIMAInput(d::Int, D::Int, pMax::Int, qMax::Int, PMax::Int, QMax::Int, numCoefMax::Int, dMax::Int, DMax::Int, pStart::Int, qStart::Int, PStart::Int, QStart::Int, stationary::Bool, seasonal::Bool, ic::ASCIIString, stepwise::Bool, parallel::Bool, numCore::Int)
		(pMax < 0 || qMax < 0 || PMax < 0 || QMax < 0 || numCoefMax < 0 || dMax < 0 || DMax < 0) && error("maximum value non-negativity constraint violated in inner constructor")
		(pStart < 0 || qStart < 0 || PStart < 0 || QStart < 0) && error("step-wise starting value non-negativity constraint violated in inner constructor")
		(ic != "aicc" && ic != "aic" && ic != "bic") && error("invalid value for ic field")
		if parallel == true && stepwise == true
			println("WARNING: parallel computation incompatible with stepwise search method. parallel altered to false.")
			parallel = false
		end
		numCore < 1 && error("numCore must be strictly positive")
		new(d, D, pMax, qMax, PMax, QMax, numCoefMax, dMax, DMax, pStart, qStart, PStart, QStart, stationary, seasonal, ic, stepwise, parallel, numCore)
	end
end
#--------------------
#Constructors
#--------------------
AutoARIMAInput() = AutoARIMAInput((-1, -1, 5, 5, 2, 2, 5, 2, 1, 2, 2, 1, 1, false, false, "aic", true, false, 2))
#--------------------
#Keyword constructor
#--------------------
AutoARIMAInput(; d::Int=-1, D::Int=-1, pMax::Int=5, qMax::Int=5, PMax::Int=2, QMax::Int=2, numCoefMax::Int=5, dMax::Int=2, DMax::Int=1, pStart::Int=2, qStart::Int=2, PStart::Int=1, QStart::Int=5, stationary::Bool=false, seasonal::Bool=false, ic::ASCIIString="aic", stepwise::Bool=true, parallel::Bool=false, numCore::Int=2) = AutoARIMAInput(d, D, pMax, qMax, PMax, QMax, numCoefMax, dMax, DMax, pStart, qStart, PStart, QStart, stationary, seasonal, ic, stepwise, parallel, numCore)
#--------------------
#Methods
#--------------------
#copy, deepcopy
copy(a::AutoARIMAInput) = AutoARIMAInput(copy(a.d), copy(a.D), copy(a.pMax), copy(a.qMax), copy(a.PMax), copy(a.QMax), copy(a.numCoefMax), copy(a.dMax), copy(a.DMax), copy(a.pStart), copy(a.qStart), copy(a.PStart), copy(a.QStart), copy(a.stationary), copy(a.seasonal), copy(a.ic), copy(a.stepwise), copy(a.parallel), copy(a.numCore))
deepcopy(a::AutoARIMAInput) = AutoARIMAInput(deepcopy(a.d), deepcopy(a.D), deepcopy(a.pMax), deepcopy(a.qMax), deepcopy(a.PMax), deepcopy(a.QMax), deepcopy(a.numCoefMax), deepcopy(a.dMax), deepcopy(a.DMax), deepcopy(a.pStart), deepcopy(a.qStart), deepcopy(a.PStart), deepcopy(a.QStart), deepcopy(a.stationary), deepcopy(a.seasonal), deepcopy(a.ic), deepcopy(a.stepwise), deepcopy(a.parallel), deepcopy(a.numCore))
#show
function show(io::IO, a::AutoARIMAInput)
	println(io, "AutoARIMAInput contents:")
	a.d < 0 ? println(io, "    order of integration = not set") : println(io, "    order of integration = " * string(a.d))
	a.D < 0 ? println(io, "    seasonal order of integration = not set") : println(io, "    seasonal order of integration = " * string(a.D))
	println(io, "    pMax = " * string(a.pMax))
	println(io, "    qMax = " * string(a.qMax))
	println(io, "    PMax = " * string(a.PMax))
	println(io, "    QMax = " * string(a.QMax))
	println(io, "    numCoefMax = " * string(a.numCoefMax))
	println(io, "    dMax = " * string(a.dMax))
	println(io, "    DMax = " * string(a.DMax))
	println(io, "    pStart = " * string(a.pStart))
	println(io, "    qStart = " * string(a.qStart))
	println(io, "    PStart = " * string(a.PStart))
	println(io, "    QStart = " * string(a.QStart))
	println(io, "    restrict search to stationary models = " * string(a.stationary))
	println(io, "    include seasonal models in search path = " * string(a.seasonal))
	println(io, "    information criterion = " * a.ic)
	println(io, "    use step-wise search procedure = " * string(a.stepwise))
	println(io, "    use parallel computation (only possible when step-wise is false) = " * string(a.parallel))
	println(io, "    number of cores for parallel use = " * string(a.numCore))
end
show(a::AutoARIMAInput) = show(STDOUT, a)








#----------------------------------------------------------
#TYPE
#	ARIMAModel
#FIELDS
# 	arCoef::Vector{Float64}: Vector of AR coefficients
# 	maCoef::Vector{Float64}: Vector of MA coefficients
# 	sarCoef::Vector{Float64}: Vector of seasonal AR coefficients
# 	smaCoef::Vector{Float64}: Vector of seasonal MA coefficients
# 	intercept::Float64: Model intercept. Equal to NaN if no intercept.
# 	sigma2::Float64: Variance of residuals
# 	covCoef::Matrix{Float64}: Covariance matrix of estimated coefficients (and intercept if applicable)
# 	covCoefName::Vector{ASCIIString}: Names corresponding to the rows/columns of covCoef
# 	d::Int: Order of integration
# 	D::Int: Seasonal order of integration
# 	seasonalPeriod::Int: Seasonal period.
# 	logLik::Float64: Log-likelihood at maximum as reported by R in estimation step.
# 	aic::Float64: aic at maximum as reported by R in estimation step.
# 	resid::Vector{Float64}: Vector of residuals from estimation step.
# 	optimCode::Int: Optimisation output code as reported by R. 0=normal convergence obtained, other values imply something wrong (see R docs for more detail)
# 	rModelName::ASCIIString: Variable name assigned to model output from R estimation step in embedded R workspace. This field is used to find the R output from estimation step and feed it into R predict.Arima (if the user calls Julia forecast function)
#PURPOSE
#	This type stores a complete ARIMA specification.
#CONSTRUCTORS
#	ARIMAModel(): Default values
#KEYWORD CONSTRUCTOR
#	ARIMAModel(; arCoef::Vector{Float64}=Array(Float64, 0), maCoef::Vector{Float64}=Array(Float64, 0), sarCoef::Vector{Float64}=Array(Float64, 0), smaCoef::Vector{Float64}=Array(Float64, 0), intercept::Float64=0.0, sigma2::Float64=1.0, covCoef::Matrix{Float64}==Array(Float64, 0, 0), covCoefName::Vector{ASCIIString}=Array(ASCIIString, 0), d::Int=0, D::Int=0, seasonalPeriod::Int=1, logLik::Float64=0.0, aic::Float64=0.0, resid::Vector{Float64}=Array(Float64, 0), optimCode::Int=0, rModelName::ASCIIString="")
#METHODS
#	copy(a::ARIMAModel): Standard copy
#	deepcopy(a::ARIMAModel): Standard deepcopy
#	show(io::IO, a::ARIMAModel): Standard show
#	show(a::ARIMAModel): Show defaulting to STDOUT
#	coef(a::ARIMAModel): Returns all coefficients (and intercept if applicable) in one vector ordered [arCoef, maCoef, sarCoef, smaCoef, intecept (if applicable)]
#	vcov(a::ARIMAModel): Return the variance covariance matrix of coefficients
#	coefVar(a::ARIMAModel): Return a 2x? matrix containing all coefficients in first row (ie output of coef functions) and the corresponding variances of these estimators in the second row.
#NOTES
#	Typically this type is created when the user calls the estimate function on an ARIMAInput or AutoARIMAInput. However, the user may also want to create a specific ARIMAModel to feed into the simulate function.
#	As well as serving as output from estimate step, this type serves as input to simulate function, although obviously not all the fields will be relevant.
#----------------------------------------------------------
type ARIMAModel
	arCoef::Vector{Float64}
	maCoef::Vector{Float64}
	sarCoef::Vector{Float64}
	smaCoef::Vector{Float64}
	intercept::Float64
	sigma2::Float64
	covCoef::Matrix{Float64}
	covCoefName::Vector{ASCIIString}
	d::Int
	D::Int
	seasonalPeriod::Int
	logLik::Float64
	aic::Float64
	resid::Vector{Float64}
	optimCode::Int
	rModelName::ASCIIString
	function ARIMAModel(arCoef::Vector{Float64}, maCoef::Vector{Float64}, sarCoef::Vector{Float64}, smaCoef::Vector{Float64}, intercept::Float64, sigma2::Float64, covCoef::Matrix{Float64}, covCoefName::Vector{ASCIIString}, d::Int, D::Int, seasonalPeriod::Int, logLik::Float64, aic::Float64, resid::Vector{Float64}, optimCode::Int, rModelName::ASCIIString)
		(d < 0 || D < 0) && error("d and D must be non-negative")
		sigma2 < 0 && error("sigma2 must be non-negative")
		size(covCoef, 1) != size(covCoef, 2) && error("Coefficient covariance matrix is non-square")
		numCoefTemp = length(arCoef) + length(maCoef) + length(sarCoef) + length(smaCoef) + convert(Int, !(isnan(intercept)))
		size(covCoef, 1) != numCoefTemp && error("Logic fail. Number of coefficients is inconsistent with size of coefficient covariance matrix")
		length(covCoefName) != numCoefTemp && error("Logic fail. Number of coefficients is inconsistent with number of coefficient names")
		seasonalPeriod < 1 && error("seasonalPeriod must be strictly positive")
		funcOut = new(arCoef, maCoef, sarCoef, smaCoef, intercept, sigma2, covCoef, covCoefName, d, D, seasonalPeriod, logLik, aic, resid, optimCode, rModelName)
		return(funcOut)
	end
end


#--------------------
#Constructors
#--------------------
ARIMAModel() = ARIMAModel(Array(Float64, 0), Array(Float64, 0), Array(Float64, 0), Array(Float64, 0), 0.0, 1.0, Array(Float64, 0, 0), Array(ASCIIString, 0), 0, 0, 1, 0.0, 0.0, Array(Float64, 0), 0, "")
#ARIMAInput constructor for converting an ARIMAModel to ARIMAInput (necessary for code to be here after ARIMAModel type definition)
ARIMAInput(a::ARIMAModel) = ARIMAInput(length(a.arCoef), a.d, length(a.maCoef), length(a.sarCoef), a.D, length(a.smaCoef), a.seasonalPeriod, !(isnan(a.intercept)), true, Array(Float64, 0), Array(Float64, 0), "CSS-ML")
#--------------------
#Keyword constructor
#--------------------
function ARIMAModel(; arCoef::Vector{Float64}=Array(Float64, 0), maCoef::Vector{Float64}=Array(Float64, 0), sarCoef::Vector{Float64}=Array(Float64, 0), smaCoef::Vector{Float64}=Array(Float64, 0), intercept::Float64=0.0, sigma2::Float64=1.0, covCoef::Matrix{Float64}=Array(Float64, 0, 0), covCoefName::Vector{ASCIIString}=Array(ASCIIString, 0), d::Int=0, D::Int=0, seasonalPeriod::Int=1, logLik::Float64=0.0, aic::Float64=0.0, resid::Vector{Float64}=Array(Float64, 0), optimCode::Int=0, rModelName::ASCIIString="")
	numCoefTemp = length(arCoef) + length(maCoef) + length(sarCoef) + length(smaCoef) + convert(Int, !(isnan(intercept)))
	if size(covCoef, 1) == 0 && size(covCoef, 2) == 0
		covCoef = zeros(Float64, numCoefTemp, numCoefTemp)
	end
	if length(covCoefName) == 0
		covCoefName = Array(ASCIIString, numCoefTemp)
	end
	return(ARIMAModel(arCoef, maCoef, sarCoef, smaCoef, intercept, sigma2, covCoef, covCoefName, d, D, seasonalPeriod, logLik, aic, resid, optimCode, rModelName))
end
#--------------------
#Methods
#--------------------
#copy, deepcopy
copy(a::ARIMAModel) = ARIMAModel(copy(a.arCoef), copy(a.maCoef), copy(a.sarCoef), copy(a.smaCoef), copy(a.intercept), copy(a.sigma2), copy(a.covCoef), copy(a.d), copy(a.D), copy(a.seasonalPeriod), copy(a.logLik), copy(a.aic), copy(a.resid), copy(a.optimCode), copy(a.rModelName))
deepcopy(a::ARIMAModel) = ARIMAModel(deepcopy(a.arCoef), deepcopy(a.maCoef), deepcopy(a.sarCoef), deepcopy(a.smaCoef), deepcopy(a.intercept), deepcopy(a.sigma2), deepcopy(a.covCoef), deepcopy(a.d), deepcopy(a.D), deepcopy(a.seasonalPeriod), deepcopy(a.logLik), deepcopy(a.aic), deepcopy(a.resid), deepcopy(a.optimCode), deepcopy(a.rModelName))
#show
function show(io::IO, a::ARIMAModel)
	println("ARIMAModel contents:")
	println("    AR coefficients = " * string(a.arCoef))
	println("    MA coefficients = " * string(a.maCoef))
	println("    SAR coefficients = " * string(a.sarCoef))
	println("    SMA coefficients = " * string(a.smaCoef))
	println("    intercept = " * string(a.intercept))
	println("    residual variance = " * string(a.sigma2))
	println("    order of integration = " * string(a.d))
	println("    seasonal order of integration = " * string(a.D))
	println("    seasonal period = " * string(a.seasonalPeriod))
	println("    log-likelihood at maximum = " * string(a.logLik))
	println("    aic at maximum = " * string(a.aic))
	println("    optimisation output code = " * string(a.optimCode))
	println("    R model name = " * string(a.rModelName))
end
show(a::ARIMAModel) = show(STDOUT, a)
#coef
coef(a::ARIMAModel) = isnan(a.intercept) ? [a.arCoef, a.maCoef, a.sarCoef, a.smaCoef] : [a.arCoef, a.maCoef, a.sarCoef, a.smaCoef, [a.intercept]]
#vcov (to match R function)
vcov(a::ARIMAModel) = a.covCoef
#coefVar (get coefficients and their variances)
coefVar(a::ARIMAModel) = cat(2, coef(a), diag(a.covCoef))'




#----------------------------------------------------------
#FUNCTION
#	estimate
#	estimateARIMA
#	estimateAutoARIMA
#METHODS
#	{T<:Number}(x::Vector{T}, a::ARIMAInput, method::ASCIIString): Estimate an ARIMAModel using data x and estimation inputs in ARIMAInput. Do estimation using R arima function with method determined by method ASCIIString.
#	*estimate{T<:Number}(x::Vector{T}, a::ARIMAInput, method::Symbol): Estimate an ARIMAModel using data x and estimation inputs in ARIMAInput. Do estimation using Julia Optim package. WARNING: NOT YET IMPLEMENTED.
# 	estimate{T<:Number}(x::Vector{T}, a::ARIMAInput): Wrapper obviating the need to specify method as separate argument (since method is part of ARIMAInput)
# 	estimate{T<:Number}(a::ARIMAInput, x::Vector{T}): Wrapper ensuring order of input doesn't matter.
#	estimate{T<:Number}(x::Vector{T}, a::AutoARIMAInput): Automatically choose order and lags of ARIMAModel using inputs in a::AutoARIMAInput and data in x. This method uses the R auto.arima function from Rob Hyndman's forecast library.
#	estimate{T<:Number}(a::AutoARIMAInput, x::Vector{T}): Wrapper ensuring order of input doesn't matter.
#KEYWORD METHODS
#	estimateARIMA{T<:Number}(x::Vector{T}; p::Int=0, d::Int=0, q::Int=0, P::Int=0, D::Int=0, Q::Int=0, seasonalPeriod::Int=1, includeIntercept::Bool=true, transformParams::Bool=true, fixedParams::Vector{Float64}=Array(Float64, 0), initParams::Vector{Float64}=Array(Float64, 0), method::ASCIIString="CSS-ML")
#	estimateAutoARIMA{T<:Number}(x::Vector{T}; d::Int=-1, D::Int=-1, pMax::Int=5, qMax::Int=5, PMax::Int=2, QMax::Int=2, numCoefMax::Int=5, dMax::Int=2, DMax::Int=1, pStart::Int=2, qStart::Int=2, PStart::Int=1, QStart::Int=5, stationary::Bool=false, seasonal::Bool=false, ic::ASCIIString="aic", stepwise::Bool=true, parallel::Bool=false, numCore::Int=2)
#OUTPUT
#	ARIMAModel containing the estimated ARIMA coefficients and other variables of interest.
#PURPOSE
#	The purpose of this function is to estimate the parameters of the input model and return a complete specification of the model
#NOTES
#	If I write any custom optimization routines, make them their own type, and then add a method to estimate
#	* implies this method has not yet been implemented
#NON-EXPORTED FUNCTIONS
#	estimate_RInputString(a::ARIMAInput): Used to build the string that gets used to call the arima function in embedded R
#	estimate_RInputStringAuto(a::AutoARIMAInput): Used to build the string that gets used to call the auto.arima function in embedded R
#	estimate_RextractCoef(allCoef::Vector{Float64}, allCoefNames::Vector{ASCIIString}): Used to split up coefficient output from R into AR, MA, SAR, SMA, and intercept.
#----------------------------------------------------------
#estimate method that uses arima functionality in R (distinguished via method being an ASCIIString)
function estimate{T<:Number}(x::Vector{T}, a::ARIMAInput, method::ASCIIString)
	globalEnv[:x] = x #Send data to R
	reval("x <- Re(x)") #Convert from complex back to real (bug with sending Float64 to R)
	rModelName = "m_" * randstring(10) #Randomly allocate a model name
	if rcopy("exists(\"" * rModelName * "\")") == 1
		rModelName = "m_" * randstring(10) #If model name already exists, allocate a new one
		if rcopy("exists(\"" * rModelName * "\")") == 1 #The probability of failing to generate a unique name with two attempts is essentially zero, so throw an error.
			error("Logic fail. Unable to generate unique R model name.")
		end
	end
	try
		reval(rModelName * " <- " * estimate_RInputString(a)) #Estimate ARIMA model and store output in embedded R using rModelName
	catch
		println("WARNING: R Routine failed to estimate ARIMA model. Returning white noise ARIMAModel")
		return(ARIMAModel(optimCode=-1))
	end
	allCoef = rcopy(rModelName * "\$coef")
	if length(allCoef) > 0
		allCoefName = rcopy("names(" * rModelName * "\$coef" * ")")
	else
		allCoefName = Array(ASCIIString, 0)
	end
	(arCoef, maCoef, sarCoef, smaCoef, intercept) = estimate_RextractCoef(allCoef, allCoefName)
	sigma2 = rcopy(rModelName * "\$sigma2")[1]
	if length(allCoef) > 0
		covCoef = rcopy(rModelName * "\$var.coef")
	else
		covCoef = Array(Float64, 0, 0)
	end
	if ndims(covCoef) == 1
		covCoef = covCoef''
	end
	logLik = rcopy(rModelName * "\$loglik")[1]
	if method == "ML" || method == "CSS-ML"
		aicTemp = rcopy(rModelName * "\$aic")[1]
	else
		(aicTemp = NaN)
	end
	resid = rcopy(rModelName * "\$residuals")
	optimCode = convert(Int, rcopy(rModelName * "\$code")[1])
	optimCode != 0 && println("WARNING: Non-zero exit code in estimate function on ARIMAInput. Code = " * string(optimCode) * ". Refer to R docs for more detail.")
	return(ARIMAModel(arCoef, maCoef, sarCoef, smaCoef, intercept, sigma2, covCoef, allCoefName, a.d, a.D, a.seasonalPeriod, logLik, aicTemp, resid, optimCode, rModelName))
end
#estimate method that uses Julia Optim package (distinguished via method being a Symbol)
function estimate{T<:Number}(x::Vector{T}, a::ARIMAInput, method::Symbol)
	error("I have not coded up this method of estimate yet")
end
#POTENTIALLY include another method here that allows optimisation via minimisation of arbitrary function of model residuals
#estimate wrappers (ARIMAInput)
estimate{T<:Number}(x::Vector{T}, a::ARIMAInput) = estimate(x, a, a.method)
estimate{T<:Number}(a::ARIMAInput, x::Vector{T}) = estimate(x, a)
#Keyword versions
estimateARIMA{T<:Number}(x::Vector{T}; p::Int=0, d::Int=0, q::Int=0, P::Int=0, D::Int=0, Q::Int=0, seasonalPeriod::Int=1, includeIntercept::Bool=true, transformParams::Bool=true, fixedParams::Vector{Float64}=Array(Float64, 0), initParams::Vector{Float64}=Array(Float64, 0), method::ASCIIString="CSS-ML") = estimate(x, ARIMAInput(p=p, d=d, q=q, P=P, D=D, Q=Q, seasonalPeriod=seasonalPeriod, includeIntercept=includeIntercept, transformParams=transformParams, fixedParams=fixedParams, initParams=initParams, method=method))
estimateAutoARIMA{T<:Number}(x::Vector{T}; d::Int=-1, D::Int=-1, pMax::Int=5, qMax::Int=5, PMax::Int=2, QMax::Int=2, numCoefMax::Int=5, dMax::Int=2, DMax::Int=1, pStart::Int=2, qStart::Int=2, PStart::Int=1, QStart::Int=5, stationary::Bool=false, seasonal::Bool=false, ic::ASCIIString="aic", stepwise::Bool=true, parallel::Bool=false, numCore::Int=2) = estimate(x, AutoARIMAInput(d=d, D=D, pMax=pMax, qMax=qMax, PMax=PMax, QMax=QMax, numCoefMax=numCoefMax, dMax=dMax, DMax=DMax, pStart=pStart, qStart=qStart, PStart=PStart, QStart=QStart, stationary=stationary, seasonal=seasonal, ic=ic, stepwise=stepwise, parallel=parallel, numCore=numCore))
#Estimate ARIMA model with automatic order detection
function estimate{T<:Number}(x::Vector{T}, a::AutoARIMAInput)
	globalEnv[:x] = x #Send data to R
	reval("x <- Re(x)") #Convert from complex back to real (bug with sending Float64 to R)
	rModelName = "m_" * randstring(10) #Randomly allocate a model name
	if rcopy("exists(\"" * rModelName * "\")") == 1
		rModelName = "m_" * randstring(10) #If model name already exists, allocate a new one
		if rcopy("exists(\"" * rModelName * "\")") == 1 #The probability of failing to generate a unique name with two attempts is essentially zero, so throw an error.
			error("Logic fail. Unable to generate unique R model name.")
		end
	end
	reval(rModelName * " <- " * estimate_RInputStringAuto(a)) #automatically choose and estimate ARIMA model and store output in embedded R using rModelName
	allCoef = rcopy(rModelName * "\$coef")

	if length(allCoef) > 0
		allCoefName = rcopy("names(" * rModelName * "\$coef" * ")")
	else
		allCoefName = Array(ASCIIString, 0)
	end
	(arCoef, maCoef, sarCoef, smaCoef, intercept) = estimate_RextractCoef(allCoef, allCoefName)
	sigma2 = rcopy(rModelName * "\$sigma2")[1]
	if length(allCoef) > 0
		covCoef = rcopy(rModelName * "\$var.coef")
	else
		covCoef = Array(Float64, 0, 0)
	end
	if ndims(covCoef) == 1
		covCoef = covCoef''
	end
	armaRep = rcopy(rModelName * "\$arma")
	seasonalPeriod = convert(Int, armaRep[5])
	d = convert(Int, armaRep[6])
	D = convert(Int, armaRep[7])
	logLik = rcopy(rModelName * "\$loglik")[1]
	aicTemp = rcopy(rModelName * "\$aic")[1]
	resid = rcopy(rModelName * "\$residuals")
	optimCode = convert(Int, rcopy(rModelName * "\$code")[1])
	optimCode != 0 && println("WARNING: Non-zero optimisation exit code in R arima function. Code = " * string(optimCode) * ". Refer to R docs for more detail.")
	return(ARIMAModel(arCoef, maCoef, sarCoef, smaCoef, intercept, sigma2, covCoef, allCoefName, d, D, seasonalPeriod, logLik, aicTemp, resid, optimCode, rModelName))
end
estimate{T<:Number}(a::AutoARIMAInput, x::Vector{T}) = estimate(x, a)
#Non-exported function called by estimate on ARIMAInput
function estimate_RInputString(a::ARIMAInput)
	s1 = "arima(x, "
	s2 = "order=c(" * string(a.p) * "," * string(a.d) * "," * string(a.q) * "), "
	s3 = "seasonal=list(order=c(" * string(a.P) * "," * string(a.D) * "," * string(a.Q) * "), period=" * string(a.seasonalPeriod) * "), "
	if a.includeIntercept == false || a.d > 0
		s4 = "include.mean=FALSE, "
	else
		s4 = "include.mean=TRUE, "
	end
	s5 = "transform.pars=" * uppercase(string(a.transformParams)) * ", "
	if length(a.fixedParams) > 0
		s6 = "fixed=c(" * string(a.fixedParams)[2:end-1] * "), "
		s6 = replace(s6, "NaN", "NA")
	else
		s6 = "fixed=NULL, "
	end
	if length(a.initParams) > 0
		s7 = "init=c(" * string(a.initParams)[2:end-1] * "), "
		s7 = replace(s7, "NaN", "0")
	else
		s7 = "init=NULL, "
	end
	s8 = "method=\"" * a.method * "\")"
	return(s1 * s2 * s3 * s4 * s5 * s6 * s7 * s8)
end
#Non-exported function called by estimate on AutoARIMAInput
function estimate_RInputStringAuto(a::AutoARIMAInput)
	s1 = "auto.arima(x, "
	a.d < 0 ? (s2 = "d=NA, ") : (s2 = "d=" * string(a.d) * ", ")
	a.D < 0 ? (s3 = "D=NA, ") : (s3 = "D=" * string(a.D) * ", ")
	s4 = "max.p=" * string(a.pMax) * ", max.q=" * string(a.qMax) * ", "
	s5 = "max.P=" * string(a.PMax) * ", max.Q=" * string(a.QMax) * ", "
	s6 = "max.order=" * string(a.numCoefMax) * ", "
	s7 = "max.d=" * string(a.dMax) * ", max.D=" * string(a.DMax) * ", "
	s8 = "start.p=" * string(a.pStart) * ", start.q=" * string(a.qStart) * ", "
	s9 = "start.P=" * string(a.PStart) * ", start.Q=" * string(a.QStart) * ", "
	(a.stationary == true) ? (s10 = "stationary=TRUE, ") : (s10 = "stationary=FALSE, ")
	(a.seasonal == true) ? (s11 = "seasonal=TRUE, ") : (s11 = "seasonal=FALSE, ")
	s12 = "ic=\"" * a.ic * "\", "
	(a.stepwise == true) ? (s13 = "stepwise=TRUE, ") : (s13 = "stepwise=FALSE, ")
	(a.parallel == true) ? (s14 = "parallel=TRUE, ") : (s14 = "parallel=FALSE, ")
	s15 = "num.cores=" * string(a.numCore) * ", "
	s16 = "allowdrift=FALSE)" #Currently not allowed since other functions in this module do not use forecast package
	return(s1 * s2 * s3 * s4 * s5 * s6 * s7 * s8 * s9 * s10 * s11 * s12 * s13 * s14 * s15 * s16)
end
#Non-exported function to extract coefficients from R coefficient vector and name vector
function estimate_RextractCoef(allCoef::Vector{Float64}, allCoefNames::Vector{ASCIIString})
	(length(allCoef) != length(allCoefNames)) && error("number of coefficients does not match number of names")
	arCoef = Array(Float64, 0)
	maCoef = Array(Float64, 0)
	sarCoef = Array(Float64, 0)
	smaCoef = Array(Float64, 0)
	intercept = NaN
	for n = 1:length(allCoef)
		if allCoefNames[n][1:2] == "ar"
			push!(arCoef, allCoef[n])
		elseif allCoefNames[n][1:2] == "ma"
			push!(maCoef, allCoef[n])
		elseif allCoefNames[n][1:2] == "sa"
			push!(sarCoef, allCoef[n])
		elseif allCoefNames[n][1:2] == "sm"
			push!(smaCoef, allCoef[n])
		elseif allCoefNames[n][1:2] == "in"
			intercept = allCoef[n]
		else
			error("Unidentified name in coefficient name vector")
		end
	end
	return(arCoef, maCoef, sarCoef, smaCoef, intercept)
end






#----------------------------------------------------------
#FUNCTION
#	forecast
#INPUT
#	forecast(a::ARIMAModel, numStep::Int=1): Forecast numStep ahead using ARIMAModel. Note, this method does not input the data because it relies on finding the correct ARIMA model in embedded R workspace using a.rModelName and then calling predict.Arima on this model. This method is needed because the method written entirely in Julia cannot handle bounds computations.
#	forecast{T<:Number}(x::Vector{T}, a::ARIMAModel, numStep::Int=1): Forecast numStep ahead using ARIMAModel and data in x. This method is implemented entirely in Julia, but is not able to handle seasonal ARIMA models or bounds computations.
#	forecast{T<:Number}(a::ARIMAModel, x::Vector{T}, numStep::Int=1): Wrapper so that order of inputs doesn't matter.
#OUTPUT
#	Output is tuple (pointFore::Vector{Float64}, boundFore::Matrix{Float64}), where pointFore provides the point forecast up to numStep ahead, and boundFore provides 95% lower and upper bounds for the forecasts (lower bound in first column, upper bound in second column)
#PURPOSE
#	The purpose of this function is to provide forecasts for the input dataset and ARIMA model
#NOTES
#	When d > 0 (order of integration greater than 0), the methods that do not wrap R functions use the non-exported function arima_statdiff, arima_statdiffinv, arima_statdiffinv_one. These functions difference or inverse difference data up to the specified order.
#NON-EXPORTED FUNCTIONS
#	forecast_arma{T<:Number}(x::Vector{T}, resid::Vector{T}, arCoef::Vector{Float64}, maCoef::Vector{Float64}, intercept::Float64=NaN, numStep::Int=1): This function provides an ARMA point forecast implemented entirely in Julia.
#----------------------------------------------------------
#Method for using predict.Arima function in R
function forecast(a::ARIMAModel, numStep::Int=1)
	if rcopy("exists(\"" * a.rModelName * "\")")[1] != 1
		error("Unable to find " * a.rModelName * " in embedded R workspace. Try using the function method that includes the underlying data.")
	else
		rFore = rcopy("predict(" * a.rModelName * ", n.ahead=" * string(numStep) * ")")
		pointFore = rFore[1]
		boundFore = [[pointFore - 1.96*rFore[2]] [pointFore + 1.96*rFore[2]]]
	end
	return(pointFore, boundFore)
end
#Method implemented entirely in Julia (note, bounds computations not yet implemented)
function forecast{T<:Number}(x::Vector{T}, a::ARIMAModel, numStep::Int=1)
	if length(a.sarCoef) > 0 || length(a.smaCoef) > 0 || a.D > 0
		error("This method does not currently support seasonal ARIMA models. Try using the function method that does not include underlying data (ARIMA model must be available in embedded R workspace).")
	end
	if a.d == 0
		(pointFore, boundFore) = forecast_arma(x, a.resid, a.arCoef, a.maCoef, a.intercept, numStep)
	else
		maxLag = max(length(a.arCoef), length(a.maCoef)) + a.d
		xS = x[end-maxLag:end]
		(xD, init) = arima_statdiff(xS, a.d)
		(pointFore, boundFore) = forecast_arma(xD, a.resid[end-length(xD)+1:end], a.arCoef, a.maCoef, a.intercept, numStep)
		pointFore = arima_statdiffinv([xD, pointFore], init)
		pointFore = pointFore[end-numStep+1:end]
		#bound calcuations will need to be added here if they are ever implemented in forecast_arma
	end
	return(pointFore, boundFore)
end
#Wrapper to make order of input immaterial
forecast{T<:Number}(a::ARIMAModel, x::Vector{T}, numStep::Int=1) = forecast(x, a, numStep)
#forecast_arma
function forecast_arma{T<:Number}(x::Vector{T}, resid::Vector{T}, arCoef::Vector{Float64}, maCoef::Vector{Float64}, intercept::Float64=NaN, numStep::Int=1)
	length(x) != length(resid) && error("Input x vector and residual vector must be the same length")
	!(isnan(intercept)) && (x = x - intercept) #de-mean data if we have an intercept
	maxLag = max(length(arCoef), length(maCoef))
	eFore = [resid[end-maxLag+1:end], zeros(Float64, numStep)]
	xFore = [x[end-maxLag+1:end], zeros(Float64, numStep)]
	for k = maxLag+1:maxLag+numStep
		for j = 1:length(arCoef)
			xFore[k] += arCoef[j] * xFore[k-j]
		end
		for j = 1:length(maCoef)
			xFore[k] += maCoef[j] * eFore[k-j]
		end
	end
	pointFore = xFore[end-numStep+1:end]
	!(isnan(intercept)) && (pointFore = pointFore + intercept) #add mean back in if it exists
	boundFore = NaN * Array(Float64, numStep, 2) #bound calculations could be added here
	return(pointFore, boundFore)
end





#----------------------------------------------------------
#FUNCTION
#	simulate
#	simulateARIMA
#METHODS
#	simulate(numObs::Int, a::ARIMAModel, numBurn::Int=-1, resid::Vector{Float64}=[NaN]): Simulate numObs observations using ARIMAModel a, with burn-in length numBurn and optionally use resid as the source of random numbers in simulation.
#	simulate(a::ARIMAModel, numObs::Int, numBurn::Int=-1, resid::Vector{Float64}=[NaN]): Wrapper so order doesn't matter.
#KEYWORD METHOD
#	simulateARIMA(numObs::Int; arCoef::Vector{Float64}=Array(Float64, 0), maCoef::Vector{Float64}=Array(Float64, 0), d::Int=0, intercept::Float64=NaN, sigma2::Float64=1.0, numBurn::Int=-1, resid::Vector{Float64}=[NaN]): Keyword activated simulate method. First input specifiying number of observations to simulate is non-optional. Default keywords result in white noise model.
#OUTPUT
#	Output is Vector{Float64} of data simulated using the specified ARIMA model.
#PURPOSE
#	The purpose of this function is to simulate the specified ARIMA model
#NOTES
#	This function is not currently capable of simulating seasonal ARIMA models, ie non-seasonal only.
#NON-EXPORTED FUNCTION
#	simulate_arma(numObs::Int, arCoef::Vector{Float64}, maCoef::Vector{Float64}, intercept::Float64=0, sigma2::Float64=1, numBurn::Int=-1, resid::Vector{Float64}=[NaN]): Simulates an ARMA model using specified inputs. simulate calls this function for simulating the ARMA part and uses arima_statdiff and arima_statdiffinv to take care of the order of integration
#----------------------------------------------------------
function simulate(numObs::Int, a::ARIMAModel, numBurn::Int=-1, resid::Vector{Float64}=[NaN])
	if length(a.sarCoef) > 0 || length(a.smaCoef) > 0 || a.D > 0
		error("Simluation of seasonal ARIMA models has not been implemented yet")
	end
	isnan(a.intercept) ? intercept = 0.0 : intercept = a.intercept
	xSim = simulate_arma(numObs, a.arCoef, a.maCoef, intercept, a.sigma2, numBurn, resid)
	if a.d > 0
		xSim = arima_statdiffinv(xSim, zeros(Float64, a.d))[a.d+1:end]
	end
	return(xSim)
end
#Wrapper so order doesn't matter
simulate(a::ARIMAModel, numObs::Int, numBurn::Int=-1, resid::Vector{Float64}=[NaN]) = simulate(numObs, a, numBurn, resid)
#Keyword activated function
simulateARIMA(numObs::Int; arCoef::Vector{Float64}=Array(Float64, 0), maCoef::Vector{Float64}=Array(Float64, 0), d::Int=0, intercept::Float64=NaN, sigma2::Float64=1.0, numBurn::Int=-1, resid::Vector{Float64}=[NaN]) = simulate(numObs, ARIMAModel(arCoef=arCoef, maCoef=maCoef, d=d, intercept=intercept, sigma2=sigma2), numBurn, resid)
#non-exported function that simulates an arma
function simulate_arma(numObs::Int, arCoef::Vector{Float64}, maCoef::Vector{Float64}, intercept::Float64=0, sigma2::Float64=1, numBurn::Int=-1, resid::Vector{Float64}=[NaN])
	p = length(arCoef)
	q = length(maCoef)
	if numBurn < 0
		#IMPLEMENT HERE: A better method for choosing the number of observations to burn. Should be straightforward with Polynomials package, e.g. p+q+ceil(6/log(minroots)), where minroots = min(modulus(roots(p))) where p is polynomial [-1, ar coefficients]
		numBurn = p + q + 200
	elseif numBurn < p + q
		error("Input number of burn-in observations must be of at least length p + q")
	end
	length(resid) == 0 && error("Invalid optional input resid")
	isnan(resid[1]) && (resid = sqrt(sigma2) * randn(numObs + numBurn))
	if length(resid) > numObs + numBurn
		#resid = resid[end-numObs-numBurn+1:end]
		resid = resid[1:numObs+numBurn]
	elseif length(resid) < numObs + numBurn
		error("Number of supplied residuals must be >= numObs + numBurn")
	end
	xSim = intercept * ones(Float64, numObs + numBurn)
	for n = max(p, q)+1:length(xSim)
		for k = 1:length(arCoef) #
			xSim[n] += arCoef[k] * xSim[n-k]
		end
		for k = 1:length(maCoef)
			xSim[n] += maCoef[k] * resid[n-k]
		end
		xSim[n] += resid[n]
	end
	return(xSim[end-numObs+1:end])
end







#----------------------------------------------------------
#FUNCTION
#	arima_statdiff
#	arima_statdiffinv_one
#	arima_statdiffinv
#INPUT
#	arima_statdiff{T<:Number}(x::Vector{T}, order::Int=1): x is a vector of numbers to discretely differentiate. order is the order of differentiation. If order=1, this is equivalent to calling diff from Juia base.
#	arima_statdiffinv_one{T<:Number}(x::Vector{T}, init::T): x is a vector of numbers to discretely integrate ONCE. init is the initial value needed for discrete integration
#	arima_statdiffinv{T<:Number}(x::Vector{T}, init::Vector{T}): x is a vector of numbers to discretely integrate length(init) times. init provides the initial values needed for each successive discrete integration, where the first integration of x uses init[1] and the last integration uses init[end]
#OUTPUT
#	arima_statdiff outputs the tuple (xOut, init), where xOut is the differentiated series and init is the first value of each series prior to the last differentiated series, so init[1] = x[1], init[2] = (delta x)[1], ..., and init[end] = (delta^(order-1) x)[1]
#	arima_statdiffinv_one and arima_statdiffinv both output Vector{Float64} storing the discretely integrated series.
#PURPOSE
#	These three function provide discrete differentiation and discrete integration.
#NOTES
#	These functions are not exported.
#	More efficient versions of arima_statdiff and arima_statdiffinv will almost certainly eventually become a part of Julia Base. The functions provided here can be used in the meantime.
#----------------------------------------------------------
#arima_statdiff
function arima_statdiff{T<:Number}(x::Vector{T}, order::Int=1)
	#WARNING: Temporary array allocated for each order of differencing. Highly inefficient when order is large.
	order < 1 && error("order of statistical differencing must be at least 1")
	length(x) < order + 1 && error("input vector x too short for specified order of differencing")
	xOut = diff(x)
	init = Array(T, order)
	init[1] = x[1]
	for n = 2:order
		init[n] = xOut[1]
		xOut = diff(xOut)
	end
	return(xOut, init)
end
#arima_statdiffinv_one
function arima_statdiffinv_one{T<:Number}(x::Vector{T}, init::T)
	length(x) < 1 && error("input vector x must contain at least one element")
	xOut = Array(T, length(x) + 1)
	xOut[1] = init
	for n = 2:length(xOut)
		xOut[n] = xOut[n-1] + x[n-1]
	end
	return(xOut)
end
#arima_statdiffinv
function arima_statdiffinv{T<:Number}(x::Vector{T}, init::Vector{T})
	#WARNING: Temporary array allocated for each order of integration. Highly inefficient when order is large.
	length(init) < 1 && error("order of statistical differencing must be at least 1")
	xOut = arima_statdiffinv_one(x, init[end])
	for n = 2:length(init)
		xOut = arima_statdiffinv_one(xOut, init[end-n+1])
	end
	return(xOut)
end








end # module
