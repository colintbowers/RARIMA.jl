## NOTE

This package is no longer actively maintained and will not work on the latest versions of Julia.

## RARIMA

An [ARIMA](http://en.wikipedia.org/wiki/Autoregressive_integrated_moving_average) (Autoregressive Integrated Moving Average) module for the [Julia](http://julialang.org/) language that uses libraries from the [R language](http://www.r-project.org/)


## Main Features

This module allows Julia users to estimate, forecast, and simulate ARIMA time-series models.

Some functions are implemented in pure Julia, while others call functions from the R language. The relevant R functions are accessed using the Julia package [RCall](https://github.com/JuliaStats/RCall.jl). These R functions come from either the R [stats](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/00Index.html) library or the R [forecast](http://cran.r-project.org/web/packages/forecast/index.html) library. Since these R functions are available via the GNU General Public License (v2 or v3), this module uses the GNU GPL v2.


## Setting up RCall

It is recommended that you set up `RCall` and verify it works prior to installing `RARIMA`. `RCall` is an officially registered package, so it can be installed from the Julia REPL using `Pkg.add("RCall")`. For more detail, please refer to the [`RCall` docs](https://github.com/JuliaStats/RCall.jl).

You can verify `RCall` is working correctly using the following commands at the Julia REPL:

    using RCall
    rcopy("1+2")

which will hopefully return the obvious answer. Also note, you *must* install the `forecast` library into your R library before using `RARIMA`. This can be done from an R terminal (or R IDE of choice, e.g. R-Studio) using the command:

    install.packages("forecast")

to install the library, then:

    library(forecast)

to ensure it loads.


## A Word on Julia versus R

The long-term goal is to implement ARIMA time-series modelling capability purely in Julia (preferably under an MIT license rather than GPL). This module exists only to fill the gap between the present day and the day that future goal is realised. I've tried to structure `RARIMA` so that it can serve as the basis of that future module, ie it utilises Julia's type system and multiple dispatch, and have implemented some of the functionality already in Julia.

The most significant functionality that needs to be converted to Julia is the estimation routines (regular and automatic order selection). Currently, these are implemented entirely in R. Because the estimation step typically involves solving a numerical optimisation problem, sometimes things go wrong. One of the main challenges in writing this module is dealing with all the different types of errors that can result from the estimation step in R and passing something informative back to Julia. I haven't spent a huge amount of time on this since the long-term goal is to replace that code with Julia code. Currently, the type returned from estimation functions contains a field `optimCode` that will be set to `0` if everything worked correctly. If there was a known problem from the R end, this field will take on a positive value, and the corresponding optimisation code can be looked up in the R documentation. If something went wrong that the Julia code is not able to understand, this field will be set to `-1`. Watch out for this!


## How to use RARIMA

To install `RARIMA`, first ensure `RCall` is installed and working, and that your R library contains the library `forecast`. See the above section "Setting up `RCall`" for more detail.

Since `RARIMA` is not currently an official Julia package, you can install it using:

    Pkg.clone("https://github.com/colintbowers/RARIMA.jl")

at the Julia REPL. If this is successful, you can now load the package using:

    using(RARIMA)

The nuts-and-bolts method for working with this module is via three types, `ARIMAInput`, `AutoARIMAInput`, and `ARIMAModel`, and three generic functions `estimate`, `forecast`, and `simulate`. However, keyword-based functionality is available, and this is probably how most users will choose to interact with this package. The keywords are based on the fields of the above three types. A short-list of the most important fields and the associated functionality is provided now.

#### ARIMAInput / Estimation

The most important fields in the type `ARIMAInput` are:

* p::Int=0 (autoregressive order) 
* d::Int=0 (order of integration)
* q::Int=0 (moving average order)
* P::Int=0 (seasonal autoregressive order)
* D::Int=0 (order of seasonal integration)
* Q::Int=0 (seasonal moving average order)
* seasonalPeriod::Int=1 (seasonal period - only relevant if P, D, or Q is > 0)
* includeIntercept::Bool=true (include an intercept term in the model - only relevant if d and D equal 0)

So one could construct an `ARIMAInput` for an ARIMA(2,1,0)(1,0,0) model of monthly data with no intercept using:

    aI = ARIMAInput(p=2, d=1, P=1, seasonalPeriod=12, includeIntercept=false)

Actually, in this case we could have omitted `includeIntercept` because it is automatically set to false when `d > 0` (see the R docs for `arima` for more detail).

Given a vector of observed data, `x::Vector{T}` where `T<:Number`, we can estimate the parameters of a model specified by and `ARIMAInput` using:

    aM = estimate(x, aI)

where the output is of type `ARIMAModel` (discussed later in this section). Note, `estimate(aI, x)` is equivalent (ie order doesn't matter). Alternatively, one can make the construction of the `ARIMAInput` implicit by calling the function `estimateARIMA` using keywords. For example, an equivalent command to the one above would be:

    aM = estimateARIMA(x, p=2, d=1, P=1, seasonalPeriod=12, includeIntercept=false)

#### AutoARIMAInput / Estimation

In most situations, the user doesn't actually know what order of ARIMA model best suits the data. For this situation, one can use the type `AutoARIMAInput`. The most important fields of this type are:

* d::Int=-1 (used to fix the order of integration of the ARIMA model. A negative value implies you do not wish to fix the order of integration)
* D::Int=-1 (used to fix the order of seasonal integration of the ARIMA model. A negative value implies you do not wish to fix the seasonal order of integration)
* seasonal::Bool=false (Set to false to restrict search procedure to non-seasonal ARIMA models)

Most users will be happy to leave all other fields set to their default values. So, for example, a user who wishes to estimate an ARIMA model (given data vector `x`) that they *suspect* might be seasonal and they *know* is integrated of order 1 would use:

    aAI = AutoARIMAInput(d=1, seasonal=true)

to construct the `AutoARIMAInput`, and then:

    aM = estimate(x, aAI)

to estimate the order of the model and its parameters. As before, the output is of type `ARIMAModel`. Alternatively, the construction of the `AutoARIMAInput` can be made implicit by using the `estimateAutoARIMA` with the same keywords. For example:

    aM = estimateAutoARIMA(x, d=1, seasonal=true)

is equivalent to the above command.


#### ARIMAModel / Forecast and Simulation

The output of any estimation is of type `ARIMAModel`, so it is necessary to go into a bit more depth on this type. The most important fields are:

* arCoef::Vector{Float64}=Array(Float64, 0) (Vector of AR coefficients)
* maCoef::Vector{Float64}=Array(Float64, 0) (Vector of MA coefficients)
* sarCoef::Vector{Float64}=Array(Float64, 0) (Vector of seasonal AR coefficients)
* smaCoef::Vector{Float64}=Array(Float64, 0) (Vector of seasonal MA coefficients)
* intercept::Float64=0.0 (Model intercept. Equal to NaN if no intercept)
* sigma2::Float64=1.0 (Variance of residuals)
* covCoef::Matrix{Float64}=Array(Float64, 0, 0) (Covariance matrix of coefficients (and intercept if applicable))
* covCoefName::Vector{ASCIIString}=Array(ASCIIString, 0) (Names corresponding to the rows/columns of covCoef)
* d::Int=0 (Order of integration)
* D::Int=0 (Seasonal order of integration)
* seasonalPeriod::Int=1 (seasonal period of the data)
* resid::Vector{Float64}=Array(Float64, 0) (Vector of residuals from estimation step)
* optimCode::Int=0 (Output of optimization code from estimation step. 0 implies everything went well, a positive number implies a problem with convergence, and a negative number implies an unknown problem in the estimation step)

An `ARIMAModel` can be constructed using keywords that correspond to the fields, eg:

    aM = ARIMAModel(arCoef=[0.6, 0.1], smaCoef=[0.1], D=1)

Note that in this case, the fields `covCoef` and `covCoefName` will have appropriate size properties to correspond to the number of coefficients in the model, but that the values will all be automatically initialised to 0.0 or undefined strings respectively. `aM` can then be used to as input to the `simulate` function, eg

    x = simulate(numObs, aM)

will simulate `numObs` steps of a univariate ARIMA time-series corresponding to the model in `aM`. Alternatively, one can construct the `ARIMAModel` implicitly using keyword functionality with the `simulateARIMA` function. For example:

    x = simulateARIMA(numObs, arCoef=[0.6, 0.1], smaCoef=[0.1], D=1)

is equivalent to the above call to `simulate`. Two important points worth emphasizing:

* the simulate functions are implemented entirely in Julia, and
* unfortunately the simulate functions are not currently capable of simulating seasonal ARIMA models.

In many cases, an `ARIMAModel` will be constructed naturally as the return output of an estimation function. In this case it can be used to forecast future values of an ARIMA time series. One can use:

    (pointForecast, bounds) = forecast(aM, numStep)

to forecast up to `numStep` observations ahead using the `ARIMAModel` specified in `aM` (or just `forecast(aM)` for a one-step ahead forecsat). The point forecast is returned in `pointForecast` as a `Vector{Float64}`, while 95% confidence bounds are returned in `bounds` as a `Matrix{Float64}`.

IMPORTANT: notice that the underlying data, `x`, is not passed into this method of the `forecast` function. This is because this method assumes that `aM` was generated using one of the R-based estimation functions, and so the embedded R program created by `RCall` already holds both the model and data necessary to perform the forecast. All `aM` really does here is point R to the appropriate variable that already exists in memory.

For those who want to know what is happening "under the hood" here, read this paragraph: the type `ARIMAModel` contains a field `rModelName`. This field is populated during an estimation routine with a unique string that corresponds to the variable name given to the output of an R estimation function (this output sits in the memory space of the embedded R). Thus when `forecast` is called using the above method, the function essentially provides the contents of `rModelName` to the embedded R and says "go find this variable and use it to forecast `numStep` ahead". If it can't a variable corresponding to `rModelName` then an error is thrown.

However, one does not have to use R functions here. One could also use:

    (pointForecast, bounds) = forecast(x, aM, numStep)

or `forecast(x, aM)` for one-step ahead. This method of `forecast` is implemented entirely in Julia, and uses the fields of `aM` and the underlying data `x` to create the forecasts. Unfortunately, the computations necessary to build confidence bounds for the point forecast have not yet been implemented. Thus this method is currently only useful if you're only interested in `pointForecast`. 

This concludes the basic instructions on usage of `RARIMA`. The next section contains a full description of each of the types `ARIMAInput`, `AutoARIMAInput`, and `ARIMAModel`, and is only necessary reading for those who want all the details or who plan on contributing to this package (or using the structure of this package as basis for their own package).


## MOST USERS DO NOT NEED TO READ BEYOND THIS POINT


## Type Descriptions

#### ARIMAInput

Type
* `ARIMAInput`

Fields:
* p::Int: Number of AR lags
* d::Int: Order of integration
* q::Int: Number of MA lags
* P::Int: Number of seasonal AR lags
* D::Int: Order of seasonal integration
* Q::Int: Number of seasonal MA lags
* seasonalPeriod::Int: Seasonal period of data (only relevant if P, D, or Q are > 0)
* includeIntercept::Bool: Include intercept term in model. Not relevant for d > 0.
* transformParams::Bool: AR parameters transformed to ensure they remain in region of stationarity
* fixedParams::Vector{Float64}: Array of size p+q+P+Q containing fixed values for coefficients. Set value to NaN if you don't want to fix a particular parameter. Set to empty array to not fix any parameters during estimation.
* initParams::Vector{Float64}: Array of size p+q+P+Q of initial values for estimation. Set value to NaN if you don't have a initial value for a particular parameter (0 will be used). Set to empty array to not supply any initial values.
* method<:validMethodType: Variable that indicates method of estimating parameters. If ASCIIString, set to "CSS-ML", "CSS", or "ML". See R docs for arima function for more detail. Other types for this variable not coded up yet.
Purpose: This type stores all of the input, excluding the data, needed to estimate the parameters of an ARIMA model.

Constructors:
* ARIMAInput{T<:validMethodType}(p::Int, d::Int, q::Int, P::Int, D::Int, Q::Int, seasonalPeriod::Int, includeIntercept::Bool, transformParams::Bool, fixedParams::Vector{Float64}, initParams::Vector{Float64}, method::T): Constructor for every field
* ARIMAInput(): White noise constructor

Keyword constructor:
* ARIMAInput(; p::Int=0, d::Int=0, q::Int=0, P::Int=0, D::Int=0, Q::Int=0, seasonalPeriod::Int=1, includeIntercept::Bool=true, transformParams::Bool=true, fixedParams::Vector{Float64}=Array(Float64, 0), initParams::Vector{Float64}=Array(Float64, 0), method::ASCIIString="CSS-ML")

Methods:
* copy(a::ARIMAInput): Standard copy
* deepcopy(a::ARIMAInput): Standard deepcopy
* show(io::IO, a::ARIMAInput): Standard show
* show(a::ARIMAInput): Show defaulting to STDOUT


#### AutoARIMAInput

Type:
* `AutoARIMAInput`

Fields:
* d::Int: order of integration (use to fix a specific value). Set to negative value if d not fixed.
* D::Int: order of seasonal integration (use to fix a specific value). Set to negative value if D not fixed.
* pMax::Int: maximum number of AR lags allowed in search procedure
* qMax::Int: maximum number of MA lags allowed in search procedure
* PMax::Int: maximum number of seasonal AR lags allowed in search procedure
* QMax::Int: maximum number of seasonal MA lags allowed in search procedure
* numCoefMax::Int: maximum number of coefficients total allowed in search procedure (not relevant for step-wise search procedure)
* dMax::Int: maximum order of integration allowed in search procedure (over-ridden by d)
* DMax::Int: maximum seasonal order of integration allowed in search procedure (over-ridden by D)
* pStart::Int: Initial value for p in step-wise search procedure
* qStart::Int: Initial value for q in step-wise search procedure
* PStart::Int: Initial value for P in step-wise search procedure
* QStart::Int: Initial value for Q in step-wise search procedure
* stationary::Bool: Set to true to restrict search procedure to stationary models
* seasonal::Bool: Set to false to restrict search procedure to non-seasonal models
* ic::ASCIIString: Information criteria to use in search procedure. Use "aic", "aicc", or "bic" (see R auto-arima docs for more detail)
* stepwise::Bool: If true, use step-wise search procedure. If false, check every feasible model (slower).
* parallel::Bool: If true and stepwise=false, then do checks in parallel if possible.
* numCore::Int: Number of cores to use if doing checks in parallel.

Purpose
* This type stores all of the input, excluding the data, needed to automatically select an ARIMA model and estimate it

Constructors:
* AutoARIMAInput(d::Int, D::Int, pMax::Int, qMax::Int, PMax::Int, QMax::Int, numCoefMax::Int, dMax::Int, DMax::Int, pStart::Int, qStart::Int, PStart::Int, QStart::Int, stationary::Bool, seasonal::Bool): Default methodology options
* AutoARIMAInput(): No seasonality, all method options set to default, all order and step-wise start values set to default.

Keyword constructor:
* AutoARIMAInput(; d::Int=-1, D::Int=-1, pMax::Int=5, qMax::Int=5, PMax::Int=2, QMax::Int=2, numCoefMax::Int=5, dMax::Int=2, DMax::Int=1, pStart::Int=2, qStart::Int=2, PStart::Int=1, QStart::Int=5, stationary::Bool=false, seasonal::Bool=false, ic::ASCIIString="aic", stepwise::Bool=true, parallel::Bool=false, numCore::Int=2)

Methods:
* copy(a::AutoARIMAInput): Standard copy
* deepcopy(a::AutoARIMAInput): Standard deepcopy
* show(io::IO, a::AutoARIMAInput): Standard show
* show(a::AutoARIMAInput): Show defaulting to STDOUT


#### ARIMAModel

Type:
* ARIMAModel

Fields:
* arCoef::Vector{Float64}: Vector of AR coefficients
* maCoef::Vector{Float64}: Vector of MA coefficients
* sarCoef::Vector{Float64}: Vector of seasonal AR coefficients
* smaCoef::Vector{Float64}: Vector of seasonal MA coefficients
* intercept::Float64: Model intercept. Equal to NaN if no intercept.
* sigma2::Float64: Variance of residuals
* covCoef::Matrix{Float64}: Covariance matrix of estimated coefficients (and intercept if applicable)
* covCoefName::Vector{ASCIIString}: Names corresponding to the rows/columns of covCoef
* d::Int: Order of integration
* D::Int: Seasonal order of integration
* seasonalPeriod::Int: Seasonal period.
* logLik::Float64: Log-likelihood at maximum as reported by R in estimation step.
* aic::Float64: aic at maximum as reported by R in estimation step.
* resid::Vector{Float64}: Vector of residuals from estimation step.
* optimCode::Int: Optimisation output code as reported by R. 0=normal convergence obtained, other values imply something wrong (see R docs for more detail)
* rModelName::ASCIIString: Variable name assigned to model output from R estimation step in embedded R workspace. This field is used to find the R output from estimation step and feed it into R predict.Arima (if the user calls Julia forecast function)

Purpose:
* This type stores a complete ARIMA specification.

Constructors:
* ARIMAModel(): Default values

Keyword constructor:
* ARIMAModel(; arCoef::Vector{Float64}=Array(Float64, 0), maCoef::Vector{Float64}=Array(Float64, 0), sarCoef::Vector{Float64}=Array(Float64, 0), smaCoef::Vector{Float64}=Array(Float64, 0), intercept::Float64=0.0, sigma2::Float64=1.0, covCoef::Matrix{Float64}==Array(Float64, 0, 0), covCoefName::Vector{ASCIIString}=Array(ASCIIString, 0), d::Int=0, D::Int=0, seasonalPeriod::Int=1, logLik::Float64=0.0, aic::Float64=0.0, resid::Vector{Float64}=Array(Float64, 0), optimCode::Int=0, rModelName::ASCIIString="")

Methods:
* copy(a::ARIMAModel): Standard copy
* deepcopy(a::ARIMAModel): Standard deepcopy
* show(io::IO, a::ARIMAModel): Standard show
* show(a::ARIMAModel): Show defaulting to STDOUT
* coef(a::ARIMAModel): Returns all coefficients (and intercept if applicable) in one vector ordered [arCoef, maCoef, sarCoef, smaCoef, intecept (if applicable)]
* vcov(a::ARIMAModel): Return the variance covariance matrix of coefficients
* coefVar(a::ARIMAModel): Return a 2x? matrix containing all coefficients in first row (ie output of coef functions) and the corresponding variances of these estimators in the second row.


## Function descriptions

##### estimation functions

Functions:
* estimate
* estimateARIMA
* estimateAutoARIMA

Methods:
* {T<:Number}(x::Vector{T}, a::ARIMAInput, method::ASCIIString): Estimate an ARIMAModel using data x and estimation inputs in ARIMAInput. Do estimation using R arima function with method determined by method ASCIIString.
* estimate{T<:Number}(x::Vector{T}, a::ARIMAInput, method::Symbol): Estimate an ARIMAModel using data x and estimation inputs in ARIMAInput. Do estimation using Julia Optim package. WARNING: NOT YET IMPLEMENTED.
* estimate{T<:Number}(x::Vector{T}, a::ARIMAInput): Wrapper obviating the need to specify method as separate argument (since method is part of ARIMAInput)
* estimate{T<:Number}(a::ARIMAInput, x::Vector{T}): Wrapper ensuring order of input doesn't matter.
* estimate{T<:Number}(x::Vector{T}, a::AutoARIMAInput): Automatically choose order and lags of ARIMAModel using inputs in a::AutoARIMAInput and data in x. This method uses the R auto.arima function from Rob Hyndman's forecast library.
* estimate{T<:Number}(a::AutoARIMAInput, x::Vector{T}): Wrapper ensuring order of input doesn't matter.

Keyword methods:
* estimateARIMA{T<:Number}(x::Vector{T}; p::Int=0, d::Int=0, q::Int=0, P::Int=0, D::Int=0, Q::Int=0, seasonalPeriod::Int=1, includeIntercept::Bool=true, transformParams::Bool=true, fixedParams::Vector{Float64}=Array(Float64, 0), initParams::Vector{Float64}=Array(Float64, 0), method::ASCIIString="CSS-ML")
* estimateAutoARIMA{T<:Number}(x::Vector{T}; d::Int=-1, D::Int=-1, pMax::Int=5, qMax::Int=5, PMax::Int=2, QMax::Int=2, numCoefMax::Int=5, dMax::Int=2, DMax::Int=1, pStart::Int=2, qStart::Int=2, PStart::Int=1, QStart::Int=5, stationary::Bool=false, seasonal::Bool=false, ic::ASCIIString="aic", stepwise::Bool=true, parallel::Bool=false, numCore::Int=2)

Output:
* ARIMAModel containing the estimated ARIMA coefficients and other variables of interest.

Purpose:
* The purpose of this function is to estimate the parameters of the input model and return a complete specification of the model



#### forecast methods

Functions:
* forecast

Methods:
* forecast(a::ARIMAModel, numStep::Int=1): Forecast numStep ahead using ARIMAModel. Note, this method does not input the data because it relies on finding the correct ARIMA model in embedded R workspace using a.rModelName and then calling predict.Arima on this model. This method is needed because the method written entirely in Julia cannot handle bounds computations.
* forecast{T<:Number}(x::Vector{T}, a::ARIMAModel, numStep::Int=1): Forecast numStep ahead using ARIMAModel and data in x. This method is implemented entirely in Julia, but is not able to handle seasonal ARIMA models or bounds computations.
* forecast{T<:Number}(a::ARIMAModel, x::Vector{T}, numStep::Int=1): Wrapper so that order of inputs doesn't matter.

Output:
* Output is tuple (pointFore::Vector{Float64}, boundFore::Matrix{Float64}), where pointFore provides the point forecast up to numStep ahead, and boundFore provides 95% lower and upper bounds for the forecasts (lower bound in first column, upper bound in second column)

Purpose:
* The purpose of this function is to provide forecasts for the input dataset and ARIMA model


#### simulation functions

Functions:
* simulate
* simulateARIMA

Methods:
* simulate(numObs::Int, a::ARIMAModel, numBurn::Int=-1, resid::Vector{Float64}=[NaN]): Simulate numObs observations using ARIMAModel a, with burn-in length numBurn and optionally use resid as the source of random numbers in simulation.
* simulate(a::ARIMAModel, numObs::Int, numBurn::Int=-1, resid::Vector{Float64}=[NaN]): Wrapper so order doesn't matter.

Keyword method:
* simulateARIMA(numObs::Int; arCoef::Vector{Float64}=Array(Float64, 0), maCoef::Vector{Float64}=Array(Float64, 0), d::Int=0, intercept::Float64=NaN, sigma2::Float64=1.0, numBurn::Int=-1, resid::Vector{Float64}=[NaN]): Keyword activated simulate method. First input specifiying number of observations to simulate is non-optional. Default keywords result in white noise model.

Output:
* Output is Vector{Float64} of data simulated using the specified ARIMA model.

Purpose:
* The purpose of this function is to simulate the specified ARIMA model

Notes:
* This function is not currently capable of simulating seasonal ARIMA models, ie non-seasonal only.
