## RARIMA

An [ARIMA](http://en.wikipedia.org/wiki/Autoregressive_integrated_moving_average) (Autoregressive Integrated Moving Average) module for the [Julia](http://julialang.org/) language that uses libraries from the [R language](http://www.r-project.org/)


## Main Features

This module allows Julia users to estimate, forecast, and simulate ARIMA time-series models.

Some functions are implemented in pure Julia, while others call functions from the R language. The relevant R functions are accessed using the Julia package [*RCall*](https://github.com/JuliaStats/RCall.jl). These R functions come from either the R [stats](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/00Index.html) library or the R [forecast](http://cran.r-project.org/web/packages/forecast/index.html) library. Since these R functions are available via the GNU General Public License (v2 or v3), this module uses the GNU GPL v2.


## Setting up *RCall*

It is recommended that you set up *RCall* and verify it works prior to installing *RARIMA*. *RCall* is an officially registered package, so it can be installed from the Julia REPL using `Pkg.add("RCall")`. For more detail, please refer to the [*RCall* docs](https://github.com/JuliaStats/RCall.jl).

You can verify *RCall* is working correctly using the following commands at the Julia REPL:

    using RCall
    rcopy("1+2")

which will hopefully return the obvious answer. Also note, you *must* install the *forecast* library into your R library before using *RARIMA*. This can be done from an R terminal (or R IDE of choice, e.g. R-Studio) using the command:

    install.packages("forecast")

to install the library, then:

    library(forecast)

to ensure it loads.


## A Word on Julia versus R

The long-term goal is to implement ARIMA time-series modelling capability purely in Julia (preferably under an MIT license rather than GPL). This module exists only to fill the gap between the present day and the day that future goal is realised. I've tried to structure *RARIMA* so that it can serve as the basis of that future module, ie it utilises Julia's type system and multiple dispatch, and have implemented some of the functionality already in Julia.

The most significant functionality that needs to be converted to Julia is the estimation routines (regular and automatic order selection). Currently, these are implemented entirely in R. Because the estimation step typically involves solving a numerical optimisation problem, sometimes things go wrong. One of the main challenges in writing this module is dealing with all the different types of errors that can result from the estimation step in R and passing something informative back to Julia. I haven't spent a huge amount of time on this since the long-term goal is to replace that code with Julia code. Currently, the type returned from estimation functions contains a field `optimCode` that will be set to `0` if everything worked correctly. If there was a known problem from the R end, this field will take on a positive value, and the corresponding optimisation code can be looked up in the R documentation. If something went wrong that the Julia code is not able to understand, this field will be set to `-1`. Watch out for this!


## Quick Start

To install *RARIMA*, first ensure *RCall* is installed and working, and that your R library contains the library *forecast*. See the above section "Setting up *RCall*" for more detail.

Since *RARIMA* is not currently an official Julia package, you can install it using:

    Pkg.clone("https://github.com/colintbowers/RARIMA.jl")

at the Julia REPL. If this is successful, you can now load the package using:

    using(RARIMA)

The nuts-and-bolts method for working with this module is via three types, `ARIMAInput`, `AutoARIMAInput`, and `ARIMAModel`, and three generic functions `estimate`, `forecast`, and `simulate`. However, keyword-based functionality is available, and this is probably how most users will choose to interact with this package. The keywords are based on the fields of the above three types. A short-list of the most important fields and the associated functionality is provided now.

#### ARIMAInput / Estimation

The most important field in the type `ARIMAInput` are:

* p::Int=0 (autoregressive order) 
* d::Int=0 (order of integration)
* q::Int=0 (moving average order)
* P::Int=0 (seasonal autoregressive order)
* D::Int=0 (order of seasonal integration)
* Q::Int=0 (seasonal moving average order)
* seasonalPeriod::Int=1 (seasonal period - only relevant if P, D, or Q is > 0)
* includeIntercept::Bool=true (include an intercept term in the model - only relevant if d and D equal 0)

So one could construct an `ARIMAInput` for an ARIMA(2,1,0,1,0,0) model of monthly data with no intercept using:

    aI = ARIMAInput(p=2, d=1, P=1, seasonalPeriod=12, includeIntercept=false)

Actually, in this case we could have omitted `includeIntercept` because it is automatically set to false when `d > 0` (see the R docs for `arima` for more detail).

Given a vector of observed data, `x::Vector{T}` where `T<:Number`, we can estimate the parameters of a model specified by and `ARIMAInput` using:

    aM = estimate(x, aI)

where the output is of type `ARIMAModel` (discussed later in this section). Note, `estimate(aI, x)` is equivalent. Alternatively, one can make the construction of the `ARIMAInput` implicit by calling the function `estimateARIMA` using keywords. For example, an equivalent command to the one above would be:

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

Note that in this case, the fields `covCoef` and `covCoefName` will have appropriate size properties to correspond to the number of coefficients in the model, but that the values will all be automatically initialised to 0.0 or an undefined string respectively. It can then be used to as input to the `simulate` function, eg

    x = simulate(numObs, aM)

will simulate `numObs` steps of a univariate ARIMA time-series corresponding to the model in `aM`. Alternatively, one can construct the `ARIMAModel` implicitly using keyword functionality with the `simulateARIMA` function. For example:

    x = simulateARIMA(numObs, arCoef=[0.6, 0.1], smaCoef=[0.1], D=1)

is equivalent to the above call to `simulate`.

In many cases, an `ARIMAModel` will be constructed naturally as the return output of an estimation function. In this case it can be used to forecast future values of an ARIMA time series.

Usually an `ARIMAModel` will be constructed naturally as the output of an estimation function. 





#### Estimation

The first non-optional input to `estimateARIMA` must be your data expressed as type `Vector{T}`, where `T<:Number`. The remaining inputs are keyword arguments, some of which follow (all possible keyword arguments can be found by examining the fields of the type `ARIMAInput`): 

For example, given vector of seasonal data `x` with seasonal period of 4, an ARIMA(2,1,1,0,0,1) could be estimated using:

`M = estimateARIMA(x, p=2, d=1, q=1, Q=1, seasonalPeriod=4)`

Alternatively, if the user is not sure of the order and lags to use, these can be automatically detected using , which also requires the data as first input, and then exhibits the following keyword arguments (again this is a short list - for all keyword arguments see the fields of the type `AutoARIMAInput`):


For example, given vector of data `x` that is suspected to exhibit a seasonal component use:

`M = estimateAutoARIMA(x, seasonal=true)`

For the above two functions, the output `M` is of type `ARIMAModel`. This type is a complete specification of an ARIMA model, including coefficients values etc. A short-list of the fields of this type includes (for a full list of fields see the type `ARIMAModel`):


#### Forecast

One can generate the 5-step ahead forecast from an ARIMAModel using

`(pointForecast, forecastBound) = forecast(M, 5)`

or simply `forecast(M)` for a one-step ahead forecast. The point forecasts are returned in `pointForecast` as type `Vector{Float64}`. 95% confidence bounds for the point forecasts are provided in the matrix `forecastBound` which is of type `Matrix{Float64}`. The first column contains lower bounds and the second column contains upper bounds.

The `forecast` function is also implemented in pure Julia and can be called by including the data in the function call, e.g. `forecast(x, M, 5)`. However, this method is only implemented for non-seasonal models, and it also is not yet able to compute confidence bounds.

#### Simulation

The `simulateARIMA` function simulates the specified ARIMA model. The first (non-optional) argument is the number of steps to simulate, expressed as type `Int`. The remaining arguments are keyword arguments. A full list follows:

* arCoef::Vector{Float64}=Array(Float64, 0) (Vector of autoregressive coefficients of model to simulate)
* maCoef::Vector{Float64}=Array(Float64, 0) (Vector of moving average coefficients of model to simulate)
* d::Int=0 (Order of integration of model to simulate)
* intercept::Float64=NaN (Intercept term of model to simulate. Set to `NaN` for no intercept)
* sigma2::Float64=1.0 (Variance of residuals in model to simulate)
* numBurn::Int=-1 (Number of observations in the "burn-in" period - these observations are used to warm-up the procedure and are discarded at the conclusion of the function). Default is to currently use p+q+100. This is rudimentary and hopefully a model-driven default choice will be implemented soon.
* resid::Vector{Float64}=\[NaN\] (Vector of numbers to use as the source of randomness in the model. Default value of `[NaN]` implies that `rnorm` will be used to generate random numbers)

This function returns the simulated data as type `Vector{Float64}`. For example, once could simulate 1000 observations of an ARIMA(2,1,1) model with no intercept using:

`xSimulated = simulateARIMA(1000, arCoef=[0.5, 0.2], maCoef=[0.1], d=1)`

The simulate functions in RARIMA are implemented entirely in Julia. Unfortunately the ability to simulate seasonal ARIMA models has not yet been added.

 This concludes the quick start. As discussed, keyword arguments are not the preferred way of working with the functions in this module. Ideally the user will engage with the three types `ARIMAInput`, `AutoARIMAInput`, and `ARIMAModel` and how they can be used with the generic functions `estimate`, `forecast`, and `simulate`.
 
## Module Types
 
 More detailed info to come. In the meantime feel free to examine the source code in /src/RARIMA.jl. Each type and function is extensively documented within the source code.
 

[![Build Status](https://travis-ci.org/colintbowers/RARIMA.jl.svg?branch=master)](https://travis-ci.org/colintbowers/RARIMA.jl)
