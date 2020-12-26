[![Maintenance](https://img.shields.io/badge/Maintained%3F-no-red.svg)](https://bitbucket.org/lbesson/ansi-colors)

## Barotropic QG code

A pseudospectral code for solving the barotropic vorticity equation on a beta-plane.

## Time-stepping Scheme

For time-stepping the exponential time-differencing Runge-Kutta 4th order scheme is used. More details in [these notes](https://github.com/navidcy/ETDRK4_notes/blob/master/ETDRK4-timestepping.pdf)

## Author

Navid C. Constantinou. Scripps Institution of Oceanography, UC San Diego. 2016

## Disclaimer

This repository is not being maintained anymore. 

An alternative is the [Julia](https://www.julialang.org) package [GeophysicalFlows.jl](http://github.com/FourierFlows/GeophysicalFlows.jl) and in particular the module [SingleLayerQG](https://fourierflows.github.io/GeophysicalFlowsDocumentation/stable/modules/singlelayerqg/). Check out [this example](https://fourierflows.github.io/GeophysicalFlowsDocumentation/stable/generated/singlelayerqg_betaforced/).
