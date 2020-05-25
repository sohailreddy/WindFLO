# WindFLO - A Framework for Fast Wind Farm Layout Optimization
Copyright, 2019 Sohail R. Reddy  (sredd001@fiu.edu)

## Version v2.0.0

<img src="/images/Turbine.png" width="210">             <img src="/images/WindFarmLayout.png" width="400">

The WindFLO model incorporates the following advancements:

* Local terrain elevation is modeled.
* The atmospheric wind velocity is modeled using the Log-Law, the Power-Law and Deaves-Harris model.
* Six different analytical wake models are available
* Four different wake merging schemes are available
* The properties of each turbine (diameter, height, power curve, etc.) can be different. 
* Each turbine can use a different wake model and wake superposition scheme.
* Partial influence of the wake on a turbine is also considered.
* Rotor diameter and tower height dependent cost models are available.
* Wind farm area is computed as the area of convex hull.

<img src="/images/TurbineOverLand.png" width="400"> <img src="/images/WakeCone.png" width="450">

## Getting Started

The repository contains two options for running WindFLO. 

* The executubles can be downloaded for the respecitive operating system (OSX and Linux).
* The precompiled library can be download


### Prerequisites

The framework only requires libgfortran. It is compiled and linked using the g++ compiler.


### Installing

* If the executables are used, there is no installation procedure. These executable only contains the models described in the paper

* If the precompiled library is being used, the user-defined source files also need to be compiled and linked. Simply running the make command in the downloaded directories will compile, link and create the WindFLO executable. No further installation is needed. 

```
cd release; make OS=OSX MAIN=main
```
```
cd release; make OS=LINUX MAIN=main
```

### Running WindFLO

WindFLO requires an input file containing all the parameters needed to run the model. The parameters that need to be specified are given in the sample input file (input/WindFLO.dat) and the manual. The input directory contains the input files needed to run WindFLO. The executable can be ran using the command

```
./WindFLO_OSX /path/to/the/input/file /path/where/you/want/the/result/file
```

## How to Cite

If you use WindFLO, please cite the following paper
```
S.R. Reddy, "Wind Farm Layout Optimization (WindFLO) : An Advanced Framework for Fast Wind Farm Analysis and Optimization," Applied Energy, Vol. 269, pp 115090, 2020
```

## Contributing

The WindFLO framework was developed to accelerate the dissemination of results and development of models for wind farm analysis. Using the precompiled library allows users to develop and incorporate their own wake models into the WindFLO framework. The user defined wake models and wake merge schemes can be added to the *userDefined.cpp* and linked with the WindFLO library. A user-defined cost model can also be included in the same source file. The *userDefined.cpp* in this repository provides an example on how to incorporate a new wake model.

To include your wake model in the next release, please contact me at sredd001@fiu.edu 


## Authors

**Sohail R. Reddy**


## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details
