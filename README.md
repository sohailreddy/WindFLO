# WindFLO - A Framework for Fast Wind Farm Layout Optimization
Copyright, 2019 Sohail R. Reddy  (sredd001@fiu.edu)

## Version v2.1.0 with source code

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

The repository contains the source code and the python API for WindFLO.


### Prerequisites

The framework requires libgfortran and f90nml (for Python API). It is compiled and linked using the g++ compiler.


### Installing

* Enter the 'src/' directory and update the c++ and Fortran compilers in the Makefile, then run make with your OS. This will create the executable, the static and dynamic libraries in the 'release/' and 'API/' directories for your OS. 

```
cd src; make OS=OSX
```

### Running WindFLO

WindFLO requires an input file containing all the parameters needed to run the model. The parameters that need to be specified are given in the sample input file (input/WindFLO.dat) and the manual. The input directory contains the input files needed to run WindFLO. The executable can be ran using the command

```
./WindFLO /path/to/the/input/file /path/where/you/want/the/result/file
```

## How to Cite

If you use WindFLO, please cite the following papers
```
S.R. Reddy, "Wind Farm Layout Optimization (WindFLO) : An Advanced Framework for Fast Wind Farm Analysis and Optimization," Applied Energy, Vol. 269, pp 115090, 2020
```
```
S.R. Reddy, “An Efficient Method for Modeling Terrain and Complex Terrain Boundaries in Constrained Wind Farm Layout Optimization,” Renewable Energy, vol. 165, pp. 162–173, Mar. 2021
```
```
S.R. Reddy, “A Machine Learning Approach for Modeling Irregular Regions with Multiple Owners in Wind Farm Layout Design,” Energy, vol. 220, Apr. 2021
```
## Contributing

The WindFLO is an open source framework and was developed to accelerate the dissemination of results and development of models for wind farm analysis. The framework allows users to develop and incorporate their own wake models into the WindFLO framework. The user defined wake models and wake merge schemes can be added to the *userDefined.cpp*. A user-defined cost model can also be included in the same source file. The *userDefined.cpp* in this repository provides an example on how to incorporate a new wake model.

To include your wake model in the repo, please fork the repo and create a pull request.


## Authors

**Sohail R. Reddy**


## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details
