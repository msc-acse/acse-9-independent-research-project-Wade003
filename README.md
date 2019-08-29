## Introduction
Reduced order modelling is a powerful technique for rapidly modelling high dimensional fluid dynamics systems.
Its speed could enable real-time decision making and operational modelling.

In this repository you can find the tools for reduced order modelling.

Example problems are presented and solved, including:
 - A simple channel flow navier-stokes equation (see: demo/flow_past_cylinder_demo.ipynb)
 - A poisson equation (see: demo/poisson_temp_demo.ipynb)
 - A three component chemical reactions  equation coupled to the navier-stokes equation (see: demo/cylinder_rxn_demo.ipynb)
 - A radionuclide transport equation coupled to the navier-stokes equation (see: demo/radio_transport_demo.ipynb)
 - A radionuclide transport equation coupled to the shallow water equation (see: demo/hydrodynamics_demo.ipynb)
 - A method of manufactured solutions verification on the radionuclide transport problem (see: demo/MMS_transport/demo.ipynb)
 
## System requirement
Linux.

## Installation instructions
To download the repository to your local machine.
```bash
  git clone https://github.com/msc-acse/acse-9-independent-research-project-Wade003.git
```
Installing all packages below.
```bash
   pip install -U numpy
   pip install sklearn 
   pip install keras
   pip install tensorflow
   pip install matplotlib
   pip install torch==1.2.0+cpu torchvision==0.4.0+cpu -f https://download.pytorch.org/whl/torch_stable.html
   pip install vtk
```
Adding environment path for Opal and IC-Ferst
```bash
   export PYTHONPATH='/data/wade/test/multifluids_icferst-master/python:$PYTHONPATH'
   export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/data/wade/test/Opal-master/spud
   export PATH="/data/wade/test/Opal-master/spud:$PATH"

```


## Dependencies
 The external libraries:

 - numpy >= 1.16.4
 - scipy >= 1.2.2
 - matplotlib >= 3.1.0
 - jupyter >= 1.0.0 (to run demo notebooks)
 - Sphinx >= 1.8.5 (to compile documentation)
 - Keras >= 2.2.4
 - scikit-learn >= 0.20.3
 - tensorflow >= 1.14.0   
 - torch >= 1.1.0    
 - torchvision >= 0.3.0
 - vtk >= 8.1.2
## Repository Information
* __animation__		- recorded animation video for several test case
* __images__		- report related images
* __software__		- contains Opal (which is the major software I developed) and IC-Ferst (support to run the Opal)
* __report__		- contains the final report, detailing this project's motivations, software design, analysis, and conclusions 
* __test_cases__		- 1D square wave and 2D flow past a cylinder test cases

## Author and Course Information
__Author:__ Wade Song
__Github:__ Wade003
__CID:__ 01569843

This project is completed for Imperial College's MSc in Applied Computational Science and Engineering program,
as part of the final course module ACSE9. This project was completed under the supervision of Professor Christopher Pain. 
## License  
Licensed under the MIT [license](https://github.com/msc-acse/acse-9-independent-research-project-Wade003/blob/master/LICENSE)