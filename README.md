# AMM: Analytic Modal Method for Diffraction Grating Efficiency Analysis Scripts

##### Programming language: MATLAB

[![Version](https://img.shields.io/badge/version-1.0-green.svg)](README.md) [![License](https://img.shields.io/badge/license-MIT-blue.svg)](http://opensource.org/licenses/MIT)

---
## Program goals
The aim of the AMM scripts is to provide examples of the AMM method for diffraction efficiency calculations.  Complete examples are given for Duty Cycle sweep, Theta (angle of incidence) sweep, and Wavelength (lambda) sweep. Both TE and TM modes are addressed. The scripts execute serially and can take significant time to execute. The exact analytic modes are computed (up to some specified number). The algorithms involved lend themselves to parallelization, offering significant speed-up opportunity. The goal is to motivate development of AMM for parallel (cluster or gpu) processing.

## Solution method
AMM is implemented following the lead of I. Botten, M Craiag, R. McPhedran, J. Adams, J. Andrewartha, ‘The dielectric lamellar diffraction grating,’ Optica Acta, 28(3) 413-428 (1981). Piecewise-analytic solution of Helmholtz equations are found for gratings modeled in the piecwise constant multi-layer model. 

## Scientific work
1. D. Fluckiger, 'The Analyltic Modal Method Algorithm,'

---
## Manual
1. Three 'stand-alone' scripts implement AMM, one for each of a Duty Cycle, Theta, Lambda parameter sweep.
    * AMM_mainDCSweep
    * AMM_mainThetaSweep
    * AMM_mainLambdaSweep
2. Free parameters (tuning) are found in the set-up of the rectangular root search region in FindRootsGRPF.m
3. Free parameters (tuning) are found in root-polishing rootYasmin.m which control root resolution, number of iterations, root location failure limits. These limits are currently hard-wired, but experimentation is encouraged to find optimal values, which likely depend on grating layer parameters.
 
## Short description of the scripts and files
-[DC sweepCWA_AMM_51modes_A.mp4]DC sweepCWA_AMM_51modes_X.mp4 - Animation of couple wave (CWA) DC sweep with N=51 orders, CWA roots show as red stars, AMM true roots as green circles. There are four segments (A, B, C, D) to cover the full sweep

-[DC sweepCWA_AMM_201modes_A.mp4]DC sweepCWA_AMM_201modes_X.mp4 - Animation of couple wave (CWA) DC sweep compared with AMM with N = 2011, shows both anomalies and omissions of roots (X = A,B,C,D; break up sweep to 4 segments)

-[AMM_mainDCSweep.m]AMM_mainDCSweep.m - AMM DC sweep of single layer, two region grating

-[AMM_mainLambdaSweep.m]AMM_mainLambdaSweep.m - AMM wavelength sweep, multilayer-blaze profile in aluminum (high MP required for first few roots)

-[AMM_mainThetaSweep.m]AMM_mainThetaSweep.m - AMM theta sweep (angle of incidence), multilayer-blaze profile in glass

-[BoundaryMatrix.m]BoundaryMatrix.m - Construct the boundary condition matrix for all layers, as each layer is solved for

-[dEigEq.m]dEigEq.m - Solve for the derivative (with respect to mu) of the eigenvalue equation

-[Drude.m]Drude.m -	Drude model of the index of refraction for several metals

-[EigEq.m]EigEq.m -	Solve the eigenvalue equation at some mu give the layer grating parameters

-[FindNextNode.m]FindNextNode.m - Part of the GRPF root finding algorithm,  finds the next node in the candidate boundary creation process

-[FindRootsGRPF.m]FindRootsGRPF.m - Main driver to find N roots of a layer grating structure

-[GRPF.m]GRPF.m - Main entry point of the GRPF algorithm

-[plotEigEq.m]plotEigEq.m - Generates a contour plot of the log magnitude of the eigenvalue equation

-[QY2x2.m]QY2x2.m -	Solve for the theta and psi, and Field eigenfunctions for a given grating layer structure

-[rectdom.m]rectdom.m - Generate the starting nodes of a rectangular region of the complex plane (for GRPF)

-[rootYasmin.m]rootYasmin.m -	Solve for a root of the eigeq given a starting point ‘guess’ using 4th order method

-[Schott.m]Schott.m - Give index of refraction of a catalogue material according to the Schott formula

-[Sellmeier.m]Sellmeier.m - Ggive index of refraction of a catalogue material according to the Sellmeier formula

-[TableNK.m]TableNK.m -	Give index of refraction for materials in table form using cubic interpolation between entries

-[vinq.m]vinq.m - Rreturns phase quadrant of complex number (used in GRPF)

-[zShift.m]zShift.m - Constrains a complex number to remain within a defined rectangle

## Additional comments
The code involves MATLAB function [delaunayTriangulation](https://uk.mathworks.com/help/matlab/ref/delaunaytriangulation.html) which was introduced in R2013a version. The muti-precision toolbox used is from Advanpix, https://www.advanpix.com  

## Author
The project was developed by David Fluckiger, author and sole proprietor of GSolver, a grating analysis code first introduced for commercial consumption in 1994 based on Coupled Wave Analysis algorithms. https://www.gsolver.com

## License
These scripts are provided under open-source Matlab code licensed using the [MIT license](LICENSE.md).
