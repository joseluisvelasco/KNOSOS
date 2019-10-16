# KNOSOS
================

The KiNetic Orbit-averaging SOlver for Stellarators (KNOSOS)[J. L. Velasco, I. Calvo, F. I. Parra and J. M. García-Regaña, KNOSOS: a fast orbit-averaging neoclassical code for arbitrary stellarator geometry, submitted to J. Comp. Phys, arXiv:1908.11615 [physics.plasm-ph]](https://arxiv.org/abs/1908.11615) is a freely available, open-source code that calculates neoclassical transport in low-collisionality plasmas of three-dimensional magnetic confinement devices by solving the radially local drift-kinetic and quasineutrality equations. The main feature of KNOSOS is that it relies on orbit-averaging to solve the drift-kinetic equation very fast. KNOSOS treats rigorously the effect of the component of the magnetic drift that is tangent to magnetic surfaces, and of the component of the electrostatic potential that varies on the flux surface, φ1. Furthermore, the equation solved is linear in φ1, which permits an efficient solution of the quasineutrality equation. As long as the radially local approach is valid, KNOSOS can be applied to the calculation of neoclassical transport in stellarators (helias, heliotrons, heliacs, etc.) and tokamaks with broken axisymmetry.

------------------------------------------------------------------------

### Theory

The theory behind the code is described in section 2 of the paper, as well as in chapter 3 of the user manual. More details can be found in the paper [I. Calvo, F. I. Parra, J. L. Velasco and J. A. Alonso, The effect of tangential drifts on neoclassical transport in stellarators close to omnigeneity, Plasma Phys. Control. Fusion 59, 055014 (2017)](https://doi.org/10.1088/1361-6587/aa63ce).

------------------------------------------------------------------------

### Installation and compilation

Chapter 5 of the user manual explains how to install and compile the code.

------------------------------------------------------------------------

### Input Data Format

Chapter 7 of the user manual explains how the input data is organized. Chapters 8 and 9 provide more details.

------------------------------------------------------------------------

### Execution

Chapter 6 of the user manual explains how to execute the code and what it produces.

------------------------------------------------------------------------

### Output Data Format

Chapter 10 of the user manual explains how the output data is formatted.

------------------------------------------------------------------------

### Visualization

The data can be visualized with the gnuplot scripts available in folder KNOSOS/PLT/.

------------------------------------------------------------------------

### Tutorials

The user manual can be found in folder KNOSOS/MANUAL/.
