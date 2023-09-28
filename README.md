# Model Reduction in Soft Robotics Using Locally Volume-Preserving Primitives

## Introduction
This is the MATLAB implementation for the reduced-order modeling of soft robots using closed-form locally volume-preserving deformation primitives. We focus on developing and applying primitive deformations that each observe the locally volume preserving constraint. By composing a wide enough variety of such deformations, many of the most common behaviors observed in soft robots can be replicated. 

- Paper: [IEEE Robotics and Automation Letters (RA-L)](https://ieeexplore.ieee.org/abstract/document/10197570)

## Scripts
The codes calculate the mode parameters that satisfy the given boundary conditions for each case.

- [`Chamber.m`](Chamber.m): For Chambers, as __Fig. 4(1)__ in the paper.
- [`PlanarRod.m`](PlanarRod.m): For rods bending in 2D on a horizontal plane (no gravity), as __Fig. 4(2)__ in the paper.
- [`ThreeDRod.m`](ThreeDRod.m): For rods bending in 3D space (neglect gravity), as __Fig. 4(3)__ in the paper.
- [`RodGravity.m`](RodGravity.m): For rods bending in 2D on a vertical plane (with gravity), as __Fig. 4(4)__ in the paper. 
- [`BlockandFilm.m`](BlockandFilm.m): For blocks and thin films, as __Fig. 4(5)__ and __Fig. 4(6)__ in the paper. 
- [`StrainEnergy_Analytical.m`](StrainEnergy_Analytical.m): This computes the analytical strain energy density function.
- [`GetWeightMatrix.m`](GetWeightMatrix.m): This computes the weight matrix in pseudoinverse, as described in paper section 2D. 
