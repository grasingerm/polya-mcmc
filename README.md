# polya-mcmc
Symmetry, discrete groups, and Markov chain Monte Carlo methods.
See https://arxiv.org/abs/2205.00028 for the manuscript.

## Introduction
Conceptually, the idea here is to leverage discrete symmetries, 
or Hamiltonians which have a kind of semi-discrete symmetry, in order to perform 
efficient Markov chain Monte Carlo sampling of energy landscapes with multiple,
separated energy wells--a kind of hybrid computational statistical physics and 
group theory study.
Instead of writing a single modular, user- and developer-friendly code,
I've opted to put together many quick and simple examples.
Practically, the main goal here is to further research, not develop software. 
There are plenty of Monte Carlo codes out there. The primary goal of this code is to 
develop and illustrate principles which can be incorporated into better, 
more practical, and more mature MCMC codes.

## Getting started
This research code is written in the [julia language](https://julialang.org).
Once julia is installed, you can install the modules that this code depends on by
running

    julia install_dependencies.jl

in the root directory of this project.

## Examples
### Simulation parameters
Each example has a command-line interface.
To obtain a list of program and simulation parameters, use

    julia <example-file-name.jl> --help

### ex1-refl.jl
This example is a 1 DOF system with a reflection semi-symmetry (when _f_ small enough).
It is a double well potential given by _U = a x<sup>4</sup> - b x<sup>2</sup> - f x_.
This example is discussed in detail in section III of the manuscript.

### ex2-trans.jl
This example is a 1 DOF system with discrete translational semi-symmetry (_f_ small enough).
It is a periodic potential consisting of many wells given by _U = a_ cos(_n x_) _- f x_.
This example is discussed in detail in section IV.A of the manuscript.

### ex3-D2.jl
This example is a 2 DOF system with D<sub>2</sub> semi-symmetry.

### ex4-rosenbrock.jl
This example is a toy 2 DOF system which has a notoriously difficult energy landscape.
It has a semi-symmetry x &#8594; -x, which is strictly symmetric when a &#8594; 0.

### ex5-harddisks-refl.jl
This is the first example involving interacting particles.
Here we have multiple particles in the same energy landscape with reflection symmetry as seen in refl-ex1.jl.
The particles each have a single degree of freedom and are hard disks (interact via an excluded volume potential).
We use the same group as in refl-ex1.jl.
This example shows the ways in which this method can be generalized to more complex systems.
It is also the first example which resides in a high dimensional space (scales linearly with the number of disks), which motivates the necessity for stochastic integration.
That is, all of the previous examples can be readily solved using quadrature, but quadrature would struggle for this system when there is a large number of particles.

### ex6-D2h.jl
This example represents another push toward more realistic physical systems.
Here we imagine an ion in a crystalline lattice.
It interacts with atoms in the lattice (which are idealized as fixed) via an electrostatic Coloumb potential.
The ion can move in 3 dimensional Euclidean space and the electrostatic energy landscape is described by a D<sub>2h</sub> symmetry (see Tinkham's text on group theory).
This example is discussed in detail in section IV.B of the manuscript.

### ex7-D2h-LJ.jl
Here we imagine a particle in a crystalline lattice.
It interacts with atoms in the lattice (which are idealized as fixed) via a Lennard-Jones potential.
The particle can move in 3 dimensional Euclidean space and the LJ energy landscape is described by a D<sub>2h</sub> symmetry (see Tinkham's text on group theory).

## TODO
- [x] Plot convergence rates using something more publication friendly--like gnuplot
- [x] Derive and test exact/approximate solutions; test low-dimensional cases against quadrature
- [x] Implement noninteracting polymer chain with orientational energy
- [x] Implement interacting polymer chain with orientational energy; use clustering?
- [x] Write descriptions for each of the examples; add page numbers when (if?) manuscript is ever published/posted
- [ ] Create an example where the group action involves rescaling, e.g. U<sub>G</sub> = a sin(b / x)
