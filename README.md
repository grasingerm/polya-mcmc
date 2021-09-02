# polya-mcmc
Symmetry, finite groups, and Markov chain Monte Carlo methods.

## Introduction
Conceptually, the basic idea here is to leverage discrete symmetries, 
or Hamiltonians which have a kind of weak discrete symmetry, in order to perform 
efficient Markov cain Monte Carlo sampling of energy landscapes with multiple,
separated energy wells--a kind of hybrid computational statistical physics and 
group theory study, if you will.
Instead of writing a single modular, user- and developer-friendly code, 
(which would require time, patience, skill, and planning--none of which I have)
I've opted to put together many quick and simple examples.
Practically, the main goal here is to further research, not develop software. 
There are plenty of Monte Carlo codes out there. I have no intention of 
competing with any of them and couldn't even if I tried. Instead, my hope is to 
develop and illustrate principles which can be incorporated into better, 
more practical, and more mature MCMC codes.

## Getting started
This research code is written in the [julia language](https://julialang.org).
Once julia is installed, you can install the modules that this code depends on by
running

    julia install_dependencies.jl

in the root directory of this project.

## Examples
### refl-ex1.jl
This example is a 1 DOF system with a reflection symmetry.

### trans-ex2.jl
This example is a 1 DOF system with discrete translational symmetries.

### D2-ex3.jl
This example is a 2 DOF system with D<sub>2</sub> symmetry.

### rosenbrock-ex4.jl
This example is a toy 2 DOF system which has a notoriously difficult energy landscape.
It has a "weak" x &#8594; -x symmetry, which is strictly symmetric when a &#8594; 0.

### harddisks-refl-ex5.jl
This is the first example involving interacting particles.
Here we have multiple particles in the same energy landscape with reflection symmetry as seen in refl-ex1.jl.
The particles each have a single degree of freedom and are hard disks (interact via an excluded volume potential).
We use the same group as in refl-ex1.jl.
This example shows the ways in which this method can be generalized to more complex systems.
It is also the first example which resides in a high dimensional space (scales linearly with the number of disks), which motivates the necessity for stochastic integration.
That is, all of the previous examples can be readily solved using quadrature, but quadrature would struggle for this system when there is a large number of particles.

### D2h-ex6.jl
This example represents another push toward more realistic physical systems.
Here we imagine an ion in a crystalline lattice.
It interacts with atoms in the lattice (which are idealized as fixed) via an electrostatic Coloumb potential.
The ion can move in 3 dimensional Euclidean space and the electrostatic energy landscape is described by a D<sub>2h</sub> symmetry (see Tinkham's text on group theory).

## TODO
- [ ] Plot convergence rates using something more publication friendly--like gnuplot
- [ ] Derive and test exact/approximate solutions; test low-dimensional cases against quadrature
- [ ] Implement noninteracting polymer chain with orientational energy
- [ ] Implement interacting polymer chain with orientational energy; use clustering?
- [ ] Question: will clustering type algorithms work for fluid--solid transitions? There must be something in the literature
- [ ] Write descriptions for each of the examples; add page numbers when (if?) manuscript is ever published/posted
