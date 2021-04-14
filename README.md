# polya-mcmc
Symmetry, finite groups, and Markov chain Monte Carlo methods.

# Introduction
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

# Examples
## refl-ex1.jl
This example is a 1 DOF system with a reflection symmetry.

## trans-ex2.jl
This example is a 1 DOF system with discrete translational symmetries.

## D2-ex3.jl
This example is a 2 DOF system with D<sub>2</sub> symmetry.

# TODO
- [ ] Plot convergence rates using something more publication friendly--like gnuplot
- [ ] Derive and test exact/approximate solutions; test low-dimensional cases against quadrature
- [ ] Implement noninteracting polymer chain with orientational energy
- [ ] Implement interacting polymer chain with orientational energy; use clustering?
- [ ] Question: will clustering type algorithms work for fluid--solid transitions? There must be something in the literature
- [ ] Write descriptions for each of the examples; add page numbers when (if?) manuscript is ever published/posted
