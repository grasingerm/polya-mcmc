import Pkg;

for name in [
             "ArgParse",
             "Distributions",
             "ProfileView",
             "Plots",
             "DecFP",
             "Quadmath",
             "Logging",
             "Cubature",
             "Logging",
             "ParallelDataTransfer"
            ]

  Pkg.add(name);

end

Pkg.update();
