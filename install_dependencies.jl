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
             "ParallelDataTransfer",
             "GLM",
             "DataFrames"
            ]

  Pkg.add(name);

end

Pkg.update();
