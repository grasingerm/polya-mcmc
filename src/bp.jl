using Distributed;
if false
  using ProfileView;
end
@everywhere using DelimitedFiles;
@everywhere using Plots;
@everywhere using ArgParse;
@everywhere using ParallelDataTransfer;
@everywhere using Distributions;
@everywhere using LinearAlgebra;
@everywhere using Logging;
@everywhere using DecFP;
@everywhere using Quadmath;

src_include(fname) = include(joinpath(@__DIR__, fname));
src_evalfile(fname) = evalfile(joinpath(@__DIR__, fname));
