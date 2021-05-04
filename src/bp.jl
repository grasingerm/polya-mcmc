using ArgParse;
using ProfileView;
using DelimitedFiles;
using Plots;
using Distributed;
@everywhere using Distributions;
@everywhere using LinearAlgebra;
@everywhere using Logging;
@everywhere using DecFP;
@everywhere using Quadmath;

src_include(fname) = include(joinpath(@__DIR__, fname));
src_evalfile(fname) = evalfile(joinpath(@__DIR__, fname));
