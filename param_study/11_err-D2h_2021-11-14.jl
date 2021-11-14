using Distributed;
using DelimitedFiles;
using Glob;
@everywhere using Printf;

@show workdir = if length(ARGS) > 0
  ARGS[1];
else
  ".";
end

HOME = "/home/mgrasing"
dryrun = false;
fmt(x) = @sprintf("%06d", round(Int, 1e3*x));
nsteps = convert(Int, 1e6);
nruns = convert(Int, 1e2);
@show nsubprocs = if length(ARGS) > 1
  eval(Meta.parse(ARGS[2]));
else
  1;
end

#############################################################################
println("Example 6: D2h");

cases = [
         Dict(
              :q => q,
              :f => f,
             )
         for q in 1.0:0.5:10.0,
             f in 0.0:0.5:9.0
        ];
cases = reshape(cases, length(cases));

ks = ["std", "umb", "polya", "gu"];
colidxs = Dict("std" => 2, "polya" => 3, "umb" => 4, "gu" => 5);
println("total number of cases to run for ex1: $(length(cases))");
mkpath(workdir);
open(joinpath(workdir, "D2h-err.csv"), "w") do w
  rows = pmap(case -> begin;
           subdir = "q-$(round(Int, case[:q]*1000))_f-$(round(Int, case[:f]*1000))";
           mkpath(joinpath(workdir, subdir));
           command = `julia -O 3 -p $nsubprocs "D2h-ex6.jl" -R $nruns -N $nsteps -a 1.0 -b 1.0 -c 1.0 --kT 1.0 -f "[$(case[:f]),$(case[:f]),$(case[:f])] / sqrt(3)" -q $(case[:q]) --verbose 2 --do-conv-rates --do-csvs --outdir $(joinpath(workdir, subdir))`
           output = read(command, String);
           errs = Dict();
           for datafile in readdir(glob"L1*.csv", joinpath(workdir, subdir))
             data = readdlm(datafile, ',');
             for (k, v) in colidxs
               if haskey(errs, k)
                 push!(errs[k], data[end, v]);
               else
                 errs[k] = [data[end, v]];
               end
             end
           end
           return hcat([case[:q] case[:f]], transpose(map(k -> maximum(errs[k]), ks)));
       end, cases);
  writedlm(w, vcat(hcat(["q" "f"], permutedims(ks)), rows...), ',');
end
