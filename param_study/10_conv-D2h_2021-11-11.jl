using Distributed;
using DelimitedFiles;
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
println("total number of cases to run for ex1: $(length(cases))");
mkpath(workdir);
open(joinpath(workdir, "D2h-conv.csv"), "w") do w
  rows = pmap(case -> begin;
           command = `julia -O 3 -p $nsubprocs "ex6-D2h.jl" -R $nruns -N $nsteps -a 1.0 -b 1.0 -c 1.0 --kT 1.0 -f "[$(case[:f]),$(case[:f]),$(case[:f])] / sqrt(3)" -q $(case[:q]) --verbose 2 --do-conv-rates`
           output = read(command, String);
           max_αs = Dict();
           for line in split(output, '\n')
             if contains(line, '=')
               sides = map(strip, split(line, '='; limit=2));
               if (
                   startswith(sides[1], "αs") &&
                   maximum(map(k -> contains(sides[1], k), ks))
                  )
                 αs = eval(Meta.parse(sides[2]));
                 max_αs[split(sides[1], ' ')[end]] = maximum(values(αs));
               end
             end
           end
           return hcat([case[:q] case[:f]], transpose(map(k -> max_αs[k], ks)));
       end, cases);
  writedlm(w, vcat(hcat(["q" "f"], permutedims(ks)), rows...), ',');
end
