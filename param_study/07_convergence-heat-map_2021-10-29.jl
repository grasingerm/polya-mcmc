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
nsteps = convert(Int, 1e5);
nruns = convert(Int, 1e1);
@show nsubprocs = if length(ARGS) > 1
  eval(Meta.parse(ARGS[2]));
else
  1;
end

#############################################################################
println("Example 1: reflection symmetry");

cases = [
         Dict(
              :b => b,
              :f => f,
             )
         for b in 1.0:1.0:32.0,
             f in 0.0:0.5:15.0
        ];
cases = reshape(cases, length(cases));

ks = ["std", "umb", "polya", "gu"];
println("total number of cases to run for ex1: $(length(cases))");
mkpath(workdir);
open(joinpath(workdir, "refl-conv.csv"), "w") do w
  rows = pmap(case -> begin;
           command = `julia -O 3 -p $nsubprocs "refl-ex1.jl" -R $nruns -N $nsteps -a 1.0 -b $(case[:b]) --kT 1.0 -f $(case[:f]) --verbose 2 --do-conv-rates`
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
           @show max_αs
           @show map(k -> max_αs[k], ks)
           return hcat([case[:b] case[:f]], transpose(map(k -> max_αs[k], ks)));
       end, cases);
  writedlm(w, vcat(hcat(["b" "f"], ks), rows...), ',');
end

#############################################################################
println("Example 2: translation symmetry");
cases = [
         Dict(
              :n => n,
              :f => f,
             )
         for n in 2:32,
             f in 0.0:0.5:15.0
        ];
cases = reshape(cases, length(cases));

println("total number of cases to run for ex2: $(length(cases))");
open(joinpath(workdir, "trans-conv.csv"), "w") do w
  rows = pmap(case -> begin;
         command = `julia -O 3 -p $nsubprocs "trans-ex2.jl" -R $nruns -N $nsteps -a 1.0 -n $(case[:n]) --kT 1.0 -f $(case[:f]) --verbose 2 --do-conv-rates`
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
         return hcat([case[:n] case[:f]], transpose(map(k -> max_αs[k], ks)));
       end, cases);
  writedlm(w, vcat(hcat(["n" "f"], ks), rows...), ',');
end
