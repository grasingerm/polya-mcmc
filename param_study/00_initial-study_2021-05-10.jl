using Distributed;
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
nruns = convert(Int, 2e1);
@show nsubprocs = if length(ARGS) > 1
  eval(Meta.parse(ARGS[2]));
else
  1;
end

#############################################################################
println("Example 1: reflection symmetry");

cases = [
         Dict(
              :a => a,
              :b => b,
              :x0 => x0func,
              :f => f,
             )
         for a in [0.2; 0.5; 1.0; 2.0; 5.0],
             b in [0.2; 0.5; 1.0; 2.0; 5.0],
             x0func in [
                        "(::Any) -> 0.0"; 
                        "(pargs) -> -sqrt(pargs[\"C2\"] / pargs[\"C1\"]) + pargs[\"force\"] / (2*pargs[\"C2\"])";
                        "(pargs) -> sqrt(pargs[\"C2\"] / pargs[\"C1\"]) + pargs[\"force\"] / (2*pargs[\"C2\"])";
                        "(::Any) -> rand(Uniform(-sqrt(pargs[\"C2\"] / pargs[\"C1\"]), sqrt(pargs[\"C2\"] / pargs[\"C1\"])))";
                        "(::Any) -> rand(Uniform(-3*sqrt(pargs[\"C2\"] / pargs[\"C1\"]), 3*sqrt(pargs[\"C2\"] / pargs[\"C1\"])))"
                       ],
             f in [0.0; 1e-2; 0.1; 0.2; 0.5; 1.0]
        ];
names = [
         "a-$(fmt(a))_b-$(fmt(b))_x0-$(x0func)_f-$(fmt(f))"
         for a in [0.2; 0.5; 1.0; 2.0; 5.0],
             b in [0.2; 0.5; 1.0; 2.0; 5.0],
             x0func in [
                        "origin"; 
                        "left-well";
                        "right-well";
                        "rand-l";
                        "rand-3l"
                       ],
             f in [0.0; 1e-2; 0.1; 0.2; 0.5; 1.0]
        ];
cases = reshape(cases, length(cases));

println("total number of cases to run for ex1: $(length(cases))");
exoutdir = joinpath(workdir, "refl-ex1");
mkpath(exoutdir);
pmap(pair -> begin;
       name, case = pair;
       local_outdir = joinpath(exoutdir, name);
       mkpath(local_outdir);
       outfile = joinpath(local_outdir, "out.txt");
       if !dryrun && !isfile(outfile)
         command = `$HOME/julia-1.6.1/bin/julia -O 3 -p $nsubprocs "refl-ex1.jl" -a $(case[:a]) -b $(case[:b]) --x0 $(case[:x0]) -f $(case[:f]) --verbose 2 --outdir $(local_outdir) --do-csvs`
         output = read(command, String);
         write(outfile, output);
       else
         @info "Case: $case has already been run.";
       end
     end, zip(names, cases));

#############################################################################
println("Example 2: translation symmetry");
cases = [
         Dict(
              :a => a,
              :n => n,
              :x0 => x0func,
              :f => f,
             )
         for a in [0.2; 0.5; 1.0; 2.0; 5.0; 10.0; 100.0; 1000.0],
             n in [2; 3; 5; 10; 20; 100; 1000],
             x0func in [
                        "(::Any) -> 0.0"; 
                        "(pargs) -> π - π/(2*pargs[\"C2\"])";
                        "(pargs) -> -π + π/(2*pargs[\"C2\"])";
                        "(pargs) -> rand(range(-π, π; step=π/(2*pargs[\"C2\"])))";
                        "(::Any) -> rand(Uniform(-π, π))"
                       ],
             f in [0.0; 1e-2; 0.1; 0.2; 0.5; 1.0]
        ];
names = [
         "a-$(fmt(a))_b-$(fmt(n))_x0-$(x0func)_f-$(fmt(f))"
         for a in [0.2; 0.5; 1.0; 2.0; 5.0; 10.0; 100.0; 1000.0],
             n in [2; 3; 5; 10; 20; 100; 1000],
             x0func in [
                        "origin"; 
                        "far-right-well";
                        "far-left-well";
                        "rand-well";
                        "rand-uniform"
                       ],
             f in [0.0; 1e-2; 0.1; 0.2; 0.5; 1.0]
        ];
cases = reshape(cases, length(cases));

println("total number of cases to run for ex2: $(length(cases))");
exoutdir = joinpath(workdir, "trans-ex2");
mkpath(exoutdir);
pmap(pair -> begin;
       name, case = pair;
       local_outdir = joinpath(exoutdir, name);
       mkpath(local_outdir);
       outfile = joinpath(local_outdir, "out.txt");
       if !dryrun && !isfile(outfile)
         command = `$HOME/julia-1.6.1/bin/julia -O 3 -p $nsubprocs "trans-ex2.jl" -a $(case[:a]) -n $(case[:n]) --x0 $(case[:x0]) -f $(case[:f]) --outdir $(local_outdir) --do-csvs --verbose 2`
         output = read(command, String);
         write(outfile, output);
       else
         @info "Case: $case has already been run.";
       end
     end, zip(names, cases));

#############################################################################
println("Example 3: D2 symmetry");
cases = [
         Dict(
              :k => k,
              :l => l,
              :x0 => x0func,
              :f => f,
              :g => g
             )
         for k in [0.2; 0.5; 1.0; 2.0; 5.0; 10.0],
             l in [0.2; 0.5; 1.0; 2.0; 5.0],
             x0func in [
                        "(::Any) -> [0.0; 0.0]"; 
                        "(::Any) -> [1/4; pargs[\"C2\"]/4]";
                        "(::Any) -> [-1/4; pargs[\"C2\"]/4]";
                        "(::Any) -> [-1/4; -pargs[\"C2\"]/4]";
                        "(::Any) -> [1/4; -pargs[\"C2\"]/4]";
                        "(::Any) -> [rand(Uniform(-1/2, 1/2)); rand(Uniform(-pargs[\"C2\"]/2, pargs[\"C2\"]/2))]"
                       ],
             f in [0.0; 1e-2; 0.1; 0.2; 0.5; 1.0],
             g in [0.0; 0.1; 1.0]
        ];
names = [
         "k-$(fmt(k))_l-$(fmt(l))_x0-$(x0func)_f-$(fmt(f))_g-$(fmt(g))"
         for k in [0.2; 0.5; 1.0; 2.0; 5.0; 10.0],
             l in [0.2; 0.5; 1.0; 2.0; 5.0],
             x0func in [
                        "origin";
                        "top-right-well";
                        "top-left-well";
                        "bottom-left-well";
                        "bottom-right-well";
                        "rand"
                       ],
             f in [0.0; 1e-2; 0.1; 0.2; 0.5; 1.0],
             g in [0.0; 0.1; 1.0]
        ];
cases = reshape(cases, length(cases));

println("total number of cases to run for ex3: $(length(cases))");
exoutdir = joinpath(workdir, "D2-ex3");
mkpath(exoutdir);
pmap(pair -> begin;
       name, case = pair;
       local_outdir = joinpath(exoutdir, name);
       mkpath(local_outdir);
       outfile = joinpath(local_outdir, "out.txt");
       if !dryrun && !isfile(outfile)
         command = `$HOME/julia-1.6.1/bin/julia -O 3 -p $nsubprocs "D2-ex3.jl" -k $(case[:k]) -l $(case[:l]) --x0 $(case[:x0]) -f "[$(case[:f]); $(case[:g])]" --verbose 2 --outdir $(local_outdir) --do-csvs`
         output = read(command, String);
         write(outfile, output);
       else
         @info "Case: $case has already been run.";
       end
     end, zip(names, cases));

#############################################################################
println("Example 4: rosenbock");
cases = [
         Dict(
              :a => a,
              :b => b,
              :x0 => x0func,
              :f => f,
              :g => g
             )
         for a in [0.0; 0.1; 1.0],
             b in [10.0; 100.0; 1000.0],
             x0func in [
                        "(::Any) -> [0.0; 0.0]"; 
                        "(::Any) -> rand(Normal(0.0, 1.0), 2)";
                        "(::Any) -> rand(Normal(0.0, 3.0), 2)"
                       ],
             f in [0.0; 1e-2; 0.1; 0.2; 0.5; 1.0],
             g in [0.0; 1e-2; 0.1; 0.2; 0.5; 1.0]
        ];
names = [
         "a-$(fmt(a))_b-$(fmt(b))_x0-$(x0func)_f-$(fmt(f))_g-$(fmt(g))"
         for a in [0.0; 0.1; 1.0],
             b in [10.0; 100.0; 1000.0],
             x0func in [
                        "origin"; 
                        "normal-origin-1";
                        "normal-origin-3"
                       ],
             f in [0.0; 1e-2; 0.1; 0.2; 0.5; 1.0],
             g in [0.0; 1e-2; 0.1; 0.2; 0.5; 1.0]
        ];
cases = reshape(cases, length(cases));

println("total number of cases to run for ex4: $(length(cases))");
exoutdir = joinpath(workdir, "rosenbrock-ex4");
mkpath(exoutdir);
pmap(pair -> begin;
       name, case = pair;
       local_outdir = joinpath(exoutdir, name);
       mkpath(local_outdir);
       outfile = joinpath(local_outdir, "out.txt");
       if !dryrun && !isfile(outfile)
         command = `$HOME/julia-1.6.1/bin/julia -O 3 -p $nsubprocs "rosenbrock-ex4.jl" -a $(case[:a]) -b $(case[:b]) --x0 $(case[:x0]) -f "[$(case[:f]); $(case[:g])]" --verbose 2 --outdir $(local_outdir) --do-csvs`
         output = read(command, String);
         write(outfile, output);
       else
         @info "Case: $case has already been run.";
       end
     end, zip(names, cases));
