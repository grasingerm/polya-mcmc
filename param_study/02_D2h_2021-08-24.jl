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
println("Example 6: hard disks in a potential with reflection symmetry");

cases = [
         Dict(
              :a => a,
              :b => b,
              :c => c,
              :q => q,
              :f => f,
             )
         for a in [0.2; 0.5; 1.0; 2.0; 5.0],
             b in [0.5; 1.0; 2.0],
             c in [1.0; 2.0],
             q in [-3.0; -1.0; -0.1; 0.1; 1.0; 3.0],
             f in [0.0; 1e-2; 0.1; 0.2; 0.5; 1.0]
        ];
names = [
         "a-$(fmt(a))_b-$(fmt(b))_c-$(fmt(c))_q-$(fmt(q))_f-$(fmt(f))"
         for a in [0.2; 0.5; 1.0; 2.0; 5.0],
             b in [0.5; 1.0; 2.0],
             c in [1.0; 2.0],
             q in [-3.0; -1.0; -0.1; 0.1; 1.0; 3.0],
             f in [0.0; 1e-2; 0.1; 0.2; 0.5; 1.0]
        ];
cases = reshape(cases, length(cases));

println("total number of cases to run for ex6: $(length(cases))");
exoutdir = joinpath(workdir, "D2h-ex6");
mkpath(exoutdir);
map(pair -> begin;
       name, case = pair;
       local_outdir = joinpath(exoutdir, name);
       mkpath(local_outdir);
       outfile = joinpath(local_outdir, "out.txt");
       if !dryrun && !isfile(outfile)
         command = `julia -O 3 -p $nsubprocs D2h-ex6.jl -a $(case[:a]) -b $(case[:b]) -c $(case[:c]) --charge $(case[:q]) -f "[$(case[:f]); 0.0; 0.0]" --verbose 2 --outdir $(local_outdir) --do-csvs`
         output = read(command, String);
         write(outfile, output);
       else
         @info "Case: $case has already been run.";
       end
     end, zip(names, cases));
