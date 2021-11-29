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
nruns = convert(Int, 1e2);
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
              :f => f,
              :kT => kT
             )
         for a in [0.5; 1.0],
             b in [4.0; 8.0; 16.0],
             kT in [1.0; 0.1],
             f in [0.0; 1.0; 4.0; 8.0]
        ];
names = [
         "a-$(fmt(a))_b-$(fmt(b))_kT-$(kT)_f-$(fmt(f))"
         for a in [0.5; 1.0],
             b in [4.0; 8.0; 16.0],
             kT in [1.0; 0.1],
             f in [0.0; 1.0; 4.0; 8.0]
        ];
cases = reshape(cases, length(cases));

println("total number of cases to run for ex1: $(length(cases))");
exoutdir = joinpath(workdir, "ex1-refl");
mkpath(exoutdir);
pmap(pair -> begin;
       name, case = pair;
       local_outdir = joinpath(exoutdir, name);
       mkpath(local_outdir);
       outfile = joinpath(local_outdir, "out.txt");
       if !dryrun && !isfile(outfile)
         command = `julia -O 3 -p $nsubprocs "ex1-refl.jl" -R $nruns -N $nsteps -a $(case[:a]) -b $(case[:b]) --kT $(case[:kT]) -f $(case[:f]) --verbose 2 --outdir $(local_outdir) --do-csvs`
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
              :kT => kT,
              :f => f,
             )
         for a in [1.0; 5.0; 10.0;],
             n in [2; 3; 5; 10],
             kT in [1.0; 0.1],
             f in [0.0; 1.0; 4.0; 8.0]
        ];
names = [
         "a-$(fmt(a))_n-$(fmt(n))_kT-$(kT)_f-$(fmt(f))"
         for a in [1.0; 5.0; 10.0;],
             n in [2; 3; 5; 10],
             kT in [1.0; 0.1],
             f in [0.0; 1.0; 4.0; 8.0]
        ];

println("total number of cases to run for ex2: $(length(cases))");
exoutdir = joinpath(workdir, "ex2-trans");
mkpath(exoutdir);
pmap(pair -> begin;
       name, case = pair;
       local_outdir = joinpath(exoutdir, name);
       mkpath(local_outdir);
       outfile = joinpath(local_outdir, "out.txt");
       if !dryrun && !isfile(outfile)
         command = `julia -O 3 -p $nsubprocs "ex2-trans.jl" -R $nruns -N $nsteps -a $(case[:a]) -n $(case[:n]) --kT $(case[:kT]) -f $(case[:f]) --outdir $(local_outdir) --do-csvs --verbose 2`
         output = read(command, String);
         write(outfile, output);
       else
         @info "Case: $case has already been run.";
       end
     end, zip(names, cases));

#############################################################################
println("Example 6: D2h");

cases = [
         Dict(
              :a => a,
              :b => b,
              :c => c,
              :q => q,
              :kT => kT,
              :f => f,
             )
         for a in [1.0; 2.0; 5.0],
             b in [1.0; 2.0],
             c in [1.0;],
             q in [-3.0; -1.0; 1.0; 3.0],
             kT in [1.0; 0.1; 0.05],
             f in [0.0; 1.0; 4.0]
        ];
names = [
         "a-$(fmt(a))_b-$(fmt(b))_c-$(fmt(c))_q-$(fmt(q))_kT-$(fmt(kT))_f-$(fmt(f))"
         for a in [1.0; 2.0; 5.0],
             b in [1.0; 2.0],
             c in [1.0;],
             q in [-3.0; -1.0; 1.0; 3.0],
             kT in [1.0; 0.1; 0.05],
             f in [0.0; 1.0; 4.0]
        ];
cases = reshape(cases, length(cases));

println("total number of cases to run for ex6: $(length(cases))");
exoutdir = joinpath(workdir, "ex6-D2h");
mkpath(exoutdir);
map(pair -> begin;
       name, case = pair;
       local_outdir = joinpath(exoutdir, name);
       mkpath(local_outdir);
       outfile = joinpath(local_outdir, "out.txt");
       if !dryrun && !isfile(outfile)
         command = `julia -O 3 -p $nsubprocs ex6-D2h.jl --x0 "(args) -> rand(Uniform(-1, 1), 3)" -R $nruns -N $nsteps -a $(case[:a]) -b $(case[:b]) -c $(case[:c]) --charge $(case[:q]) -f "[$(case[:f]); 0.0; 0.0]" --verbose 2 --outdir $(local_outdir) --do-csvs`
         output = read(command, String);
         write(outfile, output);
       else
         @info "Case: $case has already been run.";
       end
     end, zip(names, cases));

#############################################################################
println("Example 7: D2h-LJ");

cases = [
         Dict(
              :a => a,
              :b => b,
              :c => c,
              :eps => eps,
              :sigma => sigma,
              :kT => kT,
              :f => f,
             )
         for a in [1.0; 2.0; 5.0],
             b in [1.0; 2.0],
             c in [1.0;],
             eps in [0.1; 1.0; 10.0],
             sigma in [0.05; 0.1; 0.25; 0.5],
             kT in [1.0; 0.1; 0.05],
             f in [0.0; 1.0; 4.0]
        ];
names = [
         "a-$(fmt(a))_b-$(fmt(b))_c-$(fmt(c))_eps-$(fmt(eps))_sigma-$(fmt(sigma))_kT-$(fmt(kT))_f-$(fmt(f))"
         for a in [1.0; 2.0; 5.0],
             b in [1.0; 2.0],
             c in [1.0;],
             eps in [0.1; 1.0; 10.0],
             sigma in [0.05; 0.1; 0.25; 0.5],
             kT in [1.0; 0.1; 0.05],
             f in [0.0; 1.0; 4.0]
        ];
cases = reshape(cases, length(cases));

println("total number of cases to run for ex7: $(length(cases))");
exoutdir = joinpath(workdir, "D2h-ex7");
mkpath(exoutdir);
map(pair -> begin;
       name, case = pair;
       local_outdir = joinpath(exoutdir, name);
       mkpath(local_outdir);
       outfile = joinpath(local_outdir, "out.txt");
       if !dryrun && !isfile(outfile)
         command = `julia -O 3 -p $nsubprocs ex7-D2h-LJ.jl -R $nruns -N $nsteps -a $(case[:a]) -b $(case[:b]) -c $(case[:c]) --eps $(case[:eps]) --sigma $(case[:sigma]) -f "[$(case[:f]); 0.0; 0.0]" --verbose 2 --outdir $(local_outdir) --do-csvs`
         output = read(command, String);
         write(outfile, output);
       else
         @info "Case: $case has already been run.";
       end
     end, zip(names, cases));
