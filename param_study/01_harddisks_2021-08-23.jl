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
println("Example 5: hard disks in a potential with reflection symmetry");

cases = [
         Dict(
              :a => a,
              :b => b,
              :nd => num_disks,
              :dd => disk_diam,
              :f => f,
             )
         for a in [0.2; 0.5; 1.0; 2.0; 5.0],
             b in [0.2; 0.5; 1.0; 2.0; 5.0],
             num_disks in [4; 10; 100],
             disk_diam in [0.01; 0.1; 1.0],
             f in [0.0; 1e-2; 0.1; 0.2; 0.5; 1.0]
        ];
names = [
         "a-$(fmt(a))_b-$(fmt(b))_f-$(fmt(f))_num-disks-$(num_disks)_disk-diam_$(fmt(disk_diam))"
         for a in [0.2; 0.5; 1.0; 2.0; 5.0],
             b in [0.2; 0.5; 1.0; 2.0; 5.0],
             num_disks in [4; 10; 100],
             disk_diam in [0.01; 0.1; 1.0],
             f in [0.0; 1e-2; 0.1; 0.2; 0.5; 1.0]
        ];
cases = reshape(cases, length(cases));

println("total number of cases to run for ex5: $(length(cases))");
exoutdir = joinpath(workdir, "harddisks-refl-ex5");
mkpath(exoutdir);
map(pair -> begin;
       name, case = pair;
       local_outdir = joinpath(exoutdir, name);
       mkpath(local_outdir);
       outfile = joinpath(local_outdir, "out.txt");
       if !dryrun && !isfile(outfile)
         command = `julia -O 3 -p $nsubprocs harddisks-refl-ex5.jl -a $(case[:a]) -b $(case[:b]) --num-disks $(case[:nd]) --disk-diameter $(case[:dd]) -f $(case[:f]) --verbose 2 --outdir $(local_outdir) --do-csvs`
         output = read(command, String);
         write(outfile, output);
       else
         @info "Case: $case has already been run.";
       end
     end, zip(names, cases));
