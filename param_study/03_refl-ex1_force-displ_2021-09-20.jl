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

exoutdir = workdir;
mkpath(exoutdir);
for f in vcat(collect(0.001:0.001:0.005), collect(0.01:0.01:0.05),
              collect(0.1:0.1:0.5), collect(1.0:1.0:5.0),
              collect(10.0:10.0:50.0))
  local_outdir = joinpath(exoutdir, "f-$(fmt(f))");
  mkpath(local_outdir);
  outfile = joinpath(local_outdir, "out.txt");
  if !dryrun && !isfile(outfile)
    command = `julia -O 3 -p $nsubprocs "refl-ex1.jl" -a 1.0 -b 8.0 --x0 "(::Any) -> 0.0" -f $f --umbrella-b 2.0 --dx 1.0 -R 100 -N 1000000 --stepout 500 --kT 1.0 --verbose 3 --outdir $(local_outdir) --do-csvs`
    output = read(command, String);
    write(outfile, output);
  else
    @info "$f already run.";
  end
end
