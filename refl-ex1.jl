include(joinpath(@__DIR__, "src", "bp.jl"));

s = ArgParseSettings();
@add_arg_table! s begin
  "--C1", "-a"
    help = "coefficient 1"
    arg_type = Float64
    default = 1.0
  "--C2", "-b"
    help = "coefficient 2"
    arg_type = Float64
    default = 3.0
  "--x0", "-X"
    help = "initial configuration"
    arg_type = String
    default = "() -> 0.0"
  "--force", "-f"
    help = "applied force"
    arg_type = Float64
    default = 0.0
end

default_options = src_include("default_options.jl");
ArgParse.import_settings!(s, default_options);

pargs = src_include("parse_args.jl");

@everywhere function mcmc(nsteps::Int, pargs)
  kT = pargs["kT"];
  a = pargs["C1"];
  b = pargs["C2"];
  f = pargs["force"];
  U = (x) -> a*x^4 - b*x^2 - f*x;
  x = pargs["x0"]();
  xstep = pargs["dx"];
  dx_dist = Uniform(-xstep, xstep);
  orbf = (pargs["orbit"]) ? x -> rand([-1; 1])*x : x -> x;

  nacc = 0;
  xtotal = x;
  xrolling = Float64[];
  Ucurr = U(x);
  Utotal = Ucurr;
  Urolling = Float64[];
  x2total = x*x;
  xstd_rolling = Float64[];
  U2total = Ucurr*Ucurr;
  Ustd_rolling = Float64[];
  stepout = pargs["stepout"];
  rolls = Int[];

  start = time();
  last_update = start;
  for s = 1:nsteps
    xtrial = orbf(x + rand(dx_dist));
    Utrial = U(xtrial);
    if (Utrial < Ucurr) || (rand() <= exp(-(Utrial - Ucurr) / kT) )
      x = xtrial;
      Ucurr = Utrial;
      nacc += 1;
    end

    xtotal += x;
    Utotal += Ucurr;
    x2total += x*x;
    U2total += Ucurr*Ucurr;
    ar = nacc / s;

    if s % stepout == 0
      push!(rolls, s);
      push!(xrolling, xtotal / s);
      push!(Urolling, Utotal / s);
      push!(xstd_rolling, sqrt(max(0.0, x2total / s - (xtotal / s)^2)));
      push!(Ustd_rolling, sqrt(max(0.0, U2total / s - (Utotal / s)^2)));
    end

    if time() - last_update > pargs["update-freq"]
      @info "elapsed: $(time() - start)";
      @info "step:    $s / $nsteps";
      last_update = time();
    end

    if (
        pargs["step-adjust-scale"] != 1.0 &&
        s % pargs["steps-per-adjust"] == 0
       ) # adjust step size?
       ar = nacc / s;
       if (ar > pargs["step-adjust-ub"])
         @info "acceptance ratio is high; increasing step size";
         xstep *= pargs["step-adjust-scale"];
       elseif ar < pargs["step-adjust-lb"]
         @info "acceptance ratio is low; decreasing step size";
         xstep /= pargs["step-adjust-scale"];
       end
       dx_dist = Uniform(-xstep, xstep);
    end

  end

  return Dict(:xavg => xtotal / nsteps, :Uavg => Utotal / nsteps,
              :xrolling => xrolling, :Urolling => Urolling,
              :x2avg => x2total / nsteps, :U2avg => U2total / nsteps,
              :xstd_rolling => xstd_rolling, :Ustd_rolling => Ustd_rolling,
              :rolls => rolls, :ar => nacc / nsteps);

end

src_include("wrap_main_runs.jl");

results_std, results_polya = wrap_main_runs(pargs);
post_process_main_runs(results_std, results_polya, pargs)

if pargs["do-plots"]
  idx = rand(1:pargs["num-runs"]);
  p = plot(results_std[idx][:rolls], results_std[idx][:xrolling]; label="std mcmc",
           xlabel = "step", ylabel = "\$\\langle x \\rangle\$");
  plot!(results_polya[idx][:rolls], results_polya[idx][:xrolling]; label="polya");
  display(p);
  println();
  println("Press RETURN to exit...");
  readline();

  p = plot(results_std[idx][:rolls], results_std[idx][:Urolling]; label="std mcmc",
           xlabel = "step", ylabel = "\$\\langle U \\rangle\$");
  plot!(results_polya[idx][:rolls], results_polya[idx][:Urolling]; label="polya");
  display(p);
  println();
  println("Press RETURN to exit...");
  readline();
end
