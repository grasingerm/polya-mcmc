include(joinpath(@__DIR__, "src", "bp.jl"));
using Cubature;

s = ArgParseSettings();
@add_arg_table! s begin
  "--x0", "-X"
    help = "initial configuration (Function, i.e. rand)"
    arg_type = String
    default = "(::Any) -> 0.0"
  "--C1", "-a"
    help = "coefficient 1"
    arg_type = Float64
    default = 1.0
  "--C2", "-n"
    help = "coefficient 2"
    arg_type = Int
    default = 4
  "--force", "-f"
    help = "applied force"
    arg_type = Float64
    default = 0.0
end

default_options = src_include("default_options.jl");
ArgParse.import_settings!(s, default_options);
pargs = src_include("parse_args.jl");

@everywhere function mcmc(nsteps::Int, pargs::Dict)
  kT = pargs["kT"];
  a = pargs["C1"];
  n = pargs["C2"];
  f = pargs["force"];
  U = (x) -> a*cos(n*x) - f*x;
  x = pargs["x0"](pargs);
  xstep = pargs["dx"];
  dx_dist = Uniform(-xstep, xstep);
  orbf = (pargs["orbit"]) ? x -> x + rand([-1; 0; 1])*(2*π / n) : x -> x;

  nacc = 0;
  xtotal = x;
  xrolling = Float64[];
  Ucurr = U(x);
  Utotal = Ucurr;
  Urolling = Float64[];
  x2total = x*x;
  x2rolling = Float64[];
  xstd_rolling = Float64[];
  U2total = Ucurr*Ucurr;
  U2rolling = Float64[];
  Ustd_rolling = Float64[];
  stepout = pargs["stepout"];
  rolls = Int[];

  start = time();
  last_update = start;
  for s = 1:nsteps
    xtrial = max(-π, min(orbf(x + rand(dx_dist)), π));
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
      push!(x2rolling, x2total / s);
      push!(U2rolling, U2total / s);
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
              :x2rolling => x2rolling, :Urolling => U2rolling,
              :xstd_rolling => xstd_rolling, :Ustd_rolling => Ustd_rolling,
              :rolls => rolls, :ar => nacc / nsteps);

end

maxevals = convert(Int, 1e6);
kT = pargs["kT"];
a = pargs["C1"];
n = pargs["C2"];
f = pargs["force"];
U = (x) -> a*cos(n*x) - f*x;
Z_quad, err = pquadrature(x -> exp(-U(x)/kT), -π, π; maxevals=maxevals);
x_quad, err = pquadrature(x -> x * exp(-U(x)/kT) / Z_quad, -π, π; maxevals=maxevals);
U_quad, err = pquadrature(x -> U(x) * exp(-U(x)/kT) / Z_quad, -π, π; maxevals=maxevals);
x2_quad, err = pquadrature(x -> x*x * exp(-U(x)/kT) / Z_quad, -π, π; maxevals=maxevals);
U2_quad, err = pquadrature(x -> (U(x))^2 * exp(-U(x)/kT) / Z_quad, -π, π; maxevals=maxevals);

@info Z_quad, x_quad, U_quad, x2_quad, U2_quad;

src_include("wrap_main_runs.jl");
results_std, results_polya = wrap_main_runs(pargs);
L1_error_std, L1_error_polya = post_process_main_runs(results_std,
                                                      results_polya,
                                                      pargs;
                                                      xrolling=x_quad,
                                                      Urolling=U_quad);

if pargs["do-plots"]
  p = plot(results_std[1][:rolls], L1_error_std[:xrolling]; label="std mcmc",
           xlabel = "step", 
           ylabel = "\\left|\\left|\$\\langle x \\rangle - x_q\\right|\\right|\$");
  plot!(results_polya[1][:rolls], L1_error_polya[:xrolling]; label="polya");
  savefig(p, joinpath(pargs["outdir"], "x_convergence.pdf"));
  display(p);
  println();
  println("Press RETURN to exit...");
  readline();

  p = plot(results_std[1][:rolls], L1_error_std[:Urolling]; label="std mcmc",
           xlabel = "step", 
           ylabel = "\\left|\\left|\$\\langle x \\rangle - x_q\\right|\\right|\$");
  plot!(results_polya[1][:rolls], L1_error_polya[:Urolling]; label="polya");
  savefig(p, joinpath(pargs["outdir"], "U_convergence.pdf"));
  display(p);
  println();
  println("Press RETURN to exit...");
  readline();
end
