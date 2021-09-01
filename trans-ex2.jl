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
  wt = pargs["weight"];
  U = (x) -> a*cos(n*x) - f*x;
  x = pargs["x0"](pargs);
  xstep = pargs["dx"];
  dx_dist = Uniform(-xstep, xstep);
  orbf = (pargs["orbit"]) ? x -> x + rand([-1; 0; 1])*(2*π / n) : x -> x;

  nacc = 0;
  wt_curr = wt(x);
  inv_wt_total = 1 / wt_curr;
  xtotal = x / wt_curr;
  xrolling = Float64[];
  Ucurr = U(x);
  Utotal = Ucurr / wt_curr;
  Urolling = Float64[];
  x2total = x*x / wt_curr;
  x2rolling = Float64[];
  xstd_rolling = Float64[];
  U2total = Ucurr*Ucurr / wt_curr;
  U2rolling = Float64[];
  Ustd_rolling = Float64[];
  stepout = pargs["stepout"];
  rolls = Int[];

  start = time();
  last_update = start;
  for s = 1:nsteps
    xtrial = orbf(x + rand(dx_dist));
    xtrial = if xtrial > π # apply periodic boundary conditions
      -π + (xtrial - π);
    elseif xtrial < -π
      π - (-π - xtrial);
    else
      xtrial;
    end
    Utrial = U(xtrial);
    wt_trial = wt(xtrial);
    if rand() <= ( exp(-(Utrial - Ucurr) / kT) * wt_trial / wt_curr )
      x = xtrial;
      Ucurr = Utrial;
      wt_curr = wt_trial;
      nacc += 1;
    end

    inv_wt_total += 1 / wt_curr;
    xtotal += x / wt_curr;
    Utotal += Ucurr / wt_curr;
    x2total += x*x / wt_curr;
    U2total += Ucurr*Ucurr / wt_curr;
    ar = nacc / s;

    if s % stepout == 0
      push!(rolls, s);
      push!(xrolling, xtotal / inv_wt_total);
      push!(Urolling, Utotal / inv_wt_total);
      push!(x2rolling, x2total / inv_wt_total);
      push!(U2rolling, U2total / inv_wt_total);
      push!(xstd_rolling, sqrt(max(0.0, x2total / inv_wt_total - 
                          (xtotal / inv_wt_total)^2)));
      push!(Ustd_rolling, sqrt(max(0.0, U2total / inv_wt_total - 
                          (Utotal / inv_wt_total)^2)));
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

  return Dict(:xavg => xtotal / inv_wt_total, 
              :Uavg => Utotal / inv_wt_total,
              :xrolling => xrolling, 
              :Urolling => Urolling,
              :x2avg => x2total / inv_wt_total,
              :U2avg => U2total / inv_wt_total,
              :x2rolling => x2rolling, 
              :U2rolling => U2rolling,
              :xstd_rolling => xstd_rolling, 
              :Ustd_rolling => Ustd_rolling,
              :rolls => rolls, :ar => nacc / nsteps);

end

maxevals = convert(Int, 1e6);
kT = pargs["kT"];
a = pargs["C1"];
n = pargs["C2"];
f = pargs["force"];
U = (x) -> a*cos(n*x) - f*x;
Zquad, Zquad_err = hquadrature(x -> exp(-U(x)/kT), -π, π; maxevals=maxevals);
xquad, xquad_err = hquadrature(x -> x * exp(-U(x)/kT) / Zquad, -π, π; maxevals=maxevals);
Uquad, Uquad_err = hquadrature(x -> U(x) * exp(-U(x)/kT) / Zquad, -π, π; maxevals=maxevals);
x2quad, x2quad_err = hquadrature(x -> x*x * exp(-U(x)/kT) / Zquad, -π, π; maxevals=maxevals);
U2quad, U2quad_err = hquadrature(x -> (U(x))^2 * exp(-U(x)/kT) / Zquad, -π, π; maxevals=maxevals);

src_include("wrap_main_runs.jl");
results_std, results_polya = wrap_main_runs(pargs);
L1_error_std, L1_error_polya = post_process_main_runs(results_std,
                                                      results_polya,
                                                      pargs;
                                                      xrolling=xquad,
                                                      Urolling=Uquad);
μs = vcat(range(0.0, -π; step=-2*π/n), range(0.0, π; step=2*π/n)[2:end]);
pargs["weight"] = x -> a * sum(map(μ -> exp(-n^2*(x-μ)^2 / 2), μs));
pargs["outdir"] = joinpath(pargs["outdir"], "umb");
results_umb_std, results_umb_polya = wrap_main_runs(pargs);
L1_error_umb_std, L1_error_umb_polya = post_process_main_runs(results_umb_std, 
                                                              results_umb_polya, 
                                                              pargs);

idx = rand(1:pargs["num-runs"]);
@info "Z via quadrature = $Zquad; error = $Zquad_err";
@info "<x> via quadrature = $xquad; error = $xquad_err";
@info "<x> via std = $(results_std[idx][:xavg])";
@info "<x> via polya = $(results_polya[idx][:xavg])";
@info "<x> via umb std = $(results_umb_std[idx][:xavg])";
@info "<x> via umb polya = $(results_umb_polya[idx][:xavg])";
@info "<U> via quadrature = $Uquad; error = $Uquad_err";
@info "<U> via std = $(results_std[idx][:Uavg])";
@info "<U> via polya = $(results_polya[idx][:Uavg])";
@info "<U> via umb std = $(results_umb_std[idx][:Uavg])";
@info "<U> via umb polya = $(results_umb_polya[idx][:Uavg])";
@info "<x^2> via quadrature = $x2quad; error = $x2quad_err";
@info "<x^2> via std = $(results_std[idx][:x2avg])";
@info "<x^2> via polya = $(results_polya[idx][:x2avg])";
@info "<x^2> via umb std = $(results_umb_std[idx][:x2avg])";
@info "<x^2> via umb polya = $(results_umb_polya[idx][:x2avg])";

if pargs["do-plots"]
  idx = rand(1:pargs["num-runs"]);
  p = plot(results_std[idx][:rolls], results_std[idx][:xrolling]; label="std mcmc",
           xlabel = "step", ylabel = "\$\\langle x \\rangle\$");
  plot!(results_polya[idx][:rolls], results_polya[idx][:xrolling]; label="polya");
  plot!(results_umb_std[idx][:rolls], results_umb_std[idx][:xrolling]; label="umb");
  plot!(results_umb_polya[idx][:rolls], results_umb_polya[idx][:xrolling]; label="umb + polya");
  display(p);
  println();
  println("Press RETURN to exit...");
  readline();

  p = plot(results_std[idx][:rolls], results_std[idx][:Urolling]; label="std mcmc",
           xlabel = "step", ylabel = "\$\\langle U \\rangle\$");
  plot!(results_polya[idx][:rolls], results_polya[idx][:Urolling]; label="polya");
  plot!(results_umb_std[idx][:rolls], results_umb_std[idx][:Urolling]; label="umb");
  plot!(results_umb_polya[idx][:rolls], results_umb_polya[idx][:Urolling]; label="umb + polya");
  display(p);
  println();
  println("Press RETURN to exit...");
  readline();

  p = plot(results_std[idx][:rolls], results_std[idx][:x2rolling]; label="std mcmc",
           xlabel = "step", ylabel = "\$\\langle x^2 \\rangle\$");
  plot!(results_polya[idx][:rolls], results_polya[idx][:x2rolling]; label="polya");
  plot!(results_umb_std[idx][:rolls], results_umb_std[idx][:x2rolling]; label="umb");
  plot!(results_umb_polya[idx][:rolls], results_umb_polya[idx][:x2rolling]; label="umb + polya");
  display(p);
  println();
  println("Press RETURN to exit...");
  readline();

end
