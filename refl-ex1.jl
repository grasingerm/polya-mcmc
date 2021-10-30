include(joinpath(@__DIR__, "src", "bp.jl"));

using Cubature;

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
    default = "(::Any) -> 0.0"
  "--force", "-f"
    help = "applied force"
    arg_type = Float64
    default = 0.0
  "--umbrella-b"
    help = "pseudopotential 'b'"
    arg_type = Float64
    default = 0.5
  "--umbrella-C"
    help = "pseudopotential 'C'"
    arg_type = Float64
    default = 1.0
end

default_options = src_include("default_options.jl");
ArgParse.import_settings!(s, default_options);

pargs = src_include("parse_args.jl");

@everywhere function mcmc(nsteps::Int, pargs)
  kT = pargs["kT"];
  a = pargs["C1"];
  b = pargs["C2"];
  f = pargs["force"];
  wt = pargs["weight"];
  U = (x) -> a*x^4 - b*x^2 - f*x;
  x = pargs["x0"](pargs);
  xstep = pargs["dx"];
  dx_dist = Uniform(-xstep, xstep);
  orbf = (pargs["orbit"]) ? x -> rand([-1; 1])*x : x -> x;

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

  traj_file = if pargs["do-trajectory"]
    mkpath(pargs["outdir"]);
    open(joinpath(pargs["outdir"], "trajectory.csv"), "w");
  else
    nothing;
  end

  natt = 0;
  nacc = 0;
  nacc_total = 0;

  start = time();
  last_update = start;
  for s = 1:nsteps
    xtrial = orbf(x + rand(dx_dist));
    Utrial = U(xtrial);
    wt_trial = wt(xtrial);
    if rand() <= ( exp(-(Utrial - Ucurr) / kT) * wt_trial / wt_curr )
      x = xtrial;
      Ucurr = Utrial;
      wt_curr = wt_trial;
      nacc += 1;
      nacc_total += 1;
    end
    natt += 1;

    inv_wt_total += 1 / wt_curr;
    xtotal += x / wt_curr;
    Utotal += Ucurr / wt_curr;
    x2total += x*x / wt_curr;
    U2total += Ucurr*Ucurr / wt_curr;
    ar = nacc_total / s; # rolling acceptance ratio

    if s % stepout == 0
      push!(rolls, s);
      push!(xrolling, xtotal / inv_wt_total);
      push!(Urolling, Utotal / inv_wt_total);
      push!(x2rolling, x2total / inv_wt_total);
      push!(U2rolling, U2total / inv_wt_total);
      push!(xstd_rolling, sqrt(max(0.0, x2total / inv_wt_total 
                                        - (xtotal / inv_wt_total)^2)));
      push!(Ustd_rolling, sqrt(max(0.0, U2total / inv_wt_total - 
                                        (Utotal / inv_wt_total)^2)));
      if pargs["do-trajectory"]
        writedlm(traj_file, [x Ucurr x*x Ucurr*Ucurr ar s]);
      end
    end

    if time() - last_update > pargs["update-freq"]
      @info "elapsed: $(time() - start)";
      @info "step:    $s / $nsteps";
      last_update = time();
    end

    if (
        pargs["step-adjust-scale"] != 1.0   &&
        pargs["steps-per-adjust"] != 0      &&
        s % pargs["steps-per-adjust"] == 0
       ) # adjust step size?
       ar = nacc / natt;
       if (ar > pargs["step-adjust-ub"])
         nacc = 0;
         natt = 0;
         xstep *= pargs["step-adjust-scale"];
         @info "acceptance ratio is high; increasing step size", xstep;
       elseif ar < pargs["step-adjust-lb"]
         nacc = 0;
         natt = 0;
         xstep /= pargs["step-adjust-scale"];
         @info "acceptance ratio is low; decreasing step size", xstep;
       end
       dx_dist = Uniform(-xstep, xstep);
    end

  end

  if pargs["do-trajectory"]
    close(traj_file);
  end

  return Dict(:xavg => xtotal / inv_wt_total, 
              :Uavg => Utotal / inv_wt_total,
              :xrolling => copy(xrolling), :Urolling => copy(Urolling),
              :x2avg => x2total / inv_wt_total, 
              :U2avg => U2total / inv_wt_total,
              :x2rolling => copy(x2rolling), :U2rolling => copy(U2rolling),
              :xstd_rolling => copy(xstd_rolling), 
              :Ustd_rolling => copy(Ustd_rolling),
              :rolls => rolls, 
              :ar => nacc / natt,
              :AR => nacc_total / nsteps
             );

end

src_include("wrap_main_runs.jl");

results_std, results_polya = wrap_main_runs(pargs);

uC, ub = pargs["umbrella-C"], pargs["umbrella-b"];
pargs_umb = copy(pargs);
pargs_umb["weight"] = x -> uC*exp(-ub*x^2);
pargs_umb["outdir"] = joinpath(pargs["outdir"], "umb");
results_umb_std, results_umb_polya = wrap_main_runs(pargs_umb);

kT = pargs["kT"];
a = pargs["C1"];
b = pargs["C2"];
f = pargs["force"];
wt = pargs["weight"];
U = (x) -> a*x^4 - b*x^2 - f*x;
meq = 25000;
ℓ = 4*sqrt(b / a);
(Zquad, Zquad_err) = (
                      collect(hquadrature(t -> -exp(-U(1/t)/kT) / t^2, 0.0, -1/ℓ; maxevals=meq)) +
                      collect(hquadrature(x -> exp(-U(x)/kT), -ℓ, ℓ; maxevals=meq)) +
                      collect(hquadrature(t -> -exp(-U(1/t)/kT) / t^2, 1/ℓ, 0.0; maxevals=meq))
                     );
(xquad, xquad_err) = (f == 0) ? [0.0; 0.0] : (
                       (
                        collect(hquadrature(t -> -(1/t)*exp(-U(1/t)/kT) / t^2 / Zquad, 0.0, -1/ℓ; maxevals=meq)) +
                        collect(hquadrature(x -> x*exp(-U(x)/kT) / Zquad, -ℓ, ℓ; maxevals=meq)) +
                        collect(hquadrature(t -> -(1/t)*exp(-U(1/t)/kT) / t^2 / Zquad, 1/ℓ, 0.0; maxevals=meq))
                       )
                     );
(Uquad, Uquad_err) = (
                      collect(hquadrature(t -> -U(1/t)*exp(-U(1/t)/kT) / t^2 / Zquad, 0.0, -1/ℓ; maxevals=meq)) +
                      collect(hquadrature(x -> U(x)*exp(-U(x)/kT) / Zquad, -ℓ, ℓ; maxevals=meq)) +
                      collect(hquadrature(t -> -U(1/t)*exp(-U(1/t)/kT) / t^2 / Zquad, 1/ℓ, 0.0; maxevals=meq))
                     );
(x2quad, x2quad_err) = (
                      collect(hquadrature(t -> -(1/t)^2*exp(-U(1/t)/kT) / t^2 / Zquad, 0.0, -1/ℓ; maxevals=meq)) +
                      collect(hquadrature(x -> x^2*exp(-U(x)/kT) / Zquad, -ℓ, ℓ; maxevals=meq)) +
                      collect(hquadrature(t -> -(1/t)^2*exp(-U(1/t)/kT) / t^2 / Zquad, 1/ℓ, 0.0; maxevals=meq))
                     );
(U2quad, U2quad_err) = (
                        collect(hquadrature(t -> -(U(1/t))^2*exp(-U(1/t)/kT) / t^2 / Zquad, 0.0, -1/ℓ; maxevals=meq)) +
                        collect(hquadrature(x -> (U(x))^2*exp(-U(x)/kT) / Zquad, -ℓ, ℓ; maxevals=meq)) +
                        collect(hquadrature(t -> -(U(1/t))^2*exp(-U(1/t)/kT) / t^2 / Zquad, 1/ℓ, 0.0; maxevals=meq))
                     );

idx = rand(1:pargs["num-runs"]);
@info "std ar = $(results_std[idx][:ar])";
@info "polya ar = $(results_polya[idx][:ar])";
@info "umb std ar = $(results_umb_std[idx][:ar])";
@info "umb polya ar = $(results_umb_polya[idx][:ar])";
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
@info "<U^2> via quadrature = $U2quad; error = $U2quad_err";
@info "<U^2> via std = $(results_std[idx][:U2avg])";
@info "<U^2> via polya = $(results_polya[idx][:U2avg])";
@info "<U^2> via umb std = $(results_umb_std[idx][:U2avg])";
@info "<U^2> via umb polya = $(results_umb_polya[idx][:U2avg])";
#@info "<Urolling> via std = $(results_std[idx][:Urolling])";
#@info "<Urolling> via polya = $(results_polya[idx][:Urolling])";
#@info "<Urolling> via umb std = $(results_umb_std[idx][:Urolling])";
#@info "<Urolling> via umb polya = $(results_umb_polya[idx][:Urolling])";

L1_error_std, L1_error_polya = post_process_main_runs(results_std, 
                                                      results_polya, 
                                                      pargs; 
                                                      xrolling=xquad, 
                                                      Urolling=Uquad, 
                                                      x2rolling=x2quad,
                                                      U2rolling=U2quad);
L1_err_umb_std, L1_err_umb_polya = post_process_main_runs(results_umb_std, 
                                                          results_umb_polya,
                                                          pargs_umb; 
                                                          xrolling=xquad, 
                                                          Urolling=Uquad, 
                                                          x2rolling=x2quad,
                                                          U2rolling=U2quad);

if pargs["do-csvs"]
  for k in keys(L1_error_std)
    writedlm(joinpath(pargs["outdir"], "L1_$k.csv"), hcat(
                                                          results_std[1][:rolls],
                                                          L1_error_std[k],
                                                          L1_error_polya[k],
                                                          L1_err_umb_std[k],
                                                          L1_err_umb_polya[k]
                                                         ), ',');
  end
end

if pargs["do-conv-rates"]
  println("αs via std = $(convergence_rates(L1_error_std, results_std[1][:rolls]))");
  println("αs via polya = $(convergence_rates(L1_error_polya, results_std[1][:rolls]))");
  println("αs via umb = $(convergence_rates(L1_err_umb_std, results_std[1][:rolls]))");
  println("αs via gu = $(convergence_rates(L1_err_umb_polya, results_std[1][:rolls]))");
end

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
