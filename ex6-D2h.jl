include(joinpath(@__DIR__, "src", "bp.jl"));

using Cubature;

s = ArgParseSettings();
@add_arg_table! s begin
  "--len-x", "-a"
    help = "length in the x-direction"
    arg_type = Float64
    default = 1.0
  "--len-y", "-b"
    help = "length in the y-direction"
    arg_type = Float64
    default = 0.5
  "--len-z", "-c"
    help = "length in the z-direction"
    arg_type = Float64
    default = 2.0
  "--charge", "-q"
    help = "charge of ion"
    arg_type = Float64
    default = 1.0
  "--x0", "-X"
    help = "initial configuration"
    arg_type = String
    default = "(args::Any) -> [rand(Uniform(-args[\"len-x\"], args[\"len-x\"])), rand(Uniform(-args[\"len-y\"], args[\"len-y\"])), rand(Uniform(-args[\"len-z\"], args[\"len-z\"]))]"
  "--umbrella-sigma-a"
    help = "umbrella weight width"
    arg_type = Float64
    default = 1.0
  "--umbrella-sigma-b"
    help = "umbrella weight width"
    arg_type = Float64
    default = 1.0
  "--umbrella-sigma-c"
    help = "umbrella weight width"
    arg_type = Float64
    default = 1.0
  "--umbrella-C"
    help = "umbrella peak height"
    arg_type = Float64
    default = 1.0
  "--force", "-f"
    help = "applied force"
    arg_type = String
    default = "[0.0; 0.0; 0.0]"
end

default_options = src_include("default_options.jl");
ArgParse.import_settings!(s, default_options);
pargs = src_include("parse_args.jl");
pargs["force"] = eval(Meta.parse(pargs["force"]));

@everywhere function pbc!(xs::Vector{Float64}, ℓs::Vector{Float64})
  xs[:] = map(pair -> begin
         x, ℓ = pair;
         return if x <= -3*ℓ/2
           3*ℓ + x;
         elseif x >= 3*ℓ/2
           x - 3*ℓ;
         else
           x
         end
       end, zip(xs, ℓs))
end

@everywhere function mcmc(nsteps::Int, pargs)
  basis_vectors = [1 0 0; 0 1 0; 0 0 1];
  kT = pargs["kT"];
  a = pargs["len-x"];
  b = pargs["len-y"];
  c = pargs["len-z"];
  f = pargs["force"];
  q = pargs["charge"];
  wt = pargs["weight"];
  latvecs = [
             [0; 0; 0],
             [-a; 0; 0],
             [a; 0; 0],
             [0; -b; 0],
             [0; b; 0],
             [0; 0; -c],
             [0; 0; c]
            ];
  U = (x) -> q*sum(map(v -> 1/norm(x - v, 2), latvecs)) - dot(f, x)
  x = pargs["x0"](pargs);
  xstep = pargs["dx"];
  dx_dist = Uniform(-xstep, xstep);
  orbf! = if (pargs["orbit"])
    x -> begin;
      choice = rand(1:4);
      if choice == 1
        x[1] = -x[1];
      elseif choice == 2
        x[2] = -x[2];
      elseif choice == 3
        x[3] = -x[3];
      end
    end
  else
    x -> x;
  end

  wt_curr = wt(x);
  inv_wt_total = 1 / wt_curr;
  xtotal = x[:] / wt_curr;
  xrolling = Vector{Float64}[];
  Ucurr = U(x);
  Utotal = Ucurr / wt_curr;
  Urolling = Float64[];
  x2total = (x[:] .* x[:]) / wt_curr;
  x2rolling = Vector{Float64}[];
  xstd_rolling = Vector{Float64}[];
  U2total = Ucurr*Ucurr / wt_curr;
  U2rolling = Float64[];
  Ustd_rolling = Float64[];
  stepout = pargs["stepout"];
  rolls = Int[];

  natt = 0;
  nacc = 0;
  nacc_total = 0;

  start = time();
  last_update = start;
  for s = 1:nsteps
    xtrial = x[:] + rand(dx_dist, 3);
    orbf!(xtrial);
    Utrial = U(xtrial);
    wt_trial = wt(xtrial);
    if (
        (abs(xtrial[1]) <= a) &&
        (abs(xtrial[2]) <= b) &&
        (abs(xtrial[3]) <= c) &&
        (rand() <= ( exp(-(Utrial - Ucurr) / kT) * wt_trial / wt_curr ))
       )
      x[:] = xtrial[:];
      Ucurr = Utrial;
      wt_curr = wt_trial;
      nacc += 1;
      nacc_total += 1;
    end
    natt += 1;
    @assert((abs(x[1]) <= a) &&
            (abs(x[2]) <= b) &&
            (abs(x[3]) <= c));

    inv_wt_total += 1 / wt_curr;
    xtotal[:] += x[:] / wt_curr;
    Utotal += Ucurr / wt_curr;
    x2total[:] += (x[:] .* x[:]) / wt_curr;
    U2total += Ucurr*Ucurr / wt_curr;
    ar = nacc_total / s; # rolling acceptance ratio

    if s % stepout == 0
      push!(rolls, s);
      push!(xrolling, xtotal[:] / inv_wt_total);
      push!(Urolling, Utotal / inv_wt_total);
      push!(x2rolling, x2total[:] / inv_wt_total);
      push!(U2rolling, U2total / inv_wt_total);
      push!(
            xstd_rolling, 
            map(i -> sqrt(max(0.0, x2total[i] / inv_wt_total 
                              - (xtotal[i] / inv_wt_total)^2)), 1:3)
           );
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
       ar = nacc / natt;
       if (ar > pargs["step-adjust-ub"])
         @info "acceptance ratio is high; increasing step size", xstep;
         nacc = 0;
         natt = 0;
         xstep *= pargs["step-adjust-scale"];
       elseif ar < pargs["step-adjust-lb"]
         @info "acceptance ratio is low; decreasing step size", xstep;
         nacc = 0;
         natt = 0;
         xstep /= pargs["step-adjust-scale"];
       end
       dx_dist = Uniform(-xstep, xstep);
    end

  end

  return Dict(:xavg => xtotal / inv_wt_total, :Uavg => Utotal / inv_wt_total,
              :xrolling => xrolling, :Urolling => Urolling,
              :x2avg => x2total / inv_wt_total, :U2avg => U2total / inv_wt_total,
              :x2rolling => x2rolling, :U2rolling => U2rolling,
              :xstd_rolling => xstd_rolling, :Ustd_rolling => Ustd_rolling,
              :rolls => rolls,
              :ar => nacc / natt,
              :AR => nacc_total / inv_wt_total
             );

end

src_include("wrap_main_runs.jl");

results_std, results_polya = wrap_main_runs(pargs);

uC = pargs["umbrella-C"];
uΣ = 1/2*[1/(pargs["umbrella-sigma-a"]^2) 0 0;
          0 1/(pargs["umbrella-sigma-b"]^2) 0;  
          0 0 1/(pargs["umbrella-sigma-c"]^2)];  
pargs_umb = copy(pargs);
pargs_umb["weight"] = x -> uC*exp(-dot(x, uΣ*x));
pargs_umb["outdir"] = joinpath(pargs["outdir"], "umb");
results_umb_std, results_umb_polya = wrap_main_runs(pargs_umb);

kT = pargs["kT"];
a = pargs["len-x"];
b = pargs["len-y"];
c = pargs["len-z"];
q = pargs["charge"];
f = pargs["force"];
latvecs = [
             [0; 0; 0],
             [-a; 0; 0],
             [a; 0; 0],
             [0; -b; 0],
             [0; b; 0],
             [0; 0; -c],
             [0; 0; c]
            ];
U = (x) -> q*sum(map(v -> 1/norm(x - v, 2), latvecs)) - dot(f, x)
meq = 25000;
(Zquad, Zquad_err) = hcubature(x -> exp(-U(x) / kT), [-a, -b, -c], [a, b, c]; maxevals=meq);
(xquad, xquad_err) = hcubature(3, (x, v) -> (v[:] = x*exp(-U(x) / kT) / Zquad), [-a, -b, -c], [a, b, c]; maxevals=meq);
(Uquad, Uquad_err) = (
                      collect(hcubature(x -> U(x)*exp(-U(x) / kT) / Zquad, [0, -b, -c], [a, b, c]; maxevals=meq)) +
                      collect(hcubature(x -> U(x)*exp(-U(x) / kT) / Zquad, [-a, -b, -c], [0, b, c]; maxevals=meq))
                     );
(x2quad, x2quad_err) = hcubature(3, (x, v) -> (v[:] = x .* x *exp(-U(x) / kT) / Zquad), [-a, -b, -c], [a, b, c]; maxevals=meq);
(U2quad, U2quad_err) = (
                        collect(hcubature(x -> (U(x)^2)*exp(-U(x) / kT) / Zquad, [0, -b, -c], [a, b, c]; maxevals=meq)) +
                        collect(hcubature(x -> (U(x)^2)*exp(-U(x) / kT) / Zquad, [-a, -b, -c], [0, b, c]; maxevals=meq))
                     );

idx = rand(1:pargs["num-runs"]);
@info "std ar = $(results_std[idx][:ar])";
@info "polya ar = $(results_polya[idx][:ar])";
@info "umb std ar = $(results_umb_std[idx][:ar])";
@info "umb polya ar = $(results_umb_polya[idx][:ar])";
@info "std AR = $(results_std[idx][:AR])";
@info "polya AR = $(results_polya[idx][:AR])";
@info "umb std AR = $(results_umb_std[idx][:AR])";
@info "umb polya AR = $(results_umb_polya[idx][:AR])";
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
  xrolling_std = transpose(hcat(results_std[idx][:xrolling]...));
  p = plot(results_std[idx][:rolls], xrolling_std[:, 1]; label="std mcmc",
           xlabel = "step", ylabel = "\$\\langle x_1 \\rangle\$");
  xrolling_polya = transpose(hcat(results_polya[idx][:xrolling]...));
  plot!(results_polya[idx][:rolls], xrolling_polya[:, 1]; label="polya");
  display(p);
  println();
  println("Press RETURN to exit...");
  readline();

  p = plot(results_std[idx][:rolls], xrolling_std[:, 2]; label="std mcmc",
           xlabel = "step", ylabel = "\$\\langle x_2 \\rangle\$");
  plot!(results_polya[idx][:rolls], xrolling_polya[:, 2]; label="polya");
  display(p);
  println();
  println("Press RETURN to exit...");
  readline();

  p = plot(results_std[idx][:rolls], results_std[idx][:Urolling]; 
           label="std mcmc", xlabel = "step", ylabel = "\$\\langle U \\rangle\$");
  plot!(results_polya[idx][:rolls], results_polya[idx][:Urolling]; label="polya");
  display(p);
  println();
  println("Press RETURN to exit...");
  readline();
end
