include(joinpath(@__DIR__, "src", "bp.jl"));

s = ArgParseSettings();
@add_arg_table! s begin
  "--C1", "-a"
    help = "asymmetry coeff"
    arg_type = Float64
    default = 1.0
  "--C2", "-b"
    help = "valley coeff"
    arg_type = Float64
    default = 100.0
  "--x0", "-X"
    help = "initial configuration"
    arg_type = String
    default = "(::Any) -> [0.0; 0.0]"
  "--force", "-f"
    help = "applied force"
    arg_type = String
    default = "[0.0; 0.0]"
end

default_options = src_include("default_options.jl");
ArgParse.import_settings!(s, default_options);
pargs = src_include("parse_args.jl");
pargs["force"] = eval(Meta.parse(pargs["force"]));

@everywhere function mcmc(nsteps::Int, pargs)
  basis_vectors = [1 0; 0 1];
  kT = pargs["kT"];
  a = pargs["C1"];
  b = pargs["C2"];
  f = pargs["force"];
  U = (x) -> (a - x[1])^2 + b*(x[2] - x[1]^2)^2;
  x = pargs["x0"](pargs);
  xstep = pargs["dx"];
  dx_dist = Uniform(-xstep, xstep);
  orbf! = if (pargs["orbit"])
    x -> begin;
      if rand(Bool)
        x[1] = -x[1];
      end
    end
  else
    x -> begin; end;
  end

  nacc = 0;
  nacc_total = 0;
  natt = 0;
  xtotal = x;
  xrolling = Vector{Float64}[];
  Ucurr = U(x);
  Utotal = Ucurr;
  Urolling = Float64[];
  x2total = x .* x;
  x2rolling = Vector{Float64}[];
  xstd_rolling = Vector{Float64}[];
  U2total = Ucurr*Ucurr;
  U2rolling = Float64[];
  Ustd_rolling = Float64[];
  stepout = pargs["stepout"];
  rolls = Int[];

  start = time();
  last_update = start;
  for s = 1:nsteps
    xtrial = x + rand(dx_dist)*view(basis_vectors, :, rand(1:2));
    orbf!(xtrial);
    Utrial = U(xtrial);
    if (Utrial < Ucurr) || (rand() <= exp(-(Utrial - Ucurr) / kT) )
      x = xtrial;
      Ucurr = Utrial;
      nacc += 1;
      nacc_total += 1;
    end
    natt += 1;

    xtotal += x;
    Utotal += Ucurr;
    x2total += x .* x;
    U2total += Ucurr*Ucurr;
    ar = nacc_total / s;

    if s % stepout == 0
      push!(rolls, s);
      push!(xrolling, xtotal / s);
      push!(Urolling, Utotal / s);
      push!(x2rolling, x2total / s);
      push!(U2rolling, U2total / s);
      push!(
            xstd_rolling, 
            map(i -> sqrt(max(0.0, x2total[i] / s - (xtotal[i] / s)^2)), 1:2)
           );
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
       ar = nacc / natt;
       if (ar > pargs["step-adjust-ub"])
         @info "acceptance ratio is high; increasing step size";
         xstep *= pargs["step-adjust-scale"];
         nacc = 0;
         natt = 0;
       elseif ar < pargs["step-adjust-lb"]
         @info "acceptance ratio is low; decreasing step size";
         xstep /= pargs["step-adjust-scale"];
         nacc = 0;
         natt = 0;
       end
       dx_dist = Uniform(-xstep, xstep);
    end

  end

  return Dict(:xavg => xtotal / nsteps, :Uavg => Utotal / nsteps,
              :xrolling => xrolling, :Urolling => Urolling,
              :x2avg => x2total / nsteps, :U2avg => U2total / nsteps,
              :x2rolling => x2rolling, :U2rolling => U2rolling,
              :xstd_rolling => xstd_rolling, :Ustd_rolling => Ustd_rolling,
              :rolls => rolls, :ar => nacc_total / nsteps);

end

src_include("wrap_main_runs.jl");

results_std, results_polya = wrap_main_runs(pargs);
post_process_main_runs(results_std, results_polya, pargs)

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
