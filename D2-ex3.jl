include(joinpath(@__DIR__, "src", "bp.jl"));

s = ArgParseSettings();
@add_arg_table! s begin
  "--C1", "-k"
    help = "spring stiffness"
    arg_type = Float64
    default = 1.0
  "--C2", "-l"
    help = "aspect ratio"
    arg_type = Float64
    default = 4.0
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
  k = pargs["C1"];
  ℓ = pargs["C2"];
  f = pargs["force"];
  wt = pargs["weight"];
  U = (x) -> k*( (x[1]-ℓ/2)^2 + (x[1]+ℓ/2)^2 + 
                 (x[2]-1/2)^2 + (x[2]+1/2)^2   ) - dot(f, x);
  x = pargs["x0"](pargs);
  xstep = pargs["dx"];
  dx_dist = Uniform(-xstep, xstep);
  orbf! = if (pargs["orbit"])
    x -> begin;
      choice = rand(1:4);
      # choice == 4 is the identity
      if choice == 1
        x[1] = -x[1];
      elseif choice == 2
        x[2] = -x[2];
      elseif choice == 3
        x[:] = -x;
      end
    end
  else
    x -> begin; end;
  end

  nacc = 0;
  wt_curr = wt(x);
  inv_wt_total = 1 / wt_curr;
  xtotal = x * inv_wt_total;
  xrolling = Vector{Float64}[];
  Ucurr = U(x);
  Utotal = Ucurr * inv_wt_total;
  Urolling = Float64[];
  x2total = (x .* x) * inv_wt_total;
  x2rolling = Vector{Float64}[];
  xstd_rolling = Vector{Float64}[];
  U2total = Ucurr*Ucurr * inv_wt_total;
  U2rolling = Float64[];
  Ustd_rolling = Float64[];
  stepout = pargs["stepout"];
  rolls = Int[];

  start = time();
  last_update = start;
  for s = 1:nsteps
    xtrial = x + rand(dx_dist)*view(basis_vectors, :, rand(1:2));
    orbf!(xtrial);
    xtrial[1] = max(-ℓ/2, min(ℓ/2, xtrial[1]));
    xtrial[2] = max(-1/2, min(1/2, xtrial[2]));
    Utrial = U(xtrial);
    wt_trial = wt(xtrial);
    if rand() <= ( exp(-(Utrial - Ucurr) / kT) * (wt_trial / wt_curr) )
      x = xtrial;
      Ucurr = Utrial;
      wt_curr = wt_trial;
      nacc += 1;
    end

    xtotal += x / wt_curr;
    Utotal += Ucurr / wt_curr;
    x2total += (x .* x) / wt_curr;
    U2total += (Ucurr*Ucurr) / wt_curr;
    inv_wt_total += 1 / wt(x);
    ar = nacc / s;

    if s % stepout == 0
      push!(rolls, s);
      inv_wt = inv_wt_total;
      push!(xrolling, xtotal / inv_wt);
      push!(Urolling, Utotal / inv_wt);
      push!(x2rolling, x2total / inv_wt);
      push!(U2rolling, U2total / inv_wt);
      push!(
            xstd_rolling, 
            map(i -> sqrt(max(0.0, x2total[i] / inv_wt - 
                                   (xtotal[i] / inv_wt)^2)), 1:2)
           );
      push!(Ustd_rolling, sqrt(max(0.0, U2total / inv_wt - 
                                        (Utotal / inv_wt)^2)));
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

  inv_wt = inv_wt_total;
  return Dict(:xavg => xtotal / inv_wt, 
              :Uavg => Utotal / inv_wt,
              :xrolling => xrolling, :Urolling => Urolling,
              :x2avg => x2total / inv_wt, 
              :U2avg => U2total / inv_wt,
              :x2rolling => x2rolling, :U2rolling => U2rolling,
              :xstd_rolling => xstd_rolling, 
              :Ustd_rolling => Ustd_rolling,
              :rolls => rolls, :ar => nacc / nsteps);

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
