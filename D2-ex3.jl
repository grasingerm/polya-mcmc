using ArgParse;
using Distributions;
using LinearAlgebra;
using Logging;
using DelimitedFiles;
using ProfileView;
using DecFP;
using Quadmath;
using Plots;

s = ArgParseSettings();
@add_arg_table! s begin
  "--kT", "-T"
    help = "dimensionless temperature"
    arg_type = Float64
    default = 1.0
  "--num-steps", "-N"
    help = "number of steps"
    arg_type = Int
    default = convert(Int, 1e5)
  "--x0", "-X"
    help = "initial configuration"
    arg_type = String
    default = "[0.0; 0.0]"
  "--C1", "-k"
    help = "spring stiffness"
    arg_type = Float64
    default = 1.0
  "--C2", "-l"
    help = "aspect ratio"
    arg_type = Float64
    default = 4.0
  "--force", "-f"
    help = "applied force"
    arg_type = String
    default = "[0.0; 0.0]"
  "--orbit"
    help = "exploit underlying discrete symmetry by taking orbits during trial moves"
    action = :store_true
  "--dx", "-x"
    help = "maximum step length"
    arg_type = Float64;
    default = 0.25;
  "--step-adjust-lb", "-L"
    help = "adjust step sizes if acc. ratio below this threshold"
    arg_type = Float64
    default = 0.15
  "--step-adjust-ub", "-U"
    help = "adjust step sizes if acc. ratio above this threshold"
    arg_type = Float64
    default = 0.55
  "--step-adjust-scale", "-A"
    help = "scale factor for adjusting step sizes (> 1.0)"
    arg_type = Float64
    default = 1.1
  "--steps-per-adjust", "-S"
    help = "steps between step size adjustments"
    arg_type = Int
    default = 2500
#  "--acc", "-a"
#    help = "acceptance function (metropolis|kawasaki)"
#    arg_type = String
#    default = "metropolis"
  "--update-freq"
    help = "update frequency (seconds)"
    arg_type = Float64;
    default = 15.0;
  "--verbose", "-v"
    help = "verbosity level: 0-nothing, 1-errors, 2-warnings, 3-info"
    arg_type = Int
    default = 3
  "--stepout", "-s"
    help = "steps between storing rolling averages"
    arg_type = Int
    default = 500
  "--numeric-type"
    help = "numerical data type for averaging (float64|float128|dec128|big)"
    arg_type = String
    default = "float64"
  "--profile", "-Z"
    help = "profile the program"
    action = :store_true
end

pargs = parse_args(s);

if pargs["verbose"] == 3
  global_logger(ConsoleLogger(stderr, Logging.Info));
elseif pargs["verbose"] == 2
  global_logger(ConsoleLogger(stderr, Logging.Warn));
elseif pargs["verbose"] == 1
  global_logger(ConsoleLogger(stderr, Logging.Error));
else
  global_logger(Logging.NullLogger());
end

function mcmc(nsteps::Int, pargs)
  basis_vectors = [1 0; 0 1];
  kT = pargs["kT"];
  k = pargs["C1"];
  ℓ = pargs["C2"];
  f = eval(Meta.parse(pargs["force"]));
  U = (x) -> k*( (x[1]-ℓ/2)^2 + (x[1]+ℓ/2)^2 + 
                 (x[2]-1/2)^2 + (x[2]+1/2)^2   ) - dot(f, x);
  x = eval(Meta.parse(pargs["x0"]));
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
  xtotal = x;
  xrolling = Vector{Float64}[];
  Ucurr = U(x);
  Utotal = Ucurr;
  Urolling = Float64[];
  x2total = x .* x;
  xstd_rolling = Vector{Float64}[];
  U2total = Ucurr*Ucurr;
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
    if (Utrial < Ucurr) || (rand() <= exp(-(Utrial - Ucurr) / kT) )
      x = xtrial;
      Ucurr = Utrial;
      nacc += 1;
    end

    xtotal += x;
    Utotal += Ucurr;
    x2total += x .* x;
    U2total += Ucurr*Ucurr;
    ar = nacc / s;

    if s % stepout == 0
      push!(rolls, s);
      push!(xrolling, xtotal / s);
      push!(Urolling, Utotal / s);
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

pargs["orbit"] = false;
@show result_std = mcmc(pargs["num-steps"], pargs);

pargs["orbit"] = true;
@show result_polya = mcmc(pargs["num-steps"], pargs);

@show xrolling_std = transpose(hcat(result_std[:xrolling]...));
@show size(xrolling_std);
p = plot(result_std[:rolls], xrolling_std[:, 1]; label="std mcmc",
         xlabel = "step", ylabel = "\$\\langle x_1 \\rangle\$");
@show xrolling_polya = transpose(hcat(result_polya[:xrolling]...));
@show size(xrolling_polya);
plot!(result_polya[:rolls], xrolling_polya[:, 1]; label="polya");
display(p);
println();
println("Press RETURN to exit...");
readline();

p = plot(result_std[:rolls], xrolling_std[:, 2]; label="std mcmc",
         xlabel = "step", ylabel = "\$\\langle x_2 \\rangle\$");
plot!(result_polya[:rolls], xrolling_polya[:, 2]; label="polya");
display(p);
println();
println("Press RETURN to exit...");
readline();

p = plot(result_std[:rolls], result_std[:Urolling]; label="std mcmc",
         xlabel = "step", ylabel = "\$\\langle U \\rangle\$");
plot!(result_polya[:rolls], result_polya[:Urolling]; label="polya");
display(p);
println();
println("Press RETURN to exit...");
readline();
