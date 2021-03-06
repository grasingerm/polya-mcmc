using GLM
using DataFrames

@everywhere import Base.+;

@everywhere function (+)(d1::Dict, d2::Dict)
  @assert(keys(d1)==keys(d2));
  ret = Dict(d1);
  for key in keys(d1)
    ret[key] += d2[key];
  end
  return ret;
end

@everywhere function wrap_main_runs(pargs::Dict)
  nsteps, nruns = pargs["num-steps"], pargs["num-runs"];
 
  elapsed = @elapsed begin;
    pargs_std = copy(pargs);
    pargs_std["outdir"] = joinpath(pargs["outdir"], "std");
    pargs_std["orbit"] = false;
    results_std = pmap((::Int) -> mcmc(nsteps, pargs_std), 1:nruns);

    pargs_polya = copy(pargs);
    pargs_polya["outdir"] = joinpath(pargs["outdir"], "polya");
    pargs_polya["orbit"] = true;
    results_polya = pmap((::Int) -> mcmc(nsteps, pargs_polya), 1:nruns);
  end

  @info "Simulation time: $elapsed";

  return results_std, results_polya;
end

function post_process_main_runs(results_std, results_polya, pargs; kwargs...)
  L1_error_std = Dict();
  L1_error_polya = Dict();
  for (k, v) in kwargs
    L1_error_std[k] = sum(map(i -> begin;
                                map(x -> norm(x - v, 2), results_std[i][k])
                      end, 1:length(results_std))) / length(results_std);
    L1_error_polya[k] = sum(map(i -> begin;
                          map(x -> norm(x - v, 2), results_polya[i][k])
                        end, 1:length(results_polya))) / length(results_polya);
    if length(v) > 1
      for m=1:length(v)
        L1_error_std[String(k)*"-$m"] = sum(map(i -> begin;
                                 map(x -> abs(x[m] - v[m]), results_std[i][k])
                          end, 1:length(results_std))) / length(results_std);
        L1_error_polya[String(k)*"-$m"] = sum(map(i -> begin;
                                 map(x -> abs(x[m] - v[m]), results_polya[i][k])
                            end, 1:length(results_polya))) / length(results_polya);
      end
    end
  end
  outdir = pargs["outdir"];
  mkpath(outdir);
  if pargs["do-csvs"]
    for k in keys(results_std[1])
      writedlm(
               joinpath(outdir, "$(k)_std.csv"),
               hcat(map(i -> results_std[i][k], 1:length(results_std))...),
               ','
              );
      writedlm(
               joinpath(outdir, "$(k)_polya.csv"),
               hcat(map(i -> results_polya[i][k], 1:length(results_polya))...),
               ','
              );
    end
  end
  return L1_error_std, L1_error_polya;
end

function convergence_rates(error_table, rolls; stars=Dict())
  ??s = Dict();
  log_rolls = map(log, rolls);
  for (k, v) in error_table
    log_y = map(log, v);
    d = (haskey(stars, k)) ? stars[k] : 1.0;
    data = DataFrame(X = log_rolls, Y = log_y / d);
    fit = lm(@formula(Y ~ X), data);
    ??s[k] = coef(fit)[2];
  end
  return ??s;
end
