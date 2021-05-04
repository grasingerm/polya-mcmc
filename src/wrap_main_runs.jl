@everywhere import Base.+;

@everywhere function (+)(d1::Dict, d2::Dict)
  @assert(keys(d1)==keys(d2));
  ret = Dict(d1);
  for key in keys(d1)
    ret[key] += d2[key];
  end
  return ret;
end

pargs["orbit"] = false;
@show result_std = pmap(i -> mcmc(pargs["num-steps"], pargs),
                        1:pargs["num-runs"]);

pargs["orbit"] = true;
@show result_polya = pmap(i -> mcmc(pargs["num-steps"], pargs),
                          1:pargs["num-runs"]);

#=
for (k, v) in result_std
  result_std[k] = v / pargs["num-runs"];
end
for (k, v) in result_polya
  result_polya[k] = v / pargs["num-runs"];
end
=#

result_std, result_polya
