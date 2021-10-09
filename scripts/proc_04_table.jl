using DelimitedFiles
using GLM
using Glob
using DataFrames

col_idxs = ["rolls", "standard", "polya", "umbrella", "polya + umb."];

function dir_to_params(dir)
  chunks = map(piece -> split(piece, "-"), split(dir, "_"));
  param_map = Dict();
  for chunk in chunks
    param_map[chunk[1]] = Meta.parse(chunk[2]) / 1000.0;
  end
  return param_map;
end

normalizers = Dict(
                   "refl-ex1" => Dict(
                                      "xrolling" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return sqrt(pm["b"] / (2*pm["a"]));
                                      end,
                                      "Urolling" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return (pm["b"])^2 / (4*pm["a"]);
                                      end,
                                      "x2rolling" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return pm["b"] / (2*pm["a"]);
                                      end,
                                      "U2rolling" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return ((pm["b"])^2 / (4*pm["a"]))^2;
                                      end
                                     ),
                   "trans-ex2" => Dict(
                                      "xrolling" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return (2*π) / pm["n"];
                                      end,
                                      "Urolling" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return pm["a"];
                                      end,
                                      "x2rolling" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return ((2*π) / pm["n"])^2;
                                      end,
                                      "U2rolling" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return (pm["a"])^2;
                                      end
                                      ),
                  );

for d1 in readdir()
  if isdir(d1)
    cd(d1)
  else
    continue
  end
  println("Processing $d1...");
  println("=========================================");
  println("*****************************************");
  println("=========================================");
  for d2 in readdir()
    if isdir(d2)
      cd(d2)
    else
      continue
    end
    println("Processing $d2...");
    println("-----------------------------------------");
    println("*****************************************");
    @show dir_to_params(d2);
    std_αs = [];
    umb_αs = [];
    pol_αs = [];
    pub_αs = [];
    std_ϵs = [];
    umb_ϵs = [];
    pol_ϵs = [];
    pub_ϵs = [];
    for datafile in readdir(glob"L1_*.csv")
      datavar = split(split(datafile, '.')[1], '_')[2];
      #println("---------------- $datavar ----------------------");
      raw_data = readdlm(datafile, ',');
      data = DataFrame(X = map(log, raw_data[:, 1]),
                       STD = map(log, raw_data[:, 2]),
                       POL = map(log, raw_data[:, 3]),
                       UMB = map(log, raw_data[:, 4]),
                       PUB = map(log, raw_data[:, 5]));

      std_fit = lm(@formula(STD ~ X), data);
      push!(std_αs, coef(std_fit)[2]);

      umb_fit = lm(@formula(UMB ~ X), data);
      push!(umb_αs, coef(umb_fit)[2]);
      
      pol_fit = lm(@formula(POL ~ X), data);
      push!(pol_αs, coef(pol_fit)[2]);

      pub_fit = lm(@formula(PUB ~ X), data);
      push!(pub_αs, coef(pub_fit)[2]);

      nm = normalizers[d1][datavar](d2);
      push!(std_ϵs, minimum(raw_data[:, 2]) / nm);
      push!(umb_ϵs, minimum(raw_data[:, 4]) / nm);
      push!(pol_ϵs, minimum(raw_data[:, 3]) / nm);
      push!(pub_ϵs, minimum(raw_data[:, 5]) / nm);

      #println("---------------- $datavar ----------------------");
    end # L1 loop
    println("$(minimum(std_αs)) & $(maximum(std_αs)) & $(minimum(std_ϵs)) & $(maximum(std_ϵs)) \\"); 
    println("$(minimum(umb_αs)) & $(maximum(umb_αs)) & $(minimum(umb_ϵs)) & $(maximum(umb_ϵs)) \\"); 
    println("$(minimum(pol_αs)) & $(maximum(pol_αs)) & $(minimum(pol_ϵs)) & $(maximum(pol_ϵs)) \\"); 
    println("$(minimum(pub_αs)) & $(maximum(pub_αs)) & $(minimum(pub_ϵs)) & $(maximum(pub_ϵs)) \\"); 
    println("*****************************************");
    println("-----------------------------------------");
    cd("../");
  end # d2 loop
  cd("../");
  println("=========================================");
  println("*****************************************");
  println("=========================================");
end # d1 loop
