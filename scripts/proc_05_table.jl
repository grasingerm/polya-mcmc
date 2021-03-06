using DelimitedFiles
using GLM
using Glob
using DataFrames

col_idxs = ["rolls", "standard", "polya", "umbrella", "polya + umb."];

function dir_to_params(dir)
  @show chunks = map(piece -> split(piece, "-"), split(dir, "_"));
  param_map = Dict();
  for chunk in chunks
    param_map[chunk[1]] = Meta.parse(chunk[end]) / 1000.0;
  end
  return param_map;
end

normalizers = Dict(
                   "D2h-ex6" => Dict(
                                      "xrolling" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return sqrt(pm["a"]^2 +
                                                    pm["b"]^2 +
                                                    pm["c"]^2);
                                      end,
                                      "xrolling-1" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return pm["a"];
                                      end,
                                      "xrolling-2" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return pm["b"];
                                      end,
                                      "xrolling-3" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return pm["c"];
                                      end,
                                      "Urolling" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return (pm["q"]);
                                      end,
                                      "x2rolling" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return (pm["a"])^2 + (pm["b"])^2 + (pm["c"])^2;
                                      end,
                                      "x2rolling-1" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return (pm["a"])^2;
                                      end,
                                      "x2rolling-2" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return (pm["b"])^2;
                                      end,
                                      "x2rolling-3" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return (pm["c"])^2;
                                      end,
                                      "U2rolling" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return (pm["q"])^2;
                                      end
                                     ),
                   "D2h-ex7" => Dict(
                                      "xrolling" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return (2^(1/6))*pm["sigma"];
                                      end,
                                      "xrolling-1" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return (2^(1/6))*pm["sigma"];
                                      end,
                                      "xrolling-2" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return (2^(1/6))*pm["sigma"];
                                      end,
                                      "xrolling-3" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return (2^(1/6))*pm["sigma"];
                                      end,
                                      "Urolling" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return (pm["eps"]);
                                      end,
                                      "x2rolling" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return ((2^(1/6))*pm["sigma"])^2;
                                      end,
                                      "x2rolling-1" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return ((2^(1/6))*pm["sigma"])^2;
                                      end,
                                      "x2rolling-2" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return ((2^(1/6))*pm["sigma"])^2;
                                      end,
                                      "x2rolling-3" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return ((2^(1/6))*pm["sigma"])^2;
                                      end,
                                      "U2rolling" => (dir) -> begin;
                                        pm = dir_to_params(dir);
                                        return (pm["eps"])^2;
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
    std_??s = [];
    umb_??s = [];
    pol_??s = [];
    pub_??s = [];
    std_??s = [];
    umb_??s = [];
    pol_??s = [];
    pub_??s = [];
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
      push!(std_??s, coef(std_fit)[2]);

      umb_fit = lm(@formula(UMB ~ X), data);
      push!(umb_??s, coef(umb_fit)[2]);
      
      pol_fit = lm(@formula(POL ~ X), data);
      push!(pol_??s, coef(pol_fit)[2]);

      pub_fit = lm(@formula(PUB ~ X), data);
      push!(pub_??s, coef(pub_fit)[2]);

      nm = normalizers[d1][datavar](d2);
      push!(std_??s, minimum(raw_data[:, 2]) / nm);
      push!(umb_??s, minimum(raw_data[:, 4]) / nm);
      push!(pol_??s, minimum(raw_data[:, 3]) / nm);
      push!(pub_??s, minimum(raw_data[:, 5]) / nm);

      #println("---------------- $datavar ----------------------");
    end # L1 loop
    println("$(minimum(std_??s)) & $(maximum(std_??s)) & $(minimum(std_??s)) & $(maximum(std_??s)) \\"); 
    println("$(minimum(umb_??s)) & $(maximum(umb_??s)) & $(minimum(umb_??s)) & $(maximum(umb_??s)) \\"); 
    println("$(minimum(pol_??s)) & $(maximum(pol_??s)) & $(minimum(pol_??s)) & $(maximum(pol_??s)) \\"); 
    println("$(minimum(pub_??s)) & $(maximum(pub_??s)) & $(minimum(pub_??s)) & $(maximum(pub_??s)) \\"); 
    println("*****************************************");
    println("-----------------------------------------");
    cd("../");
  end # d2 loop
  cd("../");
  println("=========================================");
  println("*****************************************");
  println("=========================================");
end # d1 loop
