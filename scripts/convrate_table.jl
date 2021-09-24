using DelimitedFiles
using GLM
using Glob
using DataFrames

col_idxs = ["rolls", "standard", "polya", "umbrella", "polya + umb."];

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
    for datafile in readdir(glob"L1_*.csv")
      datavar = split(split(datafile, '.')[1], '_')[2];
      println("---------------- $datavar ----------------------");
      raw_data = readdlm(datafile, ',');
      data = DataFrame(X = map(log, raw_data[:, 1]),
                       STD = map(log, raw_data[:, 2]),
                       POL = map(log, raw_data[:, 3]),
                       UMB = map(log, raw_data[:, 4]),
                       PUB = map(log, raw_data[:, 5]));

      @show std_fit = lm(@formula(STD ~ X), data);
      @show coef(std_fit);

      @show pol_fit = lm(@formula(POL ~ X), data);
      @show coef(pol_fit);

      @show umb_fit = lm(@formula(UMB ~ X), data);
      @show coef(umb_fit);

      @show pub_fit = lm(@formula(PUB ~ X), data);
      @show coef(pub_fit);
      println("---------------- $datavar ----------------------");
    end # L1 loop
    println("*****************************************");
    println("-----------------------------------------");
    cd("../");
  end # d2 loop
  cd("../");
  println("=========================================");
  println("*****************************************");
  println("=========================================");
end # d1 loop
