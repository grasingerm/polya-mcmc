default_settings = ArgParseSettings();
@add_arg_table! default_settings begin
  "--kT", "-T"
    help = "dimensionless temperature"
    arg_type = Float64
    default = 1.0
  "--num-steps", "-N"
    help = "number of steps"
    arg_type = Int
    default = convert(Int, 1e5)
  "--num-runs", "-R"
    help = "number of steps"
    arg_type = Int
    default = 10
  "--outdir", "-O"
    help = "output directory"
    arg_type = String
    default = "."
  "--do-plots", "-P"
    help = "create interesting plots"
    action = :store_true
  "--do-csvs", "-C"
    help = "create csv files of data"
    action = :store_true
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
