pargs = parse_args(s);
pargs["x0"] = eval(Meta.parse(pargs["x0"]));
passobj(myid(), workers(), :pargs);

@everywhere begin;
  if pargs["verbose"] == 3
    global_logger(ConsoleLogger(stderr, Logging.Info));
  elseif pargs["verbose"] == 2
    global_logger(ConsoleLogger(stderr, Logging.Warn));
  elseif pargs["verbose"] == 1
    global_logger(ConsoleLogger(stderr, Logging.Error));
  else
    global_logger(Logging.NullLogger());
  end
end

pargs
