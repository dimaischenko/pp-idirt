## project setup
# list with options for for each project
# we set pathes to data files, short experiment names
# rev -- is logical value for each experiment
#  "F" means : H / (H + L), "T" : means L / (H + L) values in experiment

i.prj <- list(
  IDIRT = list(
    name = "IDIRT",
    path = "data/idirt/s_idirt.csv",
    cols = c("1a" = "Ratio H/L Normalized 1a",
             "1b" = "Ratio H/L Normalized 1b",
             "1c" = "Ratio H/L Normalized 1c",
             "2a" = "Ratio H/L Normalized 2a",
             "2b" = "Ratio H/L Normalized 2b",
             "2c" = "Ratio H/L Normalized 2c"),
    rev = c(F, F, F, F, F, F)),
  
  RNAse = list(
    name = "RNAse",
    path = "data/idirt/s_rnas.csv",
    cols = c("3b" = "Ratio H/L 3b",
             "4a" = "Ratio H/L 4a"),
    rev = c(F, T)),
  
  Tandem = list(
    name = "Tandem",
    path = "data/idirt/s_tand.csv",
    cols = c("e1" = "Ratio H/L 1f",
             "e2" = "Ratio H/L 1m"),
    rev = c(F, T)),
  
  ORF2 = list(
    name = "ORF2",
    path = "data/idirt/s_561.csv",
    cols = c("fix3" = "Ratio H/L KM040915_MAP_401_561_fix_3",
             "fix4" = "Ratio H/L KM040915_MAP_401_561_fix_4",
             "nofix3" = "Ratio H/L KM041015_MAP_401_561_nofix_3",
             "nofix4" = "Ratio H/L KM041015_MAP_401_561_nofix_4"),
    rev = c(F, T, T, F)),
  
  EndoMut = list(
    name = "EndoMut",
    path = "data/idirt/s_567.csv",
    cols = c("51f" = "Ratio H/L 51f",
             "52f" = "Ratio H/L 52f",
             "53f" = "Ratio H/L 53f"),
    rev = c(F, F, F)),
  RevMut = list(
    name = "RevMut",
    path = "data/idirt/s_567.csv",
    cols = c("61f" = "Ratio H/L 61f",
             "62f" = "Ratio H/L 62f",
             "63f" = "Ratio H/L 63f"),
    rev = c(F, F, F))
)
