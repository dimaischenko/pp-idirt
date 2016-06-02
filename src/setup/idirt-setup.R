## project setup
# list with options for for each project
# we set pathes to data files, short experiment names
# rev -- is logical value for each experiment
#  "F" means : H / (H + L), "T" : means L / (H + L) values in experiment

i.prj <- list(
  IDIRT = list(
    name = "IDIRT",
    path = "data/idirt/s_idirt.csv",
    cols = c("1a" = "1a",
             "1b" = "1b",
             "1c" = "1c",
             "2a" = "2a",
             "2b" = "2b",
             "2c" = "2c"),
    norm = T,
    rev = c(F, F, F, F, F, F)),

  RNAse = list(
    name = "RNAse",
    path = "data/idirt/s_rnas.csv",
    cols = c("3b" = "3b",
             "4a" = "4a"),
    norm = F,
    rev = c(F, T)),

  Tandem = list(
    name = "Tandem",
    path = "data/idirt/s_tand.csv",
    cols = c("e1" = "1f",
             "e2" = "1m"),
    norm = F,
    rev = c(F, T)),

  ORF2 = list(
    name = "ORF2",
    path = "data/idirt/s_561.csv",
    cols = c("fix3" = "KM040915_MAP_401_561_fix_3",
             "fix4" = "KM040915_MAP_401_561_fix_4",
             "nofix3" = "KM041015_MAP_401_561_nofix_3",
             "nofix4" = "KM041015_MAP_401_561_nofix_4"),
    norm = F,
    rev = c(F, T, T, F)),

  EndoMut = list(
    name = "EndoMut",
    path = "data/idirt/s_567.csv",
    cols = c("51f" = "51f",
             "52f" = "52f",
             "53f" = "53f"),
    norm = F,
    rev = c(F, F, F)),
  RevMut = list(
    name = "RevMut",
    path = "data/idirt/s_567.csv",
    cols = c("61f" = "61f",
             "62f" = "62f",
             "63f" = "63f"),
    norm = F,
    rev = c(F, F, F))
)
