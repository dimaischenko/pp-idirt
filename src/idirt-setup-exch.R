## project setup
# list with options for for each project
# we set pathes to data files, short experiment names
# rev -- is logical value for each experiment
#  "F" means : H / (H + L), "T" : means L / (H + L) values in experiment

l.prj <- list(
  IDIRT = list(
    name = "IDIRT",
    path = "data/part2/s_idirt.csv",
    cols = c("Ratio H/L Normalized 1c"),
    rev = c(F)),
  EX_velo = list(
    name = "ex_velo",
    path = "data/part2/s_mt302_velo.csv",
    cols = c("Ratio H/L t0",
             "Ratio H/L 3s",
             "Ratio H/L t5",
             "Ratio H/L 3m"),
    rev = c(F, F, F, F)),
  EX_qe = list(
    name = "ex_qe",
    path = "data/part2/s_mt302_qe.csv",
    cols = c("Ratio H/L t0",
             "Ratio H/L 3s",
             "Ratio H/L t5",
             "Ratio H/L 3m"),
    rev = c(F, F, F, F))
)
