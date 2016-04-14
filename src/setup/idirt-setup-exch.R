## project setup
# list with options for for each project
# we set pathes to data files, short experiment names
# rev -- is logical value for each experiment
#  "F" means : H / (H + L), "T" : means L / (H + L) values in experiment

i.prj <- list(
  IDIRT = list(
    name = "IDIRT",
    path = "data/idirt/s_idirt.csv",
    cols = c("1c" = "Ratio H/L Normalized 1c"),
    rev = c(F)),
  EX_velo = list(
    name = "ex_velo",
    path = "data/idirt/s_mt302_velo.csv",
    cols = c("v t0" = "Ratio H/L t0",
             "v 3s" = "Ratio H/L 3s",
             "v t5" = "Ratio H/L t5",
             "v 3m" = "Ratio H/L 3m"),
    rev = c(F, F, F, F)))
#
#  EX_qe = list(
#    name = "ex_qe",
#    path = "data/idirt/s_mt302_qe.csv",
#    cols = c("q t0" = "Ratio H/L t0",
#             "q 3s" = "Ratio H/L 3s",
#             "q t5" = "Ratio H/L t5",
#             "q 3m" = "Ratio H/L 3m"),
#    rev = c(F, F, F, F))
#)
