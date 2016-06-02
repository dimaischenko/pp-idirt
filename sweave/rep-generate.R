# setwd
hdir <- "/home/dima/mega/pp-idirt/sweave/"
setwd(hdir)

# load data
load("../rdat/total_data.rda")

for (v.name in names(total.prj)) {
  setwd(hdir)
  dir.create(v.name)
  for (v.exp in names(total.prj[[v.name]][["cols"]])) {
    setwd(hdir)
    wdir <- sprintf("%s/%s", v.name, v.exp)
    #save(list = c("v.name", "v.exp"), file = "exp.rda")
    dir.create(wdir)
    setwd(wdir)
    Sweave("/home/dima/mega/pp-idirt/sweave/pp-report.Rnw")
  }
}
