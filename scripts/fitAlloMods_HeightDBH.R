# fit allometric model ----

plotSets <- c("WYK", "NatureParks", "NSSF")

mods_HD <- lapply(plotSets, function(s) {
    with(treeHeights[[s]],
         modelHD(D = DBH,
                 H = height_avg,
                 method = "log1"))
})
names(mods_HD) <- plotSets