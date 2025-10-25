plotSets <- c("WYK", "NatureParks", "NSSF")
mod_HD <- lapply(plotSets, function(s) {
    with(treeHeights[[s]],
         modelHD(D = DBH,
                 H = height_avg,
                 method = "log1"))
})
names(mod_HD) <- plotSets

## trees in 20 x 20 m core plot area ----

trees$predHeights <- NA
for(i in 1:length(plotSets)) {
    trees$predHeights[trees$plotSet == plotSets[i]] <- retrieveH(D = trees$DBH[trees$plotSet == plotSets[i]],
              model = mod_HD[[plotSets[i]]])$H
}

trees$predHeights[is.na(trees$predHeights)] <- retrieveH(D = trees$DBH[is.na(trees$predHeights)],
          model = mod_HD[["WYK"]])$H

trees$genus <- sapply(trees$species, function(x) {
    strsplit(x, split = " ")[[1]][1]
})

trees$WD <- getWoodDensity(
    genus = trees$genus,
    species = trees$species,
    stand = trees$plotID
)$meanWD

trees$AGB <- computeAGB(D = trees$DBH,
                        WD = trees$WD,
                        H = trees$predHeights)

treeComm_AGB <- as.data.frame.matrix(xtabs(AGB ~ plotID + species,
                                           subset = !is.na(AGB),
                                           data = trees))

## big trees in 40 x 40 m extended plot area

bigTrees$predHeights <- NA
for(i in 1:length(plotSets)) {
    bigTrees$predHeights[bigTrees$plotSet == plotSets[i]] <- retrieveH(D = bigTrees$DBH[bigTrees$plotSet == plotSets[i]],
                                                                       model = mod_HD[[plotSets[i]]])$H
}

bigTrees$predHeights[is.na(bigTrees$predHeights)] <- retrieveH(D = bigTrees$DBH[is.na(bigTrees$predHeights)],
                                                                   model = mod_HD[["WYK"]])$H

bigTrees$genus <- sapply(bigTrees$species, function(x) {
    strsplit(x, split = " ")[[1]][1]
})

bigTrees$WD <- getWoodDensity(
    genus = bigTrees$genus,
    species = bigTrees$species,
    stand = bigTrees$plotID
)$meanWD

bigTrees$AGB <- computeAGB(D = bigTrees$DBH,
                        WD = bigTrees$WD,
                        H = bigTrees$predHeights)

bigTreeComm_AGB <- as.data.frame.matrix(xtabs(AGB ~ plotID + species,
                                           subset = !is.na(AGB),
                                           data = bigTrees))

plots_AGBsum <- data.frame(plotID = rownames(treeComm_AGB),
                           trees = rowSums(treeComm_AGB)/(20*20/10000)) %>%
    full_join(data.frame(plotID = rownames(bigTreeComm_AGB),
                         bigTrees = rowSums(bigTreeComm_AGB)/(40*40/10000)),
              by = "plotID") %>%
    left_join(select(trees, plotID, plotSet) %>% unique(), by = "plotID")