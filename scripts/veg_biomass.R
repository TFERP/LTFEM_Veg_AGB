# Load libraries ----

library(tidyverse)
library(BIOMASS)


# Read data ----

sourceFilesDir <- "./data/classified/dataProducts/LTFEM_Veg_(trial)_20250417/LTFEM_Veg/"

trees <- read.csv(paste0(sourceFilesDir, "WYK_trees_2022.csv")) %>%
    filter(pos_cirsq != "Arc" & pos_cirsq != "Outside") %>%
    mutate(species = str_replace(species, "cf. ", ""),
           plotSet = "WYK") %>%
    select(plotSet, plotID, treeID, species, DBH, date, coord_x, coord_y) %>%
    filter(!is.na(DBH)) %>%
    rbind(
        read.csv(paste0(sourceFilesDir, "NatureParks_trees_2023.csv")) %>%
            mutate(species = str_replace(species, "cf. ", ""),
                   DBH = as.numeric(DBH),
                   plotSet = "NatureParks") %>%
            select(plotSet, plotID, treeID, species, DBH, date, coord_x, coord_y) %>%
            filter(!is.na(DBH))
    ) %>%
    rbind(
        read.csv(paste0(sourceFilesDir, "NSSF_trees_2020.csv")) %>%
            mutate(species = str_replace(species, "cf. ", ""),
                   DBH = as.numeric(DBH),
                   plotSet = "NSSF") %>%
            select(plotSet, plotID, treeID, species, DBH, date, coord_x, coord_y) %>%
            filter(!is.na(DBH))
    )

bigTrees <- trees %>%
    filter(DBH >= 30) %>%
    rename(bigTreeID = treeID) %>%
    rbind(read.csv(paste0(sourceFilesDir, "WYK_bigTrees_2022.csv")) %>%
              mutate(species = str_replace(species, "cf. ", ""),
                     plotSet = "WYK") %>%
              select(plotSet, plotID, bigTreeID, species, DBH, date, coord_x, coord_y)) %>%
    rbind(read.csv(paste0(sourceFilesDir, "NatureParks_bigTrees_2023.csv")) %>%
              mutate(species = str_replace(species, "cf. ", ""),
                     plotSet = "NatureParks") %>%
              select(plotSet, plotID, bigTreeID, species, DBH, date, coord_x, coord_y)) %>%
    rbind(read.csv(paste0(sourceFilesDir, "NSSF_bigTrees_2020.csv")) %>%
              mutate(species = str_replace(species, "cf. ", ""),
                     DBH = as.numeric(DBH),
                     plotSet = "NSSF") %>%
              select(plotSet, plotID, bigTreeID, species, DBH, date, coord_x, coord_y)) %>%
    filter(!is.na(DBH))

source("./scripts/fitAlloMods_HeightDBH.R")

## trees in 20 x 20 m core plot area ----

trees$predHeights <- NA
for(i in 1:length(plotSets)) {
    trees$predHeights[trees$plotSet == plotSets[i]] <- retrieveH(D = trees$DBH[trees$plotSet == plotSets[i]],
              model = mods_HD[[plotSets[i]]])$H
}


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
                                                                       model = mods_HD[[plotSets[i]]])$H
}

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

# library(ggplot2)
# ggplot(data = plots_AGBsum) +
#     geom_point(mapping = aes(y = trees, x = bigTrees, col = plotSet)) +
#     scale_y_continuous(trans = "log10") +
#     scale_x_continuous(trans = "log10") +
#     geom_abline(intercept = 0, slope = 1)

# Comparisons with CISG values ----

plots_AGBsum <- plots_AGBsum %>%
    left_join(rbind(read.csv("./data/classified/CCNR_ACD_byplot.csv", header = TRUE),
                    read.csv("./data/classified/TP_ACD_byplot.csv", header = TRUE)) %>%
                  select(-1) %>%
                  rename(plotID = Plot,
                         CISG = ACDSum),
              by = "plotID") %>%
    left_join(read.csv("./data/classified/NParks_AGB_And_Others.csv") %>%
                  select(plot_id, Mg_ha) %>%
                  rename(plotID = plot_id,
                         Ellie = Mg_ha),
              by = "plotID")



ggplot(data = plots_AGBsum %>%
           mutate(trees_ACD = trees*0.45) %>%
           filter(plotSet != "NSSF")) +
    geom_point(aes(y = trees_ACD, x = CISG, col = plotSet)) +
        # scale_y_continuous(trans = "log10") +
        # scale_x_continuous(trans = "log10") +
    geom_abline(intercept = 0, slope = 1) +
    labs(y = "Kwek Yan's calculations", x = "from CISG team") +
    theme_classic()

ggplot(data = plots_AGBsum) +
    geom_point(aes(y = bigTrees, x = Ellie, col = plotSet)) +
    # scale_y_continuous(trans = "log10") +
    # scale_x_continuous(trans = "log10") +
    geom_abline(intercept = 0, slope = 1) +
    labs(y = "from NParks", x = "from Ellie") +
    theme_classic()

summary(lm(trees*0.45 ~ CISG, data = plots_AGBsum))
