# Load libraries ----

library(tidyverse)
library(BIOMASS)
library(ggplot2)

# Read LTFEM data ----

treeHeights <- list()
treeHeights$WYK <- 
    read.csv(paste0(sourceFilesDir, "WYK_treeHeights_2022.csv")) %>%
    mutate(height_avg = data.frame(height_1, height_2, height_3) %>%
               apply(1, mean)) %>%
    left_join(read.csv(paste0(sourceFilesDir, "WYK_trees_2022.csv")), by = c("plotID", "treeID"))
treeHeights$NatureParks <- 
    read.csv(paste0(sourceFilesDir, "NatureParks_treeHeights_2023.csv")) %>%
    mutate(height_avg = data.frame(height_1, height_2, height_3) %>%
               apply(1, mean)) %>%
    left_join(rbind(read.csv(paste0(sourceFilesDir, "NatureParks_trees_2023.csv")),
                    read.csv(paste0(sourceFilesDir, "NatureParks_bigTrees_2023.csv")) %>% rename(treeID = bigTreeID)),
              by = c("plotID", "treeID")) %>%
    mutate(DBH = as.numeric(DBH))
treeHeights$NSSF <- 
    read.csv(paste0(sourceFilesDir, "NSSF_treeHeights_2018.csv")) %>%
    mutate(height_avg = data.frame(height_1, height_2, height_3) %>%
               apply(1, mean, na.rm = TRUE)) %>%
    left_join(read.csv(paste0(sourceFilesDir, "NSSF_trees_2018.csv")), by = c("plotID", "treeID")) %>%
    mutate(DBH = as.numeric(DBH))


# fit allometric model ----

plotSets <- c("WYK", "NatureParks", "NSSF")

mods_HD <- lapply(plotSets, function(s) {
    with(treeHeights[[s]],
         modelHD(D = DBH,
                 H = height_avg,
                 method = "log1"))
})
names(mods_HD) <- plotSets

# plot graphs ----

DBH_minmax <- lapply(treeHeights, function(df) range(df$DBH, na.rm = TRUE))

mods_HD_pred <- lapply(plotSets, function(s) {
    D <- seq(DBH_minmax[[s]][1],
             DBH_minmax[[s]][2],
             length.out = 1000)
    
    cbind(D,
          apply(predict(mods_HD[[s]]$model,
                        newdata = list(D = D),
                        interval = "confidence"),
                2, exp))
})
names(mods_HD_pred) <- plotSets

ggplot() +
    theme_classic() +
    geom_point(data = treeHeights$WYK,
               mapping = aes(x = DBH, y = height_avg),
               color = "darkgreen") +
    geom_line(data = mods_HD_pred[["WYK"]],
              mapping = aes(x = D, y = fit),
              color = "darkgreen") +
    geom_point(data = treeHeights$NatureParks,
               mapping = aes(x = DBH, y = height_avg),
               color = "red") +
    geom_line(data = mods_HD_pred[["NatureParks"]],
              mapping = aes(x = D, y = fit),
              color = "red") +
    geom_point(data = treeHeights$NSSF,
               mapping = aes(x = DBH, y = height_avg),
               color = "blue") +
    geom_line(data = mods_HD_pred[["NSSF"]],
              mapping = aes(x = D, y = fit),
              color = "blue") +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(trans = "log10")
