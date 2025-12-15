library(ENMeval)
library(raster)
library(dplyr)
library(ecospat)

set.seed(100)

records<-read.table("F:/Il mio Drive/Viticoltura/Cryptoblabes distribution/Cryptoblabes gnidiella GBIF/DB_Cryptoblabes.csv",sep=";",head=T)

#getting presence coordinates fromd data
occs<-records[,c(2,3)]
envs<-selected_bios
envs.files <- list.files(path=paste(system.file(package='dismo'), '/ex', sep=''), 
                         pattern='grd', full.names=TRUE)

occs.cells <- raster::extract(envs[[1]], occs, cellnumbers = TRUE)
occs.cellDups <- duplicated(occs.cells[,1])
occs <- occs[!occs.cellDups,]
plot(envs[[1]], main="Mean diurnal range")

points(occs)
points(occs, col = 'red', pch=19)

occs.z <- raster::extract(envs, occs)
occs.sim <- similarity(envs, occs.z)
occs.mess <- occs.sim$similarity_min

occs.sp <- sp::SpatialPoints(occs)

rasterVis::levelplot(occs.mess, main = "Environmental similarity", margin = FALSE) + 
  latticeExtra::layer(sp.points(occs.sp, col="black"))

myScale <- seq(cellStats(occs.mess, min), cellStats(occs.mess, max), length.out = 100)
rasterVis::levelplot(occs.mess, main = "Environmental similarity", at = myScale, margin = FALSE) + 
  latticeExtra::layer(sp.points(occs.sp, col="black"))

occs.sf <- sf::st_as_sf(occs, coords = c("x","y"),
                        crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
crs(envs) <- raster::crs(occs.sf)

occs.buf <- sf::st_buffer(occs.sf, dist = 600000) %>% sf::st_union() %>% sf::st_sf()
plot(envs[[1]], main = names(envs)[1])
points(occs)
plot(occs.buf, border = "blue", lwd = 3, add = TRUE)

envs.bg <- raster::crop(envs, occs.buf)
plot(envs.bg[[1]], main = names(envs)[1])
points(occs)
plot(occs.buf, border = "blue", lwd = 3, add = TRUE)
bg <- dismo::randomPoints(envs.bg[[1]], n = 500) %>% as.data.frame()
colnames(bg) <- colnames(occs)
plot(envs.bg[[1]])
points(bg, pch = 20, cex = 0.2)

jack <- get.jackknife(occs, bg)

evalplot.grps(pts = occs, pts.grp = jack$occs.grp, envs = envs.bg)

e.mx <- ENMevaluate(occs = occs, envs = envs, bg = bg, 
                    algorithm = 'maxnet', partitions = 'jackknife', 
                    tune.args = list(fc = c("L","H","LQ","LH","LQH"), rm = 1:8))

res <- eval.results(e.mx)
opt.aicc <- res %>% filter(delta.AICc == 0)
opt.aicc

opt.seq <- res %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))
opt.seq

bestmod = which(e.mx@results$AICc==min(e.mx@results$AICc))
e.mx@results[bestmod,]

# Predictions
pr = predict(envs, e.mx@models[[bestmod]], type = 'cloglog')
pr_df = as.data.frame(pr, xy=T)


write.table(res,"H:/Il mio Drive/SDM-HDM Toumeyella parvicornis/Modelling_in_R/Results/model_selection_ENMeval.csv",sep = ";")

#heatmap
library(ggplot2)
ggplot() +
  geom_raster(data = pr_df, aes(x = x, y = y, fill = layer)) +
  coord_quickmap() +
  theme_bw() + 
  scale_fill_viridis_c(option="viridis",
                       na.value = "white")
