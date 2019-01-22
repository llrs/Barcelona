library("ggplot2")
library("gganimate")

models <- list.files(pattern = "models[0-9].?.RDS")
models <- lapply(models, readRDS)
models <- do.call(c, models)
out <- lapply(names(models), function(x) {
  cbind.data.frame("RNAseq" = models[[x]]$Y[[1]][, 1],
                   "Micro" = models[[x]]$Y[[2]][, 1],
                   model = x,
                   AVE_inner = models[[x]]$AVE$AVE_inner[[1]],
                   AVE_outer = models[[x]]$AVE$AVE_outer[1],
                   Sample = rownames(models[[x]]$Y[[1]]))
})

out2 <- Reduce(rbind, out)
out3 <- merge(out2, A$Meta, by.x = "Sample", by.y = "Original", all.x = TRUE,
              sort = TRUE)
out3$Interaction <- ifelse(grepl("i$", out3$model), 1, 0)
out3$model <- gsub("i$", "", out3$model)
out3$model <- gsub("model", "", out3$model)
out3$model <- gsub("b$", " best", out3$model)
theme_update(strip.background = element_blank())

out4 <- out3[out3$Interaction == 0, ]
out4 <- droplevels(out4)

basic <- ggplot(out4, aes(RNAseq, Micro, col = IBD, group = Sample)) +
  geom_point()
basic + geom_line(colour = 'grey')
anim <- basic +
  transition_states(model, wrap = FALSE) +
  view_follow(fixed_x = c(-0.35, 0.2), fixed_y = c(-0.6, 1)) +
  ggtitle('Model {closest_state}')
animate(anim)
