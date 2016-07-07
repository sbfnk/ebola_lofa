library('RBi')
library('coda')
library('dplyr')
library('RBi.helpers')

n_traj <- 20

traces <- list()
H_traj <- list()
R_traj <- list()
res <- list()

for (i in seq_len(n_traj))
{
    filename <-
        paste("~/Data/Ebola/Lofa/ebola_lofa_independent_poisson", i, sep = "_")
    res[[i]] <- readRDS(paste(filename, "rds", sep = "."))
    model <- bi_model(paste(filename, "bi", sep = "."))
 }

z_np_translate <- NULL
l <- lapply(names(res[[1]]), function(x) {
  z <- rbindlist(lapply(seq_along(res), function(y) {
    data.table(res[[y]][[x]])[, unique_np := paste(np, y, sep = "_")]
  }))
  if (is.null(z_np_translate)) {
    z_np_translate <-
      data.table(unique_np = unique(z$unique_np),
                 new_np = seq_along(unique(z$unique_np)) - 1)
  }
  z <- merge(z, z_np_translate, by = "unique_np")
  z[, np := NULL]
  z[, unique_np := NULL]
  setnames(z, "new_np", "np")
  setkey(z, np)
  z
})

names(l) <- names(res[[1]])

plot_libbi(l, model, alpha = 0.5, id = 0:99)
plot_libbi(l, model, id = 0:99, trend = NULL, quantile = NULL)

names(l) <- names(res[[1]])
saveRDS(l, "lofa_parallel_traces.rds")

mcmc <- mcmc.list(traces)

dt <- rbindlist(lapply(mcmc, data.table))

