library('RBi')
library('coda')

traces <- list()
H_traj <- list()
R_traj <- list()
for (i in 1:10)
{
    filename <-
        paste("~/Data/Ebola/Lofa/ebola_lofa_independent_5_poisson", i, sep = "_")
    res <- readRDS(paste(filename, "rds", sep = "."))
    model <- bi_model(paste(filename, "bi", sep = "."))
    traces[[i]] <- mcmc(get_traces(res, model = model))
    H_traj[[i]] <- res[["H"]] %>% mutate(value = value / 10)
    R_traj[[i]] <- res[["R0"]] %>% mutate(value = value / 10)
}

res <- list(H = bind_rows(H_traj),
            R0 = bind_rows(R_traj))

plot_libbi(res, model, steps = TRUE)

mcmc <- mcmc.list(traces)
