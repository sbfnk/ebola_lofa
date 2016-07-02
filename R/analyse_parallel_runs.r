library('RBi')
library('coda')

traces <- list()
H_traj <- list()
R_traj <- list()
for (i in 1:10)
{
    filename <-
        paste("~/Data/Ebola/Lofa/ebola_lofa_independent_poisson", i, sep = "_")
    res <- readRDS(paste(filename, "rds", sep = "."))
    model <- bi_model(paste(filename, "bi", sep = "."))
    traces[[i]] <- mcmc(get_traces(res, model = model))
    H_traj[[i]] <- res[["H"]] %>% mutate(value = value / 10)
    R_traj[[i]] <- res[["R0"]] %>% mutate(value = value / 10)
}

res <- list(H = bind_rows(H_traj),
            R0 = bind_rows(R_traj))

plot_libbi(res, model, steps = TRUE)

traj_join <- bind_rows(traj) %>% filter(np > 2499) %>% 
    group_by(time, state) %>%
    summarise(mean = mean(value),
              min.1 = quantile(value, 0.25),
              max.1 = quantile(value, 0.75),
              min.2 = quantile(value, 0.025),
              max.2 = quantile(value, 0.975)) %>%
    ungroup

p <- ggplot(traj_join, aes(x = time, y = mean)) +
    geom_step() +
    geom_ribbon(aes(ymin = min.1, ymax = max.1), alpha = 0.5) + 
    geom_ribbon(aes(ymin = min.2, ymax = max.2), alpha = 0.25) + 
    facet_wrap(~ state, scales = "free")

mcmc <- mcmc.list(traces)
