library('RBi')
library('dplyr')
library('tidyr')
library('RBi.helpers')
library('truncnorm')
library('cowplot')

output_dir <- path.expand("~/Data/Ebola/Lofa")

n_traj <- 94

res <- list()

obs <- NULL
input <- NULL
data <- NULL

for (i in seq_len(n_traj))
{
    filename <- paste0(output_dir, "/ebola_lofa_independent_", i)
    res[[i]] <- readRDS(paste0(filename, ".rds"))
    model <- bi_model(paste0(filename, ".bi"))
    if (is.null(obs)) obs <- readRDS(paste0(filename, "_obs.rds"))
    if (is.null(input)) input <- readRDS(paste0(filename, "_input.rds"))
    if (is.null(data)) data <- readRDS(paste0(filename, "_data.rds"))
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
saveRDS(l, "lofa_parallel_traces.rds")

## est_params <- l[c("p_Inf", "p_R0", "p_vol_R0", "p_early_H", "p_late_H", "p_H_alpha", "p_H_tau")]
## est_params <- lapply(names(est_params), function(x) {est_params[[x]][, state := x]})
## lep <- rbindlist(est_params)
## mep <- data.table(dcast(lep, np ~ state))
## mep$np <- NULL

## trace <- mcmc(mep)

## sample admissions

l$Admissions <- copy(l$Zh)
l$Admissions[, mean := value]
l$Admissions[, sd := sqrt(value)]
l$Admissions[, value := rtruncnorm(n = .N, a = 0, mean = mean, sd = sd)]
l$Admissions[is.na(value), value := 0]

p <- plot_libbi(l, model, data = data, date.origin = as.Date("2014-06-02") - 7, date.unit = "week", data.colour = "black", densities = "histogram")

p_obs <- plot_libbi(l, model, data = data, date.origin = as.Date("2014-06-02") - 7, date.unit = "week", data.colour = "black", densities = "histogram", states = "Admissions", params = NULL, noises = NULL)
#p <- plot_libbi(l, model, density_args = list(adjust = 2))
save_plot("lofa_fit.pdf", p_obs$states + scale_y_continuous("Weekly new admissions"))

l$p_Inf[, list(mean = mean(value),
               min.50 = quantile(value, 0.25),
               max.50 = quantile(value, 0.75),
               min.95 = quantile(value, 0.025),
               min.95 = quantile(value, 0.975))]

l$p_R0[, list(mean = mean(value),
               min.50 = quantile(value, 0.25),
               max.50 = quantile(value, 0.75),
               min.95 = quantile(value, 0.025),
               min.95 = quantile(value, 0.975))]

l$p_early_H[, list(mean = mean(value),
                   min.50 = quantile(value, 0.25),
                   max.50 = quantile(value, 0.75),
                   min.95 = quantile(value, 0.025),
                   min.95 = quantile(value, 0.975))]

l$p_late_H[, list(mean = mean(value),
                  min.50 = quantile(value, 0.25),
                  max.50 = quantile(value, 0.75),
                  min.95 = quantile(value, 0.025),
                  min.95 = quantile(value, 0.975))]

l$p_vol_R0[, list(mean = mean(value),
                  min.50 = quantile(value, 0.25),
                  max.50 = quantile(value, 0.75),
                  min.95 = quantile(value, 0.025),
                  min.95 = quantile(value, 0.975))]

late_R0[, list(mean = mean(value),
               min.50 = quantile(value, 0.25),
               max.50 = quantile(value, 0.75),
               min.95 = quantile(value, 0.025),
               max.95 = quantile(value, 0.79))]

beta_params <- c("p_R0", "p_gamma", "p_alpha", "p_cfr", "H")

beta_list <- l[beta_params]
beta_list <- lapply(names(beta_list), function(x) {if (!("value" %in% colnames(beta_list[[x]]))) {setnames(beta_list[[x]], "value", x)}})
l <- lapply(names(l), function(x) {if (!("value" %in% colnames(l[[x]]))) {setnames(l[[x]], x, "value")}})

nbl <- copy(beta_list[[1]])
for (i in seq(2, length(beta_list)))
{
  nbl <- merge(nbl, beta_list[[i]], by = intersect(colnames(nbl), colnames(beta_list[[i]])))
}

nbl[, beta := p_R0 * p_gamma * p_alpha / (p_alpha + p_cfr * (1 - H) * p_gamma)]
nbl[, cR0 := beta * (1 / p_gamma + p_cfr/p_alpha)]

l$R0$R0 <- nbl$cR0
l$p_R0$value <- nbl$cR0

for (i in names(l)){
  setnames(l[[i]], i, "value")
}

elH <- l$H[, list(early = value[1], late = value[length(value)]), by = np]
elH[, ratio := late / early]
elH[, list(median = median(ratio),
           min.50 = quantile(ratio, 0.25),
           max.50 = quantile(ratio, 0.75),
           min.95 = quantile(ratio, 0.025),
           max.95 = quantile(ratio, 0.975))]

elH[, list(median = median(early),
           min.50 = quantile(early, 0.25),
           max.50 = quantile(early, 0.75),
           min.95 = quantile(early, 0.025),
           max.95 = quantile(early, 0.975))]

elH[, list(median = median(late),
           min.50 = quantile(late, 0.25),
           max.50 = quantile(late, 0.75),
           min.95 = quantile(late, 0.025),
           max.95 = quantile(late, 0.975))]

