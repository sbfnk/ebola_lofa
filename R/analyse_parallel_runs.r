library('RBi')
library('dplyr')
library('tidyr')
library('RBi.helpers')
library('truncnorm')
library('cowplot')

output_dir <- path.expand("~/Data/Ebola/Lofa")

n_traj <- 100

res <- list()

obs <- NULL
input <- NULL
data <- NULL

for (i in seq_len(n_traj))
{
    filename <- paste0(output_dir, "/ebola_lofa_", i)
    res[[i]] <- readRDS(paste0(filename, ".rds"))
    model <- bi_model(paste0(filename, ".bi"))
    if (is.null(obs)) obs <- readRDS(paste0(filename, "_obs.rds"))
    if (is.null(input)) input <- readRDS(paste0(filename, "_input.rds"))
    if (is.null(data)) data <- readRDS(paste0(filename, "_data.rds"))
}

## combined into one chain
np_translate <- NULL
combined <- lapply(names(res[[1]]), function(x) {
  z <- rbindlist(lapply(seq_along(res), function(y) {
    data.table(res[[y]][[x]])[, unique_np := paste(np, y, sep = "_")]
  }))
  if (is.null(np_translate)) {
    np_translate <-
      data.table(unique_np = unique(z$unique_np),
                 new_np = seq_along(unique(z$unique_np)) - 1)
  }
  z <- merge(z, np_translate, by = "unique_np")
  z[, np := NULL]
  z[, unique_np := NULL]
  setnames(z, "new_np", "np")
  setkey(z, np)
  z
})

names(combined) <- names(res[[1]])
saveRDS(combined, "lofa_combined_traces.rds")

## est_params <- l[c("p_Inf", "p_R0", "p_vol_R0", "p_early_H", "p_late_H", "p_H_alpha", "p_H_tau")]
## est_params <- lapply(names(est_params), function(x) {est_params[[x]][, state := x]})
## lep <- rbindlist(est_params)
## mep <- data.table(dcast(lep, np ~ state))
## mep$np <- NULL

## trace <- mcmc(mep)

## sample admissions

combined$Admissions <- copy(combined$Zh)
combined$Admissions[, mean := value]
combined$Admissions[, sd := sqrt(value)]
combined$Admissions[, value := rtruncnorm(n = .N, a = 0, mean = mean, sd = sd)]
combined$Admissions[is.na(value), value := 0]

p <- plot_libbi(combined, model, data = data, date.origin = as.Date("2014-06-02") - 7, date.unit = "week", data.colour = "black", densities = "histogram", plot = FALSE)

combined$p_Inf[, list(mean = mean(value),
               min.50 = quantile(value, 0.25),
               max.50 = quantile(value, 0.75),
               min.95 = quantile(value, 0.025),
               min.95 = quantile(value, 0.975))]

combined$p_R0[, list(mean = mean(value),
               min.50 = quantile(value, 0.25),
               max.50 = quantile(value, 0.75),
               min.95 = quantile(value, 0.025),
               min.95 = quantile(value, 0.975))]

combined$p_early_H[, list(mean = mean(value),
                   min.50 = quantile(value, 0.25),
                   max.50 = quantile(value, 0.75),
                   min.95 = quantile(value, 0.025),
                   min.95 = quantile(value, 0.975))]

combined$p_late_H[, list(mean = mean(value),
                  min.50 = quantile(value, 0.25),
                  max.50 = quantile(value, 0.75),
                  min.95 = quantile(value, 0.025),
                  min.95 = quantile(value, 0.975))]

combined$p_vol_R0[, list(mean = mean(value),
                  min.50 = quantile(value, 0.25),
                  max.50 = quantile(value, 0.75),
                  min.95 = quantile(value, 0.025),
                  min.95 = quantile(value, 0.975))]

beta_params <- c("p_R0", "p_gamma", "p_alpha", "p_cfr", "H")

beta_list <- combined[beta_params]
beta_list <- lapply(names(beta_list), function(x) {if (!("value" %in% colnames(beta_list[[x]]))) {setnames(beta_list[[x]], "value", x)}})
combined <- lapply(names(combined), function(x) {if (!("value" %in% colnames(l[[x]]))) {setnames(l[[x]], x, "value")}})

nbl <- copy(beta_list[[1]])
for (i in seq(2, length(beta_list)))
{
  nbl <- merge(nbl, beta_list[[i]], by = intersect(colnames(nbl), colnames(beta_list[[i]])))
}

nbl[, beta := p_R0 * p_gamma * p_alpha / (p_alpha + p_cfr * (1 - H) * p_gamma)]
nbl[, cR0 := beta * (1 / p_gamma + p_cfr/p_alpha)]

combined$R0$R0 <- nbl$cR0
combined$p_R0$value <- nbl$cR0

for (i in names(combined)){
  setnames(combined[[i]], i, "value")
}

elH <- combined$H[, list(early = value[1], late = value[length(value)]), by = np]
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

