library('RBi')
library('coda')
library('dplyr')
library('tidyr')
library('RBi.helpers')
library('stringi')
library('truncnorm')
library('cowplot')

n_traj <- 94

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

#plot_libbi(l, model, alpha = 0.5, id = 0:99)
#plot_libbi(l, model, id = 0:99, trend = NULL, quantile = NULL)

names(l) <- names(res[[1]])
saveRDS(l, "lofa_parallel_traces.rds")

est_params <- l[c("p_Inf", "p_R0", "p_vol_R0", "p_early_H", "p_late_H", "p_H_alpha", "p_H_tau")]
est_params <- lapply(names(est_params), function(x) {est_params[[x]][, state := x]})
lep <- rbindlist(est_params)
mep <- data.table(dcast(lep, np ~ state))
mep$np <- NULL

trace <- mcmc(mep)

code_dir <- path.expand("~/code/ebola_lofa/")
output_dir <- path.expand("~/Data/Ebola/Lofa")
min_date <- as.Date("2014-06-01")
max_date <- as.Date("2014-10-20")

inc_filename <- "lofa_incidence.rds"
incidence <- readRDS(paste(code_dir, "data", inc_filename, sep = "/"))

rate_multiplier <- 7
admissions <- incidence[["admissions"]] %>%
  filter(classification %in% c("confirmed", "probable")) %>%
  group_by(date) %>%
  summarise(admissions = sum(admissions)) %>%
  ungroup

admission_dates <-
  data.frame(date = seq.Date(min_date, max_date, by = "week"))

admission_dates <- admission_dates %>%
  mutate(date = date - wday(date) + 2)

admissions_data <- admission_dates %>%
  left_join(admissions, by = "date") %>%
  filter(between(date, min_date, max_date)) %>%
  mutate(week = as.integer(((date - min(date)) / rate_multiplier + 1))) %>%
  mutate(admissions = ifelse(is.na(admissions), 0, admissions)) 

data <- admissions_data %>%
  select(-week) %>%
  ##      gather(state, value, admissions:deaths) %>%
  gather(state, value, admissions) %>%
  mutate(state = stri_trans_totitle(state)) %>%
  rename(time = date)

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

