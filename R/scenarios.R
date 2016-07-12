library('RBi')
library('data.table')

code_dir <- "~/code/ebola_lofa"
output_folder <- path.expand("~/Data/Ebola/Lofa")

## create working folder
working_folder <- paste(output_folder, "Scenarios", sep = "/")
unlink(working_folder, recursive = TRUE)
dir.create(working_folder)

## get posterior file
combined <-
  readRDS(paste(output_folder, "lofa_combined_traces.rds", sep = "/"))

posterior <- copy(combined$traces)

posterior$R0[, time := time - 1]
posterior$R0 <- posterior$R0[time >= 0]
posterior$H[, time := time - 1]
posterior$H <- posterior$H[time >= 0]
posterior$n_admission[, time := time - 1]
posterior$n_admission <- posterior$n_admission[time >= 0]

## read bed capacity
bed_filename <- "lofa_etc_weekly.rds"
bed_availability <- readRDS(paste(code_dir, "data", bed_filename, sep = "/"))

bed_availability <- bed_availability[date >= min(combined$data$time) &
                                     date <= max(combined$data$time)]
bed_availability[, time := seq_len(nrow(bed_availability)) - 1]

## combine original input file with R0 trajectories from the posterior
input <- c(list(R0 = posterior$R0,
                H = posterior$H,
                K = bed_availability[, list(time, value = available.ebola)],
                late_increase = 1,
                admission_delay = data.table(combined$input$admission_delay)[, list(time = week, value)],
                n_admission = posterior$n_admission))

## get model
ebola_model <- bi_model(paste(code_dir, "ebola_lofa_sim.bi", sep = "/"))

## remove R0 assign line, this is now from the input file
r0_assign_line <- grep("^[[:space:]]*R0[[:space:]]*<-", ebola_model$get_lines())
ebola_model$remove_lines(r0_assign_line)
ebola_model$remove_block("parameter")
ebola_model$fix(p_delta = 1 / (3 * 2.302684) * 7,
                p_theta = 1 / (10 * 1.279192) * 7)

## remove R0 from init file
posterior[["R0"]] <- NULL
posterior[["n_R0_walk"]] <- NULL
posterior[["H"]] <- NULL
posterior[["Zh"]] <- NULL
posterior[["n_admission"]] <- NULL

scenarios <- list()

############################################################################
## scenario 1: 10 beds in the hospital
############################################################################

scenario <- "10 beds"

scenarios[[scenario]] <- copy(input)
scenarios[[scenario]][["K"]][, value := 10]

############################################################################
## scenario 2: 40 beds in the hospital
############################################################################

scenario <- "40 beds"

scenarios[[scenario]] <- copy(input)
scenarios[[scenario]][["K"]][, value := 40]

############################################################################
## scenario 3: no troughs in R0
############################################################################

scenario <- "R0 always > 1"
scenarios[[scenario]] <- copy(input)
scenarios[[scenario]][["R0"]][value < 1,  value := 1]

############################################################################
## scenario 4: healthcare seeking behaviour as in week 1
############################################################################

scenario <- "No improvement in\nhealthcare seeking"

scenarios[[scenario]] <- copy(input)
H_week_1 <- scenarios[[scenario]][["H"]][time == 0, list(np = np, new_value = value)]
scenarios[[scenario]][["H"]] <- merge(scenarios[[scenario]][["H"]], H_week_1, by = c("np"), all.x = TRUE)
scenarios[[scenario]][["H"]] <- scenarios[[scenario]][["H"]][, list(time, np, value, new_value)]
setkey(scenarios[[scenario]][["H"]], time, np)
scenarios[[scenario]][["H"]][time > 0, value := new_value]
scenarios[[scenario]][["H"]][, new_value := NULL]

############################################################################
## scenario 5: exactly the same                                           ##
############################################################################

scenario <- "given beds"

scenarios[[scenario]] <- copy(input)

############################################################################
## scenario 6: enough beds for everyone
############################################################################

scenario <- "as observed"

scenarios[[scenario]] <- copy(input)
scenarios[[scenario]][["K"]][, value := 1e+7]

############################################################################
## simulate all scenarios                                                 ##
############################################################################

sim_scenarios <- list()
for (scenario in seq_along(scenarios))
{
    cat("Scenario ", scenario, "\n")
    global_options <-
        list("end-time" = 21,
             "start-time" = 0,
             noutputs = 21,
             nsamples = max(input$H$np) + 1,
             target = "joint")
    system.time({
        run_joint <- libbi(client = "sample", model = ebola_model,
                           global_options = global_options, run = TRUE,
                           working_folder = working_folder,
                           input = scenarios[[scenario]],
                           init = posterior)
    })
    sim_scenarios[[scenario]] <- bi_read(run_joint)
}

## sum over Erlang compartments
state_scenarios <- list()
for (i in seq_along(sim_scenarios)) {
    for (j in names(sim_scenarios[[i]])) {
      sim_scenarios[[i]][[j]] <- data.table(sim_scenarios[[i]][[j]])
        if ("time" %in% colnames(sim_scenarios[[i]][[j]])) {
            sim_scenarios[[i]][[j]] <- sim_scenarios[[i]][[j]][, list(value = sum(value)), by = list(time, np)]
            state_scenarios <- c(state_scenarios, list(sim_scenarios[[i]][[j]]))
        }
        sim_scenarios[[i]][[j]][, scenario := i]
        sim_scenarios[[i]][[j]][, state := j]
    }
}

## merge tables
all_scenarios <- rbindlist(state_scenarios)

saveRDS(all_scenarios, paste(output_folder, "lofa_scenarios.rds", sep = "/"))
all_scenarios <- readRDS(paste(output_folder, "lofa_scenarios.rds", sep = "/"))

incidences <- all_scenarios[state %in% c("Zc", "Zh", "Admissions"),
                            list(value = diff(value),
                                   time = time[-1]),
                            by = list(scenario, np, state)]
incidences[, date := as.Date("2014-06-02") + 7 * (time - 1)]

## for checking that admissions number look reasonable when running with actual bed numbers, or enough beds for everybody
admissions_no_change <-
  incidences[scenario %in% c(5,6) & state == "Admissions",
             list(mean = mean(value),
                  median = median(value),
                  min.1 = quantile(value, 0.25),
                  max.1 = quantile(value, 0.75),
                  min.2 = quantile(value, 0.025),
                  max.2 = quantile(value, 0.975)),
             by = list(scenario, date)]

## ETC admissions
admissions <- incidences[state == "Zh",
                            list(mean = mean(value),
                                 median = median(value),
                                 min.1 = quantile(value, 0.25),
                                 max.1 = quantile(value, 0.75),
                                 min.2 = quantile(value, 0.025),
                                 max.2 = quantile(value, 0.975)),
                            by = list(scenario, date)]

cases <- incidences[state == "Zc",
                       list(mean = mean(value),
                            median = median(value),
                            min.1 = quantile(value, 0.25),
                            max.1 = quantile(value, 0.75),
                            min.2 = quantile(value, 0.025),
                            max.2 = quantile(value, 0.975)),
                       by = list(scenario, date)]

final_sizes <- incidences[state == "Zc", list(value = sum(value)),
                          by = list(scenario, np)]
dfs <- dcast(final_sizes, np ~ scenario)

quantile(dfs[, `10 beds` / `as observed`], c(0.025, 0.25, 0.5, 0.75, 0.975))
quantile(dfs[, `40 beds` / `as observed`], c(0.025, 0.25, 0.5, 0.75, 0.975))
quantile(dfs[, `R0 always > 1` / `as observed`],
         c(0.025, 0.25, 0.5, 0.75, 0.975))
quantile(dfs[, `No improvement in\nhealthcare seeking` / `as observed`],
         c(0.025, 0.25, 0.5, 0.75, 0.975))

first_scenarios <- c("10 beds", "40 beds", "as observed")
second_scenarios <- c("R0 always > 1", "No improvement in\nhealthcare seeking", "as observed")

scenarios_1 <- cases[scenario %in% first_scenarios]
scenarios_1[, scenario := factor(scenario, first_scenarios)]
scenarios_2 <- cases[scenario %in% second_scenarios]
scenarios_2[, scenario := factor(scenario, second_scenarios)]

p1 <- ggplot(scenarios_1,
             aes(x = date, y = median, color = scenario, fill = scenario)) +
  geom_line() +
  geom_ribbon(aes(ymin = min.1, ymax = max.1), alpha = 0.5) +
  scale_color_brewer("", palette = "Set1") +
  scale_fill_brewer("", palette = "Set1") +
  scale_y_continuous("Incidence") +
  scale_x_date("") +
  theme(legend.position = "top")

p2 <- ggplot(scenarios_2,
            aes(x = date, y = median, color = scenario, fill = scenario)) +
  geom_line() +
  geom_ribbon(aes(ymin = min.1, ymax = max.1), alpha = 0.5) +
  scale_color_brewer("", palette = "Set1") +
  scale_fill_brewer("", palette = "Set1") +
  scale_y_continuous("Incidence") +
  scale_x_date("") +
  theme(legend.position = "top") + 
  coord_cartesian(ylim = c(0, 2000)) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

pg <- plot_grid(p1, p2, labels = c("a", "b"), align = "h")

save_plot("scenarios.pdf", pg, nrow = 1, ncol = 2, base_aspect_ratio = 1.3)
