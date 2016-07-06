library('RBi')
library('data.table')
library('dplyr')

code_dir <- path.expand("~/code/ebola_lofa/")

## create working folder
output_folder <- path.expand("~/Data/Ebola/Lofa")
master_working_folder <- paste(output_folder, "libbi", sep = "/")
suppressWarnings(dir.create(master_working_folder))
working_folder <- paste(master_working_folder, "bounded", sep = "/")
unlink(working_folder, recursive = TRUE)
dir.create(working_folder)

## get posterior file
posterior <- readRDS("lofa_parallel_traces.rds")

## get original input file
min_date <- as.Date("2014-06-01")
max_date <- as.Date("2014-10-20")
rate_multiplier <- 7

## read incidence data
inc_filename <- paste0("incidence_weekly.rds")
incidence <- readRDS(paste(code_dir, "data", inc_filename, sep = "/"))

## read admission dates
adm_filename <- paste0("admission_delays_weekly.rds")
admission_delays <- readRDS(paste(code_dir, "data", adm_filename, sep = "/"))
admission_delays <- admission_delays %>%
  mutate(admission.delay = admission.delay)

admission_dates <-
  data.frame(date = seq.Date(min_date, max_date, by = "week"))
admissions_data <- admission_dates %>%
  left_join(incidence[["admissions"]], by = "date") %>%
  left_join(incidence[["deaths"]], by = "date") %>%
  filter(between(date, min_date, max_date)) %>%
  mutate(time = as.integer(((date - min(date)) / rate_multiplier + 1)) - 1) %>%
  mutate(admissions = ifelse(is.na(admissions), 0, admissions)) %>%
  mutate(deaths = ifelse(is.na(deaths), 0, deaths))

delay_dates <-
  data_frame(date = seq.Date(min(admissions_data$date),
                             max(admissions_data$date), by = "week"))
delay_dates$value <-
  approx(x = admission_delays$date,
         y = admission_delays$admission.delay,
         xout = delay_dates$date)$y
delay_dates <- delay_dates %>% 
  mutate(time = as.integer(((date - min(date)) / rate_multiplier)))

posterior$R0[, time := time - 1]
posterior$R0 <- posterior$R0[time >= 0]
posterior$H[, time := time - 1]
posterior$H <- posterior$H[time >= 0]
posterior$n_admission[, time := time - 1]
posterior$n_admission <- posterior$n_admission[time >= 0]


## read bed capacity
bed_filename <- paste0("lofa_etc_weekly.rds")
bed_availability <- readRDS(paste(code_dir, "data", bed_filename, sep = "/"))

bed_availability <- bed_availability[date >= min_date & date <= max_date]
bed_availability[, time := seq_len(nrow(bed_availability)) - 1]

## combine original input file with R0 trajectories from the posterior
input <- list(R0 = posterior$R0,
              H = posterior$H,
              K = bed_availability[, list(time, value = available.ebola)],
              admission_delay = delay_dates %>% select(time, value),
              late_increase = 1,
              n_admission = posterior$n_admission)

## get model
ebola_model <- bi_model(paste0(code_dir, "ebola_lofa_sim.bi"))

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

scenario <- 1

scenarios[[scenario]] <- copy(input)
scenarios[[scenario]][["K"]][, value := value[1]]

############################################################################
## scenario 2: healthcare seeking behaviour as in week 1
############################################################################

scenario <- 2

scenarios[[scenario]] <- copy(input)
H_week_1 <- scenarios[[scenario]][["H"]][time == 0, list(np = np, new_value = value)]
scenarios[[scenario]][["H"]] <- merge(scenarios[[scenario]][["H"]], H_week_1, by = c("np"), all.x = TRUE)
scenarios[[scenario]][["H"]] <- scenarios[[scenario]][["H"]][, list(time, np, value, new_value)]
setkey(scenarios[[scenario]][["H"]], time, np)
scenarios[[scenario]][["H"]][time > 0, value := new_value]
scenarios[[scenario]][["H"]][, new_value := NULL]

############################################################################
## scenario 3: behaviour as in week 1                                     ##
############################################################################

scenario <- 3

scenarios[[scenario]] <- copy(input)
R0_week_1 <- scenarios[[scenario]][["R0"]][time == 1, list(np = np, new_value = value)]
scenarios[[scenario]][["R0"]] <- merge(scenarios[[scenario]][["R0"]], R0_week_1, by = "np", all.x = TRUE)
scenarios[[scenario]][["R0"]] <- scenarios[[scenario]][["R0"]][, list(time, np, value, new_value)]
setkey(scenarios[[scenario]][["R0"]], time, np)
scenarios[[scenario]][["R0"]][time > 1, value := new_value]
scenarios[[scenario]][["R0"]][, new_value := NULL]
H_week_1 <- scenarios[[scenario]][["H"]][time == 1, list(np = np, new_value = value)]
scenarios[[scenario]][["H"]] <- merge(scenarios[[scenario]][["H"]], H_week_1, by = "np", all.x = TRUE)
scenarios[[scenario]][["H"]] <- scenarios[[scenario]][["H"]][, list(time, np, value, new_value)]
setkey(scenarios[[scenario]][["H"]], time, np)
scenarios[[scenario]][["H"]][time > 1, value := new_value]
scenarios[[scenario]][["H"]][, new_value := NULL]

############################################################################
## scenario 4: behaviour as in week 9                                    ##
############################################################################

scenario <- 4

scenarios[[scenario]] <- copy(input)
R0_week_9 <- scenarios[[scenario]][["R0"]][time == 9, list(np = np, new_value = value)]
scenarios[[scenario]][["R0"]] <- merge(scenarios[[scenario]][["R0"]], R0_week_9, by = "np", all.x = TRUE)
scenarios[[scenario]][["R0"]] <- scenarios[[scenario]][["R0"]][, list(time, np, value, new_value)]
setkey(scenarios[[scenario]][["R0"]], time, np)
scenarios[[scenario]][["R0"]][time > 9, value := new_value]
scenarios[[scenario]][["R0"]][, new_value := NULL]
H_week_9 <- scenarios[[scenario]][["H"]][time == 9, list(np = np, new_value = value)]
scenarios[[scenario]][["H"]] <- merge(scenarios[[scenario]][["H"]], H_week_9, by = "np", all.x = TRUE)
scenarios[[scenario]][["H"]] <- scenarios[[scenario]][["H"]][, list(time, np, value, new_value)]
setkey(scenarios[[scenario]][["H"]], time, np)
scenarios[[scenario]][["H"]][time > 9, value := new_value]
scenarios[[scenario]][["H"]][, new_value := NULL]

############################################################################
## scenario 5: bed expansion 2 weeks later
############################################################################

scenario <- 5

scenarios[[scenario]] <- copy(input)
scenarios[[scenario]][["K"]][, value := c(value[1:2], value[seq_len(length(value) - 2)])]

############################################################################
## scenario 6: exactly the same                                           ##
############################################################################

scenario <- 6

scenarios[[scenario]] <- copy(input)

############################################################################
## scenario 7: enough beds for everyone
############################################################################

scenario <- 7

scenarios[[scenario]] <- copy(input)
scenarios[[scenario]][["K"]][, value := 1e+7]

############################################################################
## simulate all scenarios                                                 ##
############################################################################

sim_scenarios <- list()
## for (scenario in seq_along(scenarios))
for (scenario in 7)
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

all_scenarios <- rbindlist(state_scenarios)

saveRDS(all_scenarios, "lofa_scenarios.rds")

admissions_no_change <- all_scenarios[scenario %in% c(6,7) & state == "Admissions",
                                      list(mean = mean(value),
                                           min.1 = quantile(value, 0.25),
                                           max.1 = quantile(value, 0.75),
                                           min.2 = quantile(value, 0.025),
                                           max.2 = quantile(value, 0.975)),
                                      by = list(scenario, time)]

admissions <- all_scenarios[state == "Admissions",
                            list(mean = mean(value),
                                 min.1 = quantile(value, 0.25),
                                 max.1 = quantile(value, 0.75),
                                 min.2 = quantile(value, 0.025),
                                 max.2 = quantile(value, 0.975)),
                            by = list(scenario, time)]

cases <- all_scenarios[state == "Zc",
                       list(mean = mean(value),
                            min.1 = quantile(value, 0.25),
                            max.1 = quantile(value, 0.75),
                            min.2 = quantile(value, 0.025),
                            max.2 = quantile(value, 0.975)),
                       by = list(scenario, time)]

r0 <- all_scenarios[state == "R0",
                    list(mean = mean(value),
                         min.1 = quantile(value, 0.25),
                         max.1 = quantile(value, 0.75),
                         min.2 = quantile(value, 0.025),
                         max.2 = quantile(value, 0.975)),
                    by = list(scenario, time)]
