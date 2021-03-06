############################################################################
## Script for running libbi analysis                                      ##
############################################################################

library('docopt')

"Script for fitting the model to Ebola data from Lofa county.
Usage: run_libbi_lofa.r [options]
Options:
-n --nsamples=<num.samples>                  number of samples to obtain
-c --nparticles=<num.particles>              number of particles
-p --presamples=<pre.samples>                number of preparatory samples to obtain
-r --r0-trajectory=<trajectory>              R0 trajectory: exponential, bounded or independent
-t --threads=<num.threads>                   number of threads
-g --late-increase=<factor>                  late increase in transmission rate
-e --seed=<seed>                             random seed
-o --output-file=<output.file>               output file name
-y --toy                                     use toy data
-i --thin=<thin>                             thin
-q --sample-prior                            sample prior
-l --sample-observations                     sample observations
-a --deaths                                  fit to deaths
-m --model-file=<model.file>                 given model file (means there will be no adaptation step)
-k --keep                                    keep working directory
-f --force                                   force overwrite
-v --verbose                                 be verbose
-b --parallel-number=<number>                parallel number
-h --help                                    show this message" -> doc

opts <- docopt(doc)

if (opts[["help"]])
{
    print(opts)
    exit()
}

code_dir <- path.expand("~/code/ebola_lofa/")
output_dir <- path.expand("~/Data/Ebola/Lofa")

## mandatory arguments
num_samples <- as.integer(opts[["nsamples"]])
num_particles <- as.integer(opts[["nparticles"]])
pre_samples <- as.integer(opts[["presamples"]])
r0_trajectory <- opts[["r0-trajectory"]]
num_threads <- as.integer(opts[["threads"]])
late_increase <- as.integer(opts[["late-increase"]])
seed <- as.integer(opts[["seed"]])
output_file_name <- opts[["output-file"]]
thin <- as.integer(opts[["thin"]])
toy <- opts[["toy"]]
model_file <- opts[["model-file"]]
sample_obs <- opts[["sample-observations"]]
sample_prior <- opts[["sample-prior"]]
keep <- opts[["keep"]]
deaths <- opts[["deaths"]]
verbose <- opts[["verbose"]]
par_nb <- as.integer(opts[["parallel-number"]])

library('dplyr')
library('tidyr')
library('RBi')
library('RBi.helpers')
library('cowplot')
library('stringi')
library('truncnorm')

## set beginning and end of analysis
min_date <- as.Date("2014-06-01")
max_date <- as.Date("2014-10-20")

## set seed
if (length(seed) == 0) {
    seed <- ceiling(runif(1, -1, .Machine$integer.max - 1))
}

if (length(output_file_name) == 0)
{
    filebase <- "ebola_lofa"
    output_file_name <- paste0(filebase, "_", r0_trajectory, ifelse(length(par_nb) == 0, "", paste0("_", par_nb)))
}

if (length(grep("/", output_file_name)) == 0)
{
    output_file_name <- paste(output_dir, output_file_name, sep = "/")
}

set.seed(seed)

## read incidence data
inc_filename <- "lofa_incidence.rds"
incidence <- readRDS(paste(code_dir, "data", inc_filename, sep = "/"))
## read admission delays
adm_filename <- paste0("lofa_admission_delays.rds")
admission_delays <- readRDS(paste(code_dir, "data", adm_filename, sep = "/"))
admission_dates <-
  data.frame(date = seq.Date(min_date, max_date, by = "week"))

admission_dates <- admission_dates %>%
  mutate(date = date - wday(date) + 2)

admissions <- incidence[["admissions"]] %>%
  filter(classification %in% c("confirmed", "probable")) %>%
  group_by(date) %>%
  summarise(admissions = sum(admissions)) %>%
  ungroup

incident_deaths <- incidence[["deaths"]] %>%
  filter(classification %in% c("confirmed", "probable")) %>%
  group_by(date) %>%
  summarise(deaths = sum(deaths)) %>%
  ungroup

admissions_data <- admission_dates %>%
  left_join(admissions, by = "date") %>%
  left_join(incident_deaths, by = "date") %>%
  filter(between(date, min_date, max_date)) %>%
  mutate(week = as.integer(((date - min(date)) / 7 + 1))) %>%
  mutate(admissions = ifelse(is.na(admissions), 0, admissions)) %>%
  mutate(deaths = ifelse(is.na(deaths), 0, deaths))

obs <- list(Admissions = admissions_data %>% select(week, value = admissions))

saveRDS(obs, paste(output_file_name, "obs.rds", sep = "_"))

if (deaths) obs[["Deaths"]] <- admissions_data %>% select(week, value = deaths)

delay_dates <-
  data_frame(date = seq.Date(min(admissions_data$date) - 7,
                             max(admissions_data$date),
                             by = "week"))
delay_dates$value <-
  approx(x = admission_delays$date,
         y = admission_delays$median,
         xout = delay_dates$date)$y
delay_dates <- delay_dates %>% 
  mutate(week = as.integer(((date - min(date)) / 7)))

input <- list(admission_delay = delay_dates %>% select(week, value))

saveRDS(input, paste(output_file_name, "input.rds", sep = "_"))

working_dir<- paste(output_file_name, sep = "/")
unlink(working_dir, recursive = TRUE)
dir.create(working_dir)

global_options <-
    list("end-time" = max(admissions_data$week),
         "start-time" = 0,
         noutputs = max(admissions_data$week), 
         nsamples = pre_samples)

if (length(num_threads) > 0)
{
  global_options[["nthreads"]] <- num_threads
}

if (length(seed) > 0) global_options[["seed"]] <- seed

ebola_model <- bi_model(paste0(code_dir, "ebola_lofa_fit.bi"))
if (verbose) ## all states
{
  no_output_pattern <- "has_output[[:space:]]*=[[:space:]]*0"
  no_output <- grep(no_output_pattern, ebola_model$get_lines())
  updated_lines <- sub(no_output_pattern, "",
                       ebola_model$get_lines()[no_output])
  ebola_model$update_lines(no_output, updated_lines)
}

ebola_model$fix(rate_multiplier = 7)

init <- list()
if (length(late_increase) > 0)
{
    init <- c(init, list(late_increase = late_increase))
} else
{
    init <- c(init, list(late_increase = 1))
}

if (!deaths)
{
    ebola_model$fix(p_rep_d = 0)
}

## find line where we want to insert the R0 trajectory
transition_line <- grep("[[:space:]]*sub transition", ebola_model$get_lines())

if (r0_trajectory == "exponential")
{
    r0_update_lines <- c("n_R0_walk ~ wiener()",
                         "R0 <- R0 * exp(p_vol_R0 * n_R0_walk)")
} else if (r0_trajectory == "bounded")
{
    r0_update_lines <- c("n_R0_walk ~ wiener()",
                         "R0 <- max(0, R0 + p_vol_R0 * n_R0_walk)")
} else if (r0_trajectory == "independent")
{
    r0_update_lines <- c("n_R0_walk ~ gaussian()",
                         "R0 <- max(0, p_R0 + n_R0_walk * p_vol_R0)")
} else
{
    stop("'r0_trajectory' must be one of {exponential,bounded,independent)")
}

ebola_model$insert_lines(after = transition_line, r0_update_lines)

if (sample_prior)
{
    cat(date(), "Sampling from the prior distribution.\n")
    libbi_seed <- ceiling(runif(1, -1, .Machine$integer.max - 1))
    global_options[["seed"]] <- libbi_seed
    ## sample prior
    prior <- libbi(model = ebola_model, run = TRUE, target = "prior",
                   global_options = global_options, client = "sample",
                   working_folder = working_dir, time_dim = "week",
                   obs = obs, input = input, init = init,
                   verbose = verbose)
    ## reading
    res_prior <- bi_read(prior, vars = ebola_model$get_vars("param"),
                         verbose = verbose)
    saveRDS(res_prior, paste(output_file_name, "prior.rds", sep = "_"))
    prior_model_file <- paste(output_file_name, "prior.bi", sep = "_")
    prior$model$write_model_file(prior_model_file)
}


## ebola_deter <- ebola_model$fix_noise(list(n_R0_walk = 0))
## ebola_deter_prior <- ebola_deter$propose_prior()

## cat(date(), "Sampling from the posterior distribution of the deterministic model with prior = proposal.\n")
## libbi_seed <- ceiling(runif(1, -1, .Machine$integer.max - 1))
## global_options[["seed"]] <- libbi_seed
## global_options[["nsamples"]] <- pre_samples

## run_det_prior <- libbi(model = ebola_deter_prior, run = TRUE,
##                        obs = obs, input = input, 
##                        global_options = global_options, client = "sample",
##                        working_folder = working_dir, verbose = verbose)

## run_det_adapted <-
##   adapt_mcmc(run_det_prior, min = 0.1, max = 0.5, max_iter = 10, scale = 2)

## ebola_model$add_block("proposal_parameter",
##                       run_det_adapted$model$get_block("proposal_parameter"))
## ebola_model$add_block("proposal_initial",
##                       run_det_adapted$model$get_block("proposal_initial"))

min_particles <- 2**floor(log2(2 * nrow(admissions_data)))

if (length(num_particles) > 0)
{
  global_options[["nparticles"]] <- num_particles
} else
{
  global_options[["nparticles"]] <- min_particles
}

run_prior <- libbi(client = "sample", model = ebola_model,
                   global_options = global_options,
                   run = TRUE, working_folder = working_dir,
                   input = input, obs = obs, init = init,
                   time_dim = "week", verbose = verbose)

cat(date(), "Running the stochastic model.\n")
run_prior <- adapt_mcmc(run_prior, min = 0, max = 1)

if (length(num_particles) > 0) {
  run_particle_adapted <- run_prior
} else
{
  libbi_seed <- ceiling(runif(1, -1, .Machine$integer.max - 1))
  run_prior$global_options[["seed"]] <- libbi_seed
  cat(date(), "Adapting the number of particles.\n")
  run_particle_adapted <-
    adapt_particles(run_prior, min = 2 * nrow(admissions_data), max = 2**15)
}

if ("nparticles" %in% names(run_particle_adapted$global_options))
{
    nparticles <- run_particle_adapted$global_options[["nparticles"]]
} else
{
    nparticles <- 1
}

if (length(model_file) == 0)
{
  run_adapted <- adapt_mcmc(run_particle_adapted, min = 0.05, max = 0.4,
                            max_iter = 10, scale = 2, correlations = TRUE)
} else {
  run_adapted <- run_particle_adapted
}

cat(date(), "Sampling from the posterior distribution of the full model.\n")
libbi_seed <- ceiling(runif(1, -1, .Machine$integer.max - 1))
run <- run_adapted$clone()
run$run(add_options = list("init-np" = pre_samples - 1,
                           nsamples = num_samples,
                           seed = libbi_seed),
        init = run_adapted)

if (length(model_file) == 0)
{
    model_file_name <- paste(output_file_name, "bi", sep = ".")
    run$model$write_model_file(model_file_name)
}

command_file <- paste(output_file_name, "cmd", sep = ".")
cat(run$result$command, file = command_file)

final_model <- run$model

saveRDS(list(nparticles = nparticles),
        file = paste0(output_file_name, "_args.rds"))

read_options <- list(read = run, 
                     vars = c(final_model$get_vars("param"),
                              final_model$get_vars("noise"),
                              "Zh", "R0", "H",
                              "loglikelihood", "logprior"),
                     verbose = verbose)
if (length(thin) > 0) read_options[["thin"]] <- thin

res <- do.call(bi_read, read_options)

if (length(par_nb) == 0)
{
  cat(date(), "Plotting.\n")

  burn <- 0.20 * num_samples / ifelse(length(thin) > 0, thin, 1)
  plot_args <- list(read = res, model = final_model,
                    density_args = list(adjust = 2), burn = burn,
                    date.origin = as.Date("2014-06-02"),
                    hline = c(R0 = 1), steps = TRUE, date.unit = "week",
                    trend = "mean", plot = FALSE)
  if (sample_prior)
  {
    plot_args[["prior"]] <- res_prior
  }
  p_param <- do.call(plot_libbi, plot_args)
  saveRDS(p_param$data, paste0(output_file_name, "_param_fits.rds"))

  if (!is.null(p_param[["densities"]]))
  {
    ggsave(paste(output_file_name, "densities.pdf", sep = "_"), p_param$densities)
  }
  if (!is.null(p_param[["traces"]]))
  {
    ggsave(paste(output_file_name, "traces.pdf", sep = "_"), p_param$traces)
  }
  if (!is.null(p_param[["correlations"]]))
  {
    ggsave(paste(output_file_name, "correlations.pdf", sep = "_"),
           p_param$correlations)
  }
  if (!is.null(p_param[["noises"]]))
  {
    ggsave(paste(output_file_name, "noises.pdf", sep = "_"), p_param$noises)
  }
  if (!is.null(p_param[["likelihoods"]]))
  {
    ggsave(paste(output_file_name, "likelihoods.pdf", sep = "_"), p_param$likelihoods)
  }
  if (!is.null(p_param[["states"]]))
  {
    ggsave(paste(output_file_name, "r0.pdf", sep = "_"), p_param$states)
  }

  cat(date(), "..parameters.\n")
  l <- lapply(names(res), function(x) {
    res[[x]] %>% mutate(state = x)
  })

  params <- rbind_all(l)
  if ("week" %in% names(params))
  {
    params <- params %>%
      filter(is.na(week)) %>%
      select(-week)
  }

  saveRDS(params, paste0(output_file_name, "_params.rds"))

  if (sample_obs)
  {
    cat(date(), "Sampling from the joint distribution.\n")
    ## libbi_seed <- ceiling(runif(1, -1, .Machine$integer.max - 1))

    ## sample_observations(run,
    ##                     read_options = list(thin = thin, verbose = verbose),
    ##                     add_options = list(seed = libbi_seed))
    ## res_obs <- lapply(res, function(x)
    ## {
    ##   if ("nr" %in% names(x))
    ##   {
    ##     x <- x %>% mutate(nr = nr - 1) %>% filter(nr >= 0)
    ##   }
    ## })

    ## for (obs in names(res_obs))
    ## {
    ##   res[[obs]] <- res_obs[[obs]]
    ## }

    res$Admissions <- res$Zh %>%
      rename(Zh = value) %>%
      mutate(mean = Zh, sd = sqrt(Zh)) %>%
      mutate(sd = ifelse(sd < 1, 1, sd)) %>%
      mutate(value = rtruncnorm(n(), 0, mean = mean, sd = sd)) %>%
      mutate(time = time - 1) %>%
      filter(time >= 0)

    if (deaths)
    {
      res$Deaths <- res$Zd %>%
        rename(Zd = value) %>%
        left_join(res$p_rep_d %>% rename(rep_d = value), by = "np") %>%
        mutate(mean = Zd, sd = sqrt(rep_d * (1 - rep_d) * Zd)) %>%
        mutate(sd = ifelse(sd < 1, 1, sd)) %>%
        mutate(sd = ifelse(sd < 1, 1, sd)) %>%
        mutate(value = rtruncnorm(n(), 0, mean = mean, sd = sd))
    }

    data <- admissions_data %>%
      select(-week) %>%
      ##      gather(state, value, admissions:deaths) %>%
      gather(state, value, admissions) %>%
      mutate(state = stri_trans_totitle(state)) %>%
      rename(time = date)

    saveRDS(data, paste(output_file_name, "data.rds", sep = "_"))

    plot_args[["read"]] <- res
    plot_args[["data"]] <- data
    plot_args[["steps"]] <- FALSE
    plot_args[["params"]] <- c()
    plot_args[["noises"]] <- c()
    plot_args[["hline"]] <- NULL
    plot_args[["limit.to.data"]] <- TRUE

    p_obs <- do.call(plot_libbi, plot_args)

    if (!is.null(p_obs[["states"]]))
    {
      ggsave(paste(output_file_name, "states.pdf", sep = "_"), p_obs$states)
    }

    saveRDS(p_obs$data, paste0(output_file_name, "_obs_fits.rds"))
  }
} else
{
  saveRDS(res, paste0(output_file_name, ".rds"))
}

if (!keep) unlink(working_dir, recursive = TRUE)

quit()

