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
-e --seed=<seed>                             random seed
-o --output-file=<output.file>               output file name
-y --toy                                     use toy data
-i --thin=<thin>                             thin
-q --sample-prior                            sample prior
-l --sample-observations                     sample observations
-d --daily                                   daily time steps
-a --deaths                                  fit to deaths
-k --keep                                    keep working directory
-f --force                                   force overwrite
-v --verbose                                 be verbose
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
daily <- opts[["daily"]]

library('dplyr')
library('tidyr')
library('RBi')
library('RBi.helpers')
library('cowplot')
library('stringi')

## set beginning and end of analysis
min_date <- as.Date("2014-06-01")
max_date <- as.Date("2014-10-20")

## set seed
if (length(seed) == 0) {
    seed <- ceiling(runif(1, -1, .Machine$integer.max - 1))
}

set.seed(seed)

rate_multiplier <- ifelse(daily, 1, 7)

## read incidence data
inc_filename <- paste0("incidence_", ifelse(daily, "daily", "weekly"), ".rds")
incidence <- readRDS(paste(code_dir, "data", inc_filename, sep = "/"))
## read admission delays
adm_filename <- paste0("admission_delays_",
                       ifelse(daily, "daily", "weekly"), ".rds")
admission_delays <- readRDS(paste(code_dir, "data", adm_filename, sep = "/"))
admission_delays <- admission_delays %>% 
  mutate(admission.delay = admission.delay)

admission_dates <-
  data.frame(date = seq.Date(min_date, max_date,
                             by = ifelse(daily, "day", "week")))
if (!daily) {
  admission_dates <- admission_dates %>%
    mutate(date = date - wday(date) + 2)
}
admissions_data <- admission_dates %>%
  left_join(incidence[["admissions"]], by = "date") %>%
  left_join(incidence[["deaths"]], by = "date") %>%
  filter(between(date, min_date, max_date)) %>%
  mutate(nr = as.integer(((date - min(date)) / rate_multiplier + 1))) %>% 
  mutate(admissions = ifelse(is.na(admissions), 0, admissions)) %>% 
  mutate(deaths = ifelse(is.na(deaths), 0, deaths))

obs <- list(Admissions = admissions_data %>% select(nr, value = admissions))

if (deaths) obs[["Deaths"]] <- admissions_data %>% select(nr, value = deaths)

delay_dates <-
  data_frame(date = seq.Date(min(admissions_data$date) - 7,
                             max(admissions_data$date),
                             by = ifelse(daily, "day", "week")))
delay_dates$value <-
  approx(x = admission_delays$date,
         y = admission_delays$admission.delay,
         xout = delay_dates$date)$y
delay_dates <- delay_dates %>% 
  mutate(nr = as.integer(((date - min(date)) / rate_multiplier)))

input <- list(admission_delay = delay_dates %>% select(nr, value))

if (length(output_file_name) == 0)
{
    filebase <- "ebola_lofa"
    output_file_name <- paste(filebase, r0_trajectory, sep = "_")
}

if (length(grep("/", output_file_name)) == 0)
{
    output_file_name <- paste(output_dir, output_file_name, sep = "/")
}

working_dir<- paste(output_file_name, sep = "/")
unlink(working_dir, recursive = TRUE)
dir.create(working_dir)

global_options <-
    list("end-time" = max(admissions_data$nr),
         "start-time" = 0,
         noutputs = max(admissions_data$nr), 
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

ebola_model$fix(rate_multiplier = rate_multiplier)

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
                         "R0 <- max(0, p_initR0 + n_R0_walk * p_vol_R0)")
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
                   working_folder = working_dir,
                   obs = obs, input = input, verbose = verbose)
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
                   input = input, obs = obs,
                   verbose = verbose)

cat(date(), "Running the stochastic model.\n")
run_prior <- adapt_mcmc(run_prior, min = 0, max = 1)

if (length(num_particles) > 0) {
  run_particle_adapted <- run_prior
} else
{
  libbi_seed <- ceiling(runif(1, -1, .Machine$integer.max - 1))
  run_prior$global_options[["seed"]] <- libbi_seed
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

run_adapted <- adapt_mcmc(run_particle_adapted, min = 0.05, max = 0.5,
                          max_iter = 10, scale = 2)

cat(date(), "Sampling from the posterior distribution of the full model.\n")
libbi_seed <- ceiling(runif(1, -1, .Machine$integer.max - 1))
run <- run_adapted$clone()
run$run(add_options = list("init-np" = pre_samples - 1,
                           nsamples = num_samples,
                           seed = libbi_seed),
        init = run_adapted)

if (length(model_file) == 0)
{
    model_file <- paste(output_file_name, "bi", sep = ".")
    run$model$write_model_file(model_file)
}

command_file <- paste(output_file_name, "cmd", sep = ".")
cat(run$result$command, file = command_file)

final_model <- run$model

saveRDS(list(nparticles = nparticles),
        file = paste0(output_file_name, "_args.rds"))

res <- bi_read(read = run, thin = thin,
               vars = c(final_model$get_vars("param"),
                        final_model$get_vars("noise"),
                        "R0", 
                        "loglikelihood", "logprior"),
               verbose = verbose)

cat(date(), "Plotting.\n")

burn <- 0.20 * num_samples
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

cat(date(), "..parameters.\n")
l <- lapply(names(res), function(x) {
    res[[x]] %>% mutate(state = x)
})

params <- rbind_all(l)
if ("nr" %in% names(params))
{
    params <- params %>%
        filter(is.na(nr)) %>%
        select(-nr)
}

saveRDS(params, paste0(output_file_name, "_params.rds"))

if (sample_obs)
{
    cat(date(), "Sampling from the joint distribution.\n")
    libbi_seed <- ceiling(runif(1, -1, .Machine$integer.max - 1))

    res_obs <-
      sample_observations(run,
                          read_options = list(thin = thin, verbose = verbose),
                          add_options = list(seed = libbi_seed))
    res_obs <- lapply(res, function(x)
    {
      if ("nr" %in% names(x))
      {
        x <- x %>% mutate(nr = nr - 1) %>% filter(nr >= 0)
      }
    })

    for (obs in names(res_obs))
    {
      res[[obs]] <- res_obs[[obs]]
    }

    data <- admissions_data %>%
      select(-nr) %>% 
##      gather(state, value, admissions:deaths) %>%
      gather(state, value, admissions) %>%
      mutate(state = stri_trans_totitle(state)) %>%
      rename(time = date)

    plot_args[["read"]] <- res
    plot_args[["data"]] <- data
    plot_args[["steps"]] <- FALSE
    plot_args[["states"]] <- names(res_obs)
    plot_args[["params"]] <- c()
    plot_args[["noises"]] <- c()
    plot_args[["hline"]] <- NULL

    p_obs <- do.call(plot_libbi, plot_args)

    if (!is.null(p_obs[["states"]]))
    {
      ggsave(paste(output_file_name, "states.pdf", sep = "_"), p_obs$states)
    }

    saveRDS(p_obs$data, paste0(output_file_name, "_obs_fits.rds"))
}

if (!keep) unlink(working_dir, recursive = TRUE)

quit()

