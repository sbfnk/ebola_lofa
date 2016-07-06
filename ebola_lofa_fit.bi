model ebola_lofa_fit {

  const first_obs = 0
  const rate_multiplier = 1
 
  const e_rho = 2
  const e_gamma = 3
  const e_kappa = 2

  input admission_delay

  dim rho_erlang(e_rho)
  dim gamma_erlang(e_gamma)
  dim kappa_erlang(e_kappa)

  param p_N // population size
  param p_rho // rate of developing symptoms
  param p_gamma // recovery/death rate in the community
  param p_cfr // case-fatality rate
  param p_alpha // rate of burial

  param p_Inf
  param p_R0

  param p_early_H
  param p_late_H
  param p_H_alpha
  param p_H_tau

  param p_rep_d

  param p_vol_R0

  state S (has_output = 0)
  state E[rho_erlang] (has_output = 0)
  state Ic[gamma_erlang] (has_output = 0)
  state Ih[gamma_erlang,kappa_erlang] (has_output = 0)
  state Bc (has_output = 0)
  state Zc (has_output = 0)
  state Zh
  state Zd
  state R0
  state H
  state next_obs (has_output = 0)

  obs Admissions
  obs Deaths

  noise n_R0_walk
  noise n_admission

  sub transition {

    n_R0_walk ~ wiener()
    n_admission ~ gamma(shape = 100, scale = 0.01)

    // R0 trajectory: filled by R script
    // exponential: R0 <- R0 * exp(p_vol_R0 * n_R0_walk)
    // bounded: R0 <- max(0, R0 + p_vol_R0 * n_R0_walk)
    // independent: R0 <- max(0, p_vol_R0 * n_R0_walk)

    H <- p_early_H + (p_late_H * (1 - p_early_H)) / (1 + exp(p_H_alpha * (p_H_tau - t_now)))

    Zc <- (t_now <= first_obs ? 0 : (t_next_obs > next_obs ? 0 : Zc))
    Zh <- (t_now <= first_obs ? 0 : (t_next_obs > next_obs ? 0 : Zh))
    Zd <- (t_now <= first_obs ? 0 : (t_next_obs > next_obs ? 0 : Zd))
    next_obs <- (t_next_obs > next_obs ? t_next_obs : next_obs)

    inline kappa = 1 / admission_delay * rate_multiplier
    inline beta = R0 * p_gamma * p_alpha / (p_alpha + p_cfr * (1 - H) * p_gamma)
    ode {

      dS/dt =
      - beta * (Ic[0] + Ih[0,0] + Ih[0,1] + Ic[1] + Ih[1,0] + Ih[1,1] + Ic[2] + Ih[2,0] + Ih[2,1] + Bc) / p_N * S

      dE[rho_erlang]/dt =
      + (rho_erlang == 0 ? beta * (Ic[0] + Ih[0,0] + Ih[0,1] + Ic[1] + Ih[1,0] + Ih[1,1] + Ic[2] + Ih[2,0] + Ih[2,1] + Bc) / p_N * S : 0)
      + (rho_erlang > 0 ? e_rho * p_rho * E[rho_erlang - 1] : 0)
      - e_rho * p_rho * E[rho_erlang]

      dIc[gamma_erlang]/dt =
      + (gamma_erlang == 0 ? (1 - H) * e_rho * p_rho * E[e_rho - 1] : 0)
      + (gamma_erlang > 0 ? e_gamma * p_gamma * Ic[gamma_erlang - 1] : 0)
      - e_gamma * p_gamma * Ic[gamma_erlang]

      dIh[gamma_erlang,kappa_erlang]/dt =
      + (gamma_erlang + kappa_erlang == 0 ? H * e_rho * p_rho * E[e_rho - 1] : 0) // incubating proceeding to infectious and healthcare-seeking
      + (gamma_erlang > 0 ? e_gamma * p_gamma * Ih[gamma_erlang - 1,kappa_erlang] : 0) // proceeding through gamma stages
      - e_gamma * p_gamma * Ih[gamma_erlang,kappa_erlang] // proceeding through gamma stages
      + (kappa_erlang > 0 ? e_kappa * kappa * Ih[gamma_erlang, kappa_erlang - 1] : 0) // proceeding through kappa stages
      - (kappa_erlang < (e_kappa - 1) ? e_kappa * kappa * Ih[gamma_erlang,kappa_erlang] : e_kappa * kappa * Ih[gamma_erlang,kappa_erlang] * n_admission) // proceeding through kappa stages 

      dBc/dt  =
      + p_cfr * e_gamma * p_gamma * Ic[e_gamma - 1]
      - p_alpha * Bc

      dZc/dt =
      + e_rho * p_rho * E[e_rho - 1]

      dZh/dt =
      + e_kappa * kappa * Ih[0,e_kappa - 1] * n_admission
      + e_kappa * kappa * Ih[1,e_kappa - 1] * n_admission
      + e_kappa * kappa * Ih[2,e_kappa - 1] * n_admission

      dZd/dt =
      + e_gamma * p_gamma * (Ih[e_gamma - 1,0] + Ih[e_gamma - 1,1])
    }

  }

  sub parameter {
    p_N <- 270114
    p_rho <- 1 / (e_rho * 4.7) * rate_multiplier
    p_alpha <- 1 * rate_multiplier
    p_cfr <- 0.6695464
    p_gamma <- 1 / (e_gamma * 2.601496) * rate_multiplier
    p_Inf ~ uniform(0, 100)
    p_R0 ~ uniform(0, 10)
    p_vol_R0 ~ uniform(0, 1)
    p_early_H ~ uniform(0, 1)
    p_late_H ~ uniform(0, 1)
    p_H_tau ~ uniform(0, 21)
    p_H_alpha ~ uniform(0, 5)
    p_rep_d ~ uniform(0, 1)
  }

  sub initial {
    E[rho_erlang] <- p_Inf * p_gamma / (e_rho * (p_gamma + p_rho))
    H <- p_early_H + (p_late_H - p_early_H) / (1 + exp(p_H_alpha * p_H_tau))
    Ic[gamma_erlang] <- p_Inf * (1 - H) * p_rho / (e_gamma * (p_gamma + p_rho))
    Ih[gamma_erlang,kappa_erlang] <- p_Inf * H * p_rho / ((e_kappa * e_gamma) * (p_gamma + p_rho))
    S <- p_N - E[0] - E[1] - Ic[0] - Ic[1] - Ic[2] - Ih[0,0] - Ih[0,1] - Ih[1,0] - Ih[1,1] - Ih[2,0] - Ih[2,1]
    Bc <- 0
    Zc <- 0
    Zd <- 0
    Zh <- 0
    R0 <- p_R0
    next_obs <- 0
  }

  sub observation {
    Admissions ~ truncated_gaussian(Zh, max(sqrt(Zh), 1), lower = 0)
    Deaths ~ truncated_gaussian(p_rep_d * Zd, max(sqrt(p_rep_d * (1 - p_rep_d) * Zd), 1), lower = 0)
  }
}
