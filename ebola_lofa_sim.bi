model ebola_lofa_sim {

  const first_obs = 0
  const rate_multiplier = 7
 
  const e_rho = 2
  const e_gamma = 3
  const e_delta = 2
  const e_kappa = 2
  const e_theta = 10

  state admission_delay
  state K
  state late_increase

  dim rho_erlang(e_rho)
  dim gamma_erlang(e_gamma)
  dim kappa_erlang(e_kappa)
  dim theta_erlang(e_theta)
  dim delta_erlang(e_theta)

  param p_N // population size
  param p_rho // rate of developing symptoms
  param p_gamma // recovery/death rate in the community
  param p_cfr // case-fatality rate
  param p_alpha // rate of burial
  param p_theta // rate of release from hospital after recovery
  param p_delta // rate of death in hospital

  param p_Inf

  param p_rep_d

  state S
  state E[rho_erlang]
  state Ic[gamma_erlang]
  state Ih[gamma_erlang,kappa_erlang]
  state Hr[theta_erlang]
  state Hd[delta_erlang]
  state Bc
  state Zc
  state Zh
  state R0
  state H

  obs Admissions

  state n_admission

  sub transition {

    Zc <- 0
    Zh <- 0

    inline kappa = 1 / admission_delay * rate_multiplier
    inline beta = R0 * p_gamma * p_alpha / (p_alpha + p_cfr * (1 - H) * p_gamma)
    inline hospital_open = Hr[0] + Hr[1] + Hr[2] + Hr[3] + Hr[4] + Hr[5] + Hr[6] + Hr[7] + Hr[8] + Hr[9] + Hd[0] + Hd[1] < K ? 1 : 0

    ode {

      dS/dt =
      - beta * (Ic[0] + Ih[0,0] + Ih[0,1] + Ic[1] + Ih[1,0] + Ih[1,1] + late_increase * (Ic[2] + Ih[2,0] + Ih[2,1] + Bc)) / p_N * S

      dE[rho_erlang]/dt =
      + (rho_erlang == 0 ? beta * (Ic[0] + Ih[0,0] + Ih[0,1] + Ic[1] + Ih[1,0] + Ih[1,1] + late_increase * (Ic[2] + Ih[2,0] + Ih[2,1] + Bc)) / p_N * S : 0)
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
      - (kappa_erlang < (e_kappa - 1) ? (e_kappa * kappa * Ih[gamma_erlang,kappa_erlang]) : e_kappa * kappa * Ih[gamma_erlang,kappa_erlang] * n_admission * hospital_open) // proceeding through kappa stages 

      dHr[theta_erlang]/dt =
      + (theta_erlang == 0 ? ((1 - p_cfr) * e_kappa * kappa * (Ih[0,e_kappa - 1] + Ih[1,e_kappa - 1] + Ih[2,e_kappa - 1]) * n_admission * hospital_open) : 0)
      + (theta_erlang > 0 ? e_theta * p_theta * Hr[theta_erlang - 1] : 0)
      - e_theta * p_theta * Hr[theta_erlang]

      dHd[delta_erlang]/dt =
      + (delta_erlang == 0 ? p_cfr * e_kappa * kappa * (Ih[0,e_kappa - 1] + Ih[1,e_kappa - 1] + Ih[2,e_kappa - 1]) * n_admission * hospital_open : 0)
      + (delta_erlang > 0 ? e_delta * p_delta * Hd[delta_erlang - 1] : 0)
      - e_delta * p_delta * Hd[delta_erlang]

      dBc/dt  =
      + p_cfr * e_gamma * p_gamma * Ic[e_gamma - 1]
      - p_alpha * Bc

      dZc/dt =
      + e_rho * p_rho * E[e_rho - 1]

      dZh/dt =
      + e_kappa * kappa * (Ih[0,e_kappa - 1] + Ih[1,e_kappa - 1] + Ih[2,e_kappa - 1]) * n_admission * hospital_open
    }

  }

  sub parameter {
    p_N <- 270114
    p_rho <- 1 / (e_rho * 4.7) * rate_multiplier
    p_alpha <- 1 * rate_multiplier
    p_cfr <- 0.6695464
    p_gamma <- 1 / (e_gamma * 2.601496) * rate_multiplier
    p_theta <- 1 / (e_theta * 1.279192) * rate_multiplier
    p_delta <- 1 / (e_delta * 2.302684) * rate_multiplier
    p_Inf ~ uniform(0, 100)
    p_rep_d ~ uniform(0, 1)
  }

  sub initial {
    E[rho_erlang] <- p_Inf * p_gamma / (e_rho * (p_gamma + p_rho))
    Ic[gamma_erlang] <- p_Inf * (1 - H) * p_rho / (e_gamma * (p_gamma + p_rho))
    Ih[gamma_erlang,kappa_erlang] <- p_Inf * H * p_rho / ((e_kappa * e_gamma) * (p_gamma + p_rho))
    S <- p_N - E[0] - E[1] - Ic[0] - Ic[1] - Ic[2] - Ih[0,0] - Ih[0,1] - Ih[1,0] - Ih[1,1] - Ih[2,0] - Ih[2,1]
    Hr[theta_erlang] <- 0
    Hd[delta_erlang] <- 0
    Bc <- 0
    Zc <- 0
    Zh <- 0
  }

  sub observation {
    Admissions ~ truncated_gaussian(Zh, max(sqrt(Zh), 1), lower = 0)
  }
}
