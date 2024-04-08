options(mc.cores = 4)
library("rsemicompstan")
library("ggplot2")
library("cowplot")
requireNamespace("directlabels")
requireNamespace("RColorBrewer")
requireNamespace("here")
requireNamespace("cobalt")
requireNamespace("MatchIt")

USE_FAKE_DATA <- FALSE #run on synthetic data
QUICK_RUN <- FALSE #short MCMC chainSs for quick code testing
RUN_DISCREPANCY <- FALSE #calculate discrepancy metrics

# Directory setup ---------------------------------------------------------

here::i_am("scripts/01-analyze_data_application.R")
data_dir <- here::here("data")
intermediates_dir <- here::here("intermediates")
output_dir <- here::here("output")
if (!dir.exists(intermediates_dir)) dir.create(intermediates_dir, showWarnings = FALSE)
if (!dir.exists(output_dir)) dir.create(output_dir, showWarnings = FALSE)


# Read data ---------------------------------------------------------------

data_filename <- ifelse(
  USE_FAKE_DATA,
  here::here("data", "scrData-fake.rds"),
  here::here("data", "scrData.rds")
)


# Quick run setup ---------------------------------------------------------

# Function factory allowing file name tagging of quick runs
fname_splicer <- function(suffix) {
  force(suffix)
  function(filename) {
    force(filename)
    fname_pieces <- strsplit(filename, split = ".", fixed = TRUE)[[1]]
    stopifnot(length(fname_pieces) == 2)
    return(paste0(fname_pieces[1], suffix, ".", fname_pieces[2]))  
  }
}

if (QUICK_RUN) {
  mcmc_params <- list(iter = 200, warmup = 100, chains = 2)
  qr_suffix <- "_QR"
  n_pp_draws <- 10
} else {
  mcmc_params <- list(iter = 4000, warmup = 3000, chains = 4)
  qr_suffix <- ""
  n_pp_draws <- 10000
}

qrify <- fname_splicer(suffix = qr_suffix)


# Model fits --------------------------------------------------------------

da_result <- rsemicompstan::fit_data_app(
  file = data_filename,
  use_priors = TRUE,
  ps_use = "match",
  init = "random",
  init_r = 0.3,
  shared_beta = 0,
  sigma_pa = 21,
  sigma_pb = 7.1,
  iter = mcmc_params$iter,
  warmup = mcmc_params$warmup,
  chains = mcmc_params$chains,
  seed = 2019,
  mc.cores = max(mcmc_params$chains, 4))

system.time({
  saveRDS(da_result, here::here("intermediates", qrify("da_result.rds")))
})


# (Optional) MCMC exploration ---------------------------------------------

# requireNamespace(shinystan)
# shinystan::launch_shinystan(da_result$stan_fit)


# Posterior predictions ---------------------------------------------------

set.seed(42)
system.time({
  fit_iter <- length(unlist(rstan::extract(da_result$stan_fit, par = "lp__")))
  iters_to_keep <- round(seq(1, fit_iter, length.out = min(1000, fit_iter)))
  pp <- posterior_predict_sample(
    da_result$stan_fit,
    yr = da_result$dat$yr,
    yt = da_result$dat$yt,
    dyr = da_result$dat$dyr,
    dyt = da_result$dat$dyt,
    z = da_result$dat$z,
    xmat = da_result$xmat,
    frailty_type = "impute")[, , iters_to_keep]
})

system.time({
  saveRDS(pp, here::here("intermediates", qrify("pp.rds")))
})

system.time({
  pp <- readRDS(here::here("intermediates", qrify("pp.rds")))
})


# Cohort data -------------------------------------------------------------

set.seed(1513)
system.time({
  gdat <- prepare_graph_data(pp, max_t = 90, length_out = 91)
})

system.time({
  saveRDS(gdat, here::here("intermediates", qrify("gdat.rds")))
})

system.time({
  gdat <- readRDS(here::here("intermediates", qrify("gdat.rds")))
})

set.seed(3700)
system.time({
  cgdat <- prepare_cohort_graph_data(
    pp,
    cohort = c(15, 30, 45, 60, 90),
    by = 1)
})

system.time({
  saveRDS(cgdat, here::here("intermediates", qrify("cgdat.rds")))
})

system.time({
  cgdat <- readRDS(here::here("intermediates", qrify("cgdat.rds")))
})



# Cohort descriptions -----------------------------------------------------

# 69% always-alive at 30 days
# 41% always-alive at 90 days
# mean(gdat$frac_aa[round(gdat$eval_t) == 30])
# mean(gdat$frac_aa[gdat$eval_t == 90])


# Single plot versions ----------------------------------------------------

pdf(here::here("intermediates", qrify("da_cohort_tvsace.pdf")), width = 7, height = 5)
  make_cohort_tvsace_plot(cgdat, time_unit = "Day") + 
  scale_x_continuous(breaks = seq(0, 90, 15), limits = c(0, 95)) +
  theme(legend.direction = "vertical",
        legend.position = "right",
        legend.text.align = 0.5) +
  guides(color = guide_colourbar(
    title = "Percent of \npopulation \nalways-alive at t",
    title.hjust = 0.5,
    label.position = "left")) +
  scale_color_gradientn(
    colors = RColorBrewer::brewer.pal(9, "YlOrRd"),
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = c("0%", "50%", "100%")) + 
  ggtitle("Difference in cumulative incidence of readmission attributable to extra care,\nby always-alive cohort")
dev.off()

pdf(here::here("intermediates", qrify("da_tvsace_plot.pdf")), width = 7, height = 5)
  make_tvsace_plot(
    plot_dat = gdat, 
    legend_position = "right",
    time_unit = "Day") + 
  scale_x_continuous(breaks = seq(0, 90, 15), limits = c(0, 95)) +
  ggtitle("Treatment-attributable difference in readmission incidence \naccumulated by Day t among always-survivors at t")
dev.off()

pdf(here::here("intermediates", qrify("da_cohort_rmsace.pdf")), width = 7, height = 5)
  make_cohort_rmsace_plot(cgdat, time_unit = "Day") + 
    scale_x_continuous(breaks = seq(0, 90, 15), limits = c(0, 95)) +
    theme(legend.direction = "vertical",
          legend.position = "right",
          legend.text.align = 0.5) +
    guides(color = guide_colourbar(title = "Percent of \npopulation \nalways-alive at t",
                                   title.hjust = 0.5,
                                   label.position = "left")) +
    scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, "YlOrRd"),
                          limits = c(0, 1),
                          breaks = c(0, 0.5, 1),
                          labels = c("0%", "50%", "100%")) + 
    ggtitle("Accumulated readmission-free days attributable to extra care,\nby always-alive cohort")
dev.off()

pdf(here::here("intermediates", qrify("da_rmsace_plot.pdf")), width = 7, height = 5)
  make_rmsace_plot(plot_dat = gdat, legend_position = "right",
                   time_unit = "Day") +
  scale_x_continuous(breaks = seq(0, 90, 15)) +
  ggtitle("Treatment-attributable readmission-free days \naccumulated by Day t among always-survivors at t")
dev.off()


pdf(here::here("intermediates", qrify("da_aa_kmplot.pdf")), width = 6, height = 4)
  make_aa_kmplot(plot_dat = gdat, time_unit = "Day") + 
  scale_x_continuous(breaks = seq(0, 90, 15)) +
  ggtitle("Population proportion in always-alive state")
dev.off()

pdf(here::here("intermediates", qrify("da_t_kmplot.pdf")), width = 6, height = 4)
  make_pp_z_kmplot(plot_dat = gdat, frac_var = "frac_a_t") +
  scale_x_continuous(breaks = seq(0, 90, 15)) +
  ggtitle("Survival curves for entire population under extra care (z=1)")
dev.off()

pdf(here::here("intermediates", qrify("da_c_kmplot.pdf")), width = 6, height = 4)
  make_pp_z_kmplot(plot_dat = gdat, frac_var = "frac_a_c") +
  scale_x_continuous(breaks = seq(0, 90, 15)) +
  ggtitle("Survival curves for entire population under no extra care (z=0)")
dev.off()

pdf(here::here("intermediates", qrify("da_kmplot_combined.pdf")), width = 7, height = 5)
  make_kmplot_combined(plot_dat = gdat, legend_position = "right") +
    scale_color_manual("Survival Type",
                       values = c("#2E655D", "#0B353B", "#EC3F19"),
                       labels = c("Always-alive", 
                                  "Under z = 0", 
                                  "Under z = 1")) +
    scale_linetype_manual("Survival Type", 
                          values = c(1, 3, 2),
                          labels = c("Always-alive", 
                                     "Under z = 0", 
                                     "Under z = 1")) +
  scale_x_continuous(breaks = seq(0, 90, 15)) + 
  ggtitle("Posterior survival curves with (z=1) and without (z=0) extra care,\nand corresponding fraction of population in always-alive state")
dev.off()

pdf(here::here("intermediates", qrify("da_comp_plot.pdf")), width = 7, height = 5)
  make_state_composition_plot(pp, maxt = 90, length_out = 30) +
    scale_x_continuous(breaks = seq(0, 90, 15))
dev.off()

if (RUN_DISCREPANCY) {
  discrep_plotdat <- make_da_discrepancy_plot(
    stan_res = da_result, seed = 5555,
    eval_t = c(15, 30, 45, 60, 90),
    subsamp = 4000)
  saveRDS(discrep_plotdat, here::here("intermediates", qrify("da_discrep_plotdat.rds")))
  discrep_plotdat <- readRDS(here::here("intermediates", qrify("da_discrep_plotdat.rds")))
  pdf(here::here("intermediates", qrify("da_discrep_plot.pdf")), width = 7, height = 5)
    print(discrep_plotdat)
  dev.off()

  pdf(here::here("intermediates", qrify("da_frailty_densities.pdf")), width = 6, height = 4)
    make_da_frailty_density(stan_fit = da_result$stan_fit, pp) + 
    ggtitle("Theoretical posterior predictive density of frailties vs.\ndensity of posterior predictions for observed sample") + 
    theme(legend.position = "bottom")
  dev.off()
  
  calculate_km_disc(
    eval_t = 15, 
    stan_res = da_result, 
    cens_times = rep(90, NROW(da_result$dat)), 
    seed = 123, subsamp = FALSE
  )
}


# Combined plots for paper ------------------------------------------------

# TV-SACE(r,t) (cohort) + TV-SACE(t,t) not-cohort
g_tv_p <- 
  make_tvsace_plot(
    plot_dat = gdat, 
    legend_position = "right",
    time_unit = "Day") + 
  scale_x_continuous(breaks = seq(0, 90, 15), limits = c(0, 95)) + 
  guides(color = "none") + 
  ylab(expression(`TV-SACE`[snap](t)))

snap_tv_scale <- ggplot_build(g_tv_p)$layout$panel_scales_y[[1]]$range$range

c_tv_p <- 
  make_cohort_tvsace_plot(cgdat, time_unit = "Day") + 
  scale_x_continuous(breaks = seq(0, 90, 15), limits = c(0, 95)) +
  directlabels::geom_dl(aes(label = cohort_id), 
                        method = list(directlabels::dl.trans(x = x + 0.25),
                                      "last.points", 
                                      cex = 1)) + 
  scale_x_continuous(breaks = seq(0, 90, 15), limits = c(0, 95)) +
  scale_y_continuous(limits = snap_tv_scale) +
  theme(legend.direction = "vertical",
        legend.position = "right",
        legend.text.align = 0.5) +
  guides(color = guide_colourbar(title = "Percent of \npopulation \nalways-alive at t",
                                 title.hjust = 0.5,
                                 label.position = "left")) +
  scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, "YlOrRd"),
                        limits = c(0, 1),
                        breaks = c(0, 0.5, 1),
                        labels = c("0%", "50%", "100%"))
cbar <- get_legend(c_tv_p)
c_tv_p <- c_tv_p + guides(color = "none")
tv_pair <- plot_grid(c_tv_p, g_tv_p,
                     # align = "vh",
                     labels = c("A", "B"),
                     hjust = -0.3,
                     nrow = 1)
tv_pair_with_legend <- plot_grid(tv_pair, cbar, 
                                 rel_widths = c(2, 0.25))



pdf(here::here("intermediates", qrify("da_tvsace_pair.pdf")), width = 9, height = 4)
  print(tv_pair_with_legend)
dev.off()


#Same pair but for RM-SACE
g_rm_p <- 
  make_rmsace_plot(plot_dat = gdat, legend_position = "right",
                   time_unit = "Day") + 
  scale_x_continuous(breaks = seq(0, 90, 15), limits = c(0, 95)) +
  guides(color = "none") + 
  ylab(expression(`RM-SACE`[snap](t)))
snap_rm_scale <- ggplot_build(g_rm_p)$layout$panel_scales_y[[1]]$range$range
c_rm_p <- 
  make_cohort_rmsace_plot(cgdat, time_unit = "Day") +
  directlabels::geom_dl(aes(label = cohort_id), 
                        method = list(directlabels::dl.trans(x = x + 0.25),
                                      "last.points", 
                                      cex = 1)) + 
  scale_x_continuous(breaks = seq(0, 90, 15), limits = c(0, 95)) +
  scale_y_continuous(limits = snap_rm_scale) +
  theme(legend.direction = "vertical",
        legend.position = "right",
        legend.text.align = 0.5) +
  guides(color = guide_colourbar(title = "Percent of \npopulation \nalways-alive at t",
                                 title.hjust = 0.5,
                                 label.position = "left")) +
  scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, "YlOrRd"),
                        limits = c(0, 1),
                        breaks = c(0, 0.5, 1),
                        labels = c("0%", "50%", "100%"))
cbar <- get_legend(c_rm_p)
c_rm_p <- c_rm_p + guides(color = FALSE)
rm_pair <- plot_grid(c_rm_p, 
                     g_rm_p,
                     # align = "vh",
                     labels = c("A", "B"),
                     hjust = -0.3,
                     nrow = 1)
rm_pair_with_legend <- plot_grid(rm_pair, cbar, 
                                 rel_widths = c(2, 0.25))

pdf(here::here("intermediates", qrify("da_rmsace_pair.pdf")), width = 9, height = 4)
  print(rm_pair_with_legend)
dev.off()


# Cohort pairs
cohort_pair <- plot_grid(c_tv_p + ggtitle("TV-SACE"), 
                           c_rm_p + ggtitle("RM-SACE"),
                     hjust = -0.3,
                     nrow = 1)
cohort_pair_with_legend <- plot_grid(cohort_pair, cbar, 
                                 rel_widths = c(2, 0.25))
pdf(here::here("intermediates", qrify("da_cohort_pair.pdf")), width = 9, height = 4)
  print(cohort_pair_with_legend)
dev.off()


# All 4 graphs together
rm_pair2 <- plot_grid(c_rm_p, g_rm_p,
                      labels = c("C", "D"),
                      hjust = -0.3,
                      nrow = 1)
quad <- plot_grid(tv_pair, rm_pair2, nrow = 2)
quad_with_legend <- plot_grid(quad, cbar, 
                              rel_widths = c(2, 0.25))

pdf(here::here("output", qrify("da_sace_quad.pdf")), width = 11, height = 9)
  print(quad_with_legend)
dev.off()


# All-3 survival and principal state composition pair
km_all <- make_kmplot_combined(plot_dat = gdat, legend_position = "bottom",
                               time_unit = "Day") +
  scale_x_continuous(breaks = seq(0, 90, 15)) + 
  ggtitle("Posterior survival curves with (z=1) and without (z=0) extra care,\nand corresponding fraction of population in always-alive state")
comp <- make_state_composition_plot(pp, maxt = 90, length_out = 30) +
  scale_x_continuous(breaks = seq(0, 90, 15)) +
  scale_fill_manual(values = color_vals, labels = rev(state_names)) +
  guides(fill = guide_legend(reverse = TRUE)) 
km_comp_pair <- plot_grid(km_all, comp, nrow = 1, align = "vh", hjust = -0.3, labels = c("A", "B"))

pdf(here::here("output", qrify("da_kmplot_comp_pair_UPDATED.pdf")), width = 12, height = 6)
  print(km_comp_pair)
dev.off()


# Posterior prediction table ----------------------------------------------

system.time({
  set.seed(42787)
  ppred_tab <- ppred_table(stan_res = da_result, eval_t = c(30, 90),
                         frailty_q = c(0.9, 0, 0.1),
                         thin_out = 0, nrep = 10000)
                         # thin_out = 0, nrep = 10)
saveRDS(ppred_tab, here::here("intermediates", qrify("ppred_tab.rds")))

ppred_tab <- readRDS(here::here("intermediates", qrify("ppred_tab.rds")))

writeLines(
  make_ppred_table_pretty(tab = ppred_tab, digits = 3, caption = NULL),
  here::here("output", qrify("da_ppred_table.tex")))
})


# Checking paper numbers (in order) ---------------------------------------

rstan::summary(da_result$stan_fit, par = "sigma")$summary
p_est_sig <- rstan::summary(da_result$stan_fit, par = "sigma")$summary[,"mean"]
(haz_ratio <- qgamma(p = 0.9, 1 / p_est_sig, 1 / p_est_sig) / 
              qgamma(p = 0.1, 1 / p_est_sig, 1 / p_est_sig))


# Balance plots -----------------------------------------------------------

recreate_matched <- function(data_filename) {
  # Read in data
  scr_pc <- readRDS(file = data_filename)
  
  # Alias variables
  scr_pc$yr <- scr_pc$time1
  scr_pc$yt <- scr_pc$time2
  scr_pc$dyr <- scr_pc$event1
  scr_pc$dyt <- scr_pc$event2
  
  # Fuzz so that events never happen at time zero
  # Add one-half day to sojourn times of zero
  not_zero_soj   <- which(!((scr_pc$yt == scr_pc$yr) & (scr_pc$dyr == 1)))
  zero_soj <- (1:NROW(scr_pc))[-not_zero_soj]
  scr_pc$yt[zero_soj] <- scr_pc$time2[zero_soj] + 0.5
  
  # Exclude anyone discharged to something besides home or home with care
  scr_pc <- scr_pc[pmax(scr_pc$disc_snf_icf,
                        scr_pc$disc_hospice,
                        scr_pc$disc_other) == 0, ]
  
  # Adjustment covariates are race, standardized age, sex, comorbidity score, 
  # admission route, and hospital length of (initial) stay
  # "Treatment" is being discharged to home with care (vs. home with no care)
  # This is reversal of home vs. not-home from before!
  scr_pc$z <- scr_pc$disc_homecare
  
  matches <- MatchIt::matchit(
    z ~ race_noWhite + sex_female + deyo2_ + adm + 
      age_std + los_std + age_std:age_std + 
      age_std^2 + los_std^2, 
    data = scr_pc,
    method = "nearest", 
    discard = "control",
    ratio = 1,
    replace = FALSE)
  return(matches)
}

bal_tab <- cobalt::bal.tab(recreate_matched(data_filename = data_filename))
pretty_var_names <- data.frame(
  old = c("distance", "los_std", "age_std", "sex_female", "deyo2_", "race_noWhite", "adm"),
  new = c("Distance", "Hospital length of stay", "Age", "Sex", "Comorbidity score 2+", 
          "Non-white race", "Route of admission")
)

# Covariate (im)balance plot pre and post matching
pdf(here::here("output", qrify("da_postmatch_covbalance.pdf")), width = 5, height = 5)
  cobalt::love.plot(bal_tab, var.names = pretty_var_names,
            title = "Covariate balance",
            subtitle = "After matching on propensity scores",
            var.order = "unadjusted",
            abs = FALSE,
            shapes = c("triangle", "circle"),
            colors = 1:2,
            sample.names = c("Pre-matching", "Post-matching")) +
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "bottom")
dev.off()

