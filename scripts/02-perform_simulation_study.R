options(mc.cores = 4)
library("batchtools")
library("rsemicompstan")
library("ggplot2")
library("rstan")
library("cowplot")
library("dplyr")
requireNamespace("directlabels")
requireNamespace("RColorBrewer")
requireNamespace("reshape2")
requireNamespace("here")
requireNamespace("data.table")

QUICK_RUN <- TRUE #short MCMC chains and few replicates for quick code testing


# Directory setup ---------------------------------------------------------

here::i_am("scripts/02-perform_simulation_study.R")
inst_dir <- here::here("inst")
data_dir <- here::here("data")
intermediates_dir <- here::here("intermediates")
registry_dir <- here::here("intermediates", "registries")
if (!dir.exists(registry_dir)) dir.create(registry_dir, showWarnings = FALSE)
output_dir <- here::here("output")

batchtools_config <- here::here("inst", ".batchtools.conf.R")


# Quick run setup ---------------------------------------------------------

if (QUICK_RUN) {
  replicates_per_scenario <- 5
  iter <- 200
  warmup <- 100
  qr_suffix <- "_QR"
} else {
  replicates_per_scenario <- 500
  iter <- 5000
  warmup <- 3000
  qr_suffix <- ""
}

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
qrify <- fname_splicer(suffix = qr_suffix)


# Simulation study setup --------------------------------------------------



n <- 6000
chains <- 1

num_scenarios <- 3
time_each_in_min <- c(60 * 4, 60 * 4, 60 * 4)
memorys <- c(4000, 4000, 4000)
adapt_deltas <- c(0.85, 0.90, 0.90)
stepsizes <- c(2, 1.5, 1.5)
max_treedepths <- c(11, 11, 11)

regs <- as.list(rep(NA, num_scenarios))
for (s in 1:num_scenarios) {
  registry_name <- paste0("registry_s", s, qr_suffix)
  registry_loc <- paste0(registry_dir, "/", registry_name)
  if (dir.exists(registry_loc)) {
    regs[[s]] <- loadRegistry(
      file.dir = registry_loc,
      conf.file = batchtools_config,
      make.default = FALSE,
      writeable = TRUE)
  } else {
    regs[[s]] <- makeRegistry(
      file.dir = registry_loc,
      conf.file = batchtools_config,
      make.default = FALSE)
  }
}


for (s in 1:num_scenarios) {
  submit_scenario_jobs(
    registry = regs[[s]], 
    scenario = s,
    seed = (1:replicates_per_scenario) * 73,
    clear_existing = TRUE,
    n = n,
    init = get_init_truth(scenario = s, chains = chains, with_randomness = TRUE),
    shared_beta = 1,
    sigma_pa = 18/7, #increase denominator for vaguer prior on sigma
    sigma_pb = 25/7, #keep denominator synced with denominator for sigma_pa
    iter = iter,
    warmup = warmup,
    chains = chains,
    eval_t = c(30, 60, 90),
    parallelize_chains = FALSE,
    control = list(
      adapt_delta = adapt_deltas[s],
      stepsize = stepsizes[s],
      max_treedepth = max_treedepths[s]
      ),
    time_each = time_each_in_min[s],
    memory = memorys[s],
  submit = TRUE)
}


# Interactive job monitoring ----------------------------------------------

# submitJobs(findNotDone(reg = regs[[1]]), reg = regs[[1]])
# submitJobs(findNotDone(reg = regs[[2]]), reg = regs[[2]])
# submitJobs(findNotDone(reg = regs[[3]]), reg = regs[[3]])
# getJobStatus(reg = regs[[1]])
# findRunning(reg = regs[[1]])


# Frailty distribution check ----------------------------------------------

scenario_s <- 1

nR <- NROW(getJobTable(reg = regs[[scenario_s]]))
sig_true <- return_dgp_parameters(scenario_s)$sigma
sigma_replicate_n <- 2500
f_aa_ex2 <- f_aa_ex1 <- matrix(NA, nrow = sigma_replicate_n, ncol = nR)

for (result_i in 1:nR){
  res <- loadResult(id = result_i, reg = regs[[scenario_s]])
  my_const_f <- matrix(
    rgamma(
      n = NROW(res$result$xmat) * 1, 
      shape = 1 / sig_true, 
      rate = 1 / sig_true), 
    nrow = NROW(res$result$xmat), 
    ncol = sigma_replicate_n)
  ex1_pp <- posterior_predict_sample(
    stan_fit = res$result$stan_fit,
    yr = res$result$dat$yr,
    yt = res$result$dat$yt,
    dyr = res$result$dat$dyr,
    dyt = res$result$dat$dyt,
    z = res$result$dat$z,
    xmat = res$result$xmat,
    frailty_type = "given",
    frailty = my_const_f)
  ex2_pp <- pp_from_result_list(rl = res$result)
  f_aa_ex1[, result_i] <- apply(ex1_pp, MARGIN = 3, calculate_frac_aa, eval_t = 30)
  f_aa_ex2[, result_i] <- apply(ex2_pp, MARGIN = 3, calculate_frac_aa, eval_t = 30)
  
  # Progress check for slow runs
  if (!QUICK_RUN && (result_i %% 5 == 0)) {
    cat(paste0("\n", result_i, "!"))
  }
}


# calc_frac_comp <- function(eval_t = 30, pp) {
#   pp <- as.data.frame(pp)
#   pstates <- make_pstates(eval_t, pp)
#   stopifnot(anyNA(pstates) == FALSE)
#   frac_aa <- mean(pstates == "AA")
#   frac_ck <- mean(pstates == "CK")
#   frac_tk <- mean(pstates == "TK")
#   frac_dd <- mean(pstates == "DD")
#   return(c("AA" = frac_aa, "CK" = frac_ck, "TK" = frac_tk, "DD" = frac_dd))
# }

f_aa_ests <- rbind(data.frame(Estimate = c(f_aa_ex1), NumType = 0),
                   data.frame(Estimate = c(f_aa_ex2), NumType = 1))
f_aa_ests$Type <- factor(f_aa_ests$NumType,
                            levels = c(0:1),
                            labels = c("Sigma From Posterior",
                                       "Sigma Fixed At Truth"))

pdf(here::here("intermediates", qrify(paste0("frac_aa_scen", scenario_s, ".pdf")), width = 7, height = 5))
  ggplot(f_aa_ests, aes(x = Estimate, group = Type, color = Type, fill = Type)) +
    geom_density(aes(color = Type), alpha = 0.5) +
    geom_histogram(alpha = 0.5, binwidth = 0.0025, color = "black") +
    geom_vline(xintercept = ts[,"frac_aa"], size = 1.5) +
    labs(y = "Density", title = "Fraction always-alive estimates",
         subtitle = "From 100 replicates the correctly specified scenario,
  using posterior draws of other parameters to calculated the fraction always-alive,
  and with true always-alive fraction (vertical line) calculated in a N=10^7 data set") +
    theme(legend.position = "bottom",
          plot.subtitle = element_text(size = 8))
dev.off()


#' Function to prepare data for caterpillar plot using registry replicates
cat_prep <- function(s, eval_t = c(30, 60, 90), cens_time = 90, num_rep = 100, reg_list) {
  truths <- expand.grid(eval_t = eval_t, scenario = s, rep_num = 1:num_rep, 
                        tv_sace = NA, rm_sace = NA, frac_aa = NA, frac_aa2 = NA)

  res_df <- as.data.frame(matrix(NA, nrow = num_rep * length(eval_t), ncol = 13))
  colnames(res_df) <- c("scenario", "rep_num", "eval_t", 
                        "tv", "tv_ci_lb", "tv_ci_ub",
                        "rm", "rm_ci_lb", "rm_ci_ub",
                        "aa", "aa_ci_lb", "aa_ci_ub", "aa2")
  stopifnot(NROW(truths) == NROW(res_df))
  res_df$rep_num = truths$rep_num
  res_df$eval_t = truths$eval_t
  for (i in 1:num_rep) {
    ijp <- getJobPars(id = i, reg = reg_list[[s]])$job.pars[[1]]
    obs_false <- simulate_scenario(
      n = ijp$n, seed = ijp$seed, 
      scenario = ijp$scenario, 
      censor = TRUE, cens_times = rep(cens_time, n),
      observed = FALSE, add_imp = TRUE)
    truths$tv_sace[truths$rep_num == i] <- 
      v_calculate_tv_sace(eval_t = eval_t, pp = obs_false)
    truths$rm_sace[truths$rep_num == i] <-
      v_calculate_rm_sace(eval_t = eval_t, pp = obs_false)
    truths$frac_aa[truths$rep_num == i] <- 
      v_calculate_frac_aa(eval_t = eval_t, pp = obs_false)
    
    reg_res <- loadResult(i, reg = reg_list[[s]])
    tv_ci <- reg_res$oc$cis[c("tv_sace.2.5.", "tv_sace.97.5.")]
    rm_ci <- reg_res$oc$cis[c("rm_sace.2.5.", "rm_sace.97.5.")]
    aa_ci <- reg_res$oc$cis[c("frac_aa.2.5.", "frac_aa.97.5.")]
    
    res_df[res_df$rep_num == i, c("scenario")] <- ijp$scenario
    res_df[res_df$rep_num == i, c("tv_ci_lb", "tv_ci_ub")] <- tv_ci
    res_df[res_df$rep_num == i, c("rm_ci_lb", "rm_ci_ub")] <- rm_ci
    res_df[res_df$rep_num == i, c("aa_ci_lb", "aa_ci_ub")] <- aa_ci
    res_df[res_df$rep_num == i, c("tv")] <- truths$tv_sace[truths$rep_num == i]
    res_df[res_df$rep_num == i, c("rm")] <- truths$rm_sace[truths$rep_num == i]
    res_df[res_df$rep_num == i, c("aa")] <- truths$frac_aa[truths$rep_num == i]
  }
  return(res_df)
}

eval_t <- c(30, 60, 90)
wide <- cat_prep(s = 1, eval_t = eval_t, cens_time = 90, num_rep = 200, reg_list = regs)
superpop_all <- summarize_scenario_truths(scenarios = 1)[, c("eval_t", "frac_aa", "tv_sace", "rm_sace")]
superpop_aa <- superpop_all[, c("eval_t", "frac_aa")]
superpop_tv <- superpop_all[, c("eval_t", "tv_sace")]
superpop_rm <- superpop_all[, c("eval_t", "rm_sace")]
long <- reshape2::melt(wide, id.vars = c("scenario", "rep_num", "eval_t"))
tv_long <- long[long$variable %in% c("tv", "tv_ci_lb", "tv_ci_ub"), ]
rm_long <- long[long$variable %in% c("rm", "rm_ci_lb", "rm_ci_ub"), ]
aa_long <- long[long$variable %in% c("aa", "aa_ci_lb", "aa_ci_ub"), ]


# Always-alive plot -------------------------------------------------------

plot_t <- 30

aa_wide <- reshape2::dcast(
  data = aa_long, 
  formula = eval_t + rep_num + scenario ~ variable, 
  value.var = "value")

aa_30w <- aa_wide[aa_wide$eval_t == plot_t, ]
aa_30w <- aa_30w[order(aa_30w$aa), ]
aa_30w$sorted_rep <- 1:NROW(aa_30w)
aa_30w$no_cover <- as.factor(pmax((aa_30w$aa < aa_30w$aa_ci_lb),
                                  (aa_30w$aa > aa_30w$aa_ci_ub)))

pdf(here::here("intermediates", qrify(paste0("finite_sample_frac_aa_t", plot_t, ".pdf"))), width = 4, height = 8)
  cat_p1 <- 
    ggplot(data = aa_30w) + 
    geom_errorbar(mapping = aes(x = sorted_rep, ymin = aa_ci_lb, ymax = aa_ci_ub), width = 0) + 
    geom_point(mapping = aes(x = sorted_rep, y = aa, color = no_cover), size = 0.6, stroke = 0.6) + 
    geom_hline(yintercept = superpop_aa$frac_aa[superpop_aa == plot_t], color = "gray") +
    scale_color_manual("", values = c("#000000", "#FF0000"), guide = FALSE) +
    scale_x_reverse() +
    labs(x = "Replicates (sorted by underlying finite sample truth)", 
         y = paste0("Finite sample fraction always-alive at ", plot_t),
         title = "") +
    coord_flip() +
    theme(axis.ticks.y = element_blank(), 
          axis.text.y=element_blank())
  print(cat_p1)
dev.off()



# TV coverage -------------------------------------------------------

tv_wide <- reshape2::dcast(
  data = tv_long, 
  formula = eval_t + rep_num + scenario ~ variable, 
  value.var = "value")

tv_30w <- tv_wide[tv_wide$eval_t == plot_t, ]
tv_30w <- tv_30w[order(tv_30w$tv), ]
tv_30w$sorted_rep <- 1:NROW(tv_30w)
tv_30w$no_cover <- as.factor(pmax((tv_30w$tv < tv_30w$tv_ci_lb),
                                  (tv_30w$tv > tv_30w$tv_ci_ub)))
pdf(here::here("intermediates", paste0( "finite_sample_tv_sace_t", plot_t, qr_suffix, ".pdf")), width = 4, height = 8)
  cat_p2 <- ggplot(data = tv_30w) + 
    geom_errorbar(mapping = aes(x = sorted_rep, ymin = tv_ci_lb, ymax = tv_ci_ub), width = 0) + 
    geom_point(mapping = aes(x = sorted_rep, y = tv, color = no_cover), size = 0.6, stroke = 0.6) + 
    geom_hline(yintercept = superpop_tv$tv_sace[superpop_tv == plot_t], color = "gray") + 
    scale_color_manual("", values = c("#000000", "#FF0000"), guide = FALSE) +
    scale_x_reverse() +
    labs(x = "Replicates (sorted by underlying finite sample truth)", 
         y = paste0("Finite sample TV-SACE(", plot_t, ", ", plot_t, ")"),
         title = "") +
    coord_flip() +
    theme(axis.ticks.y = element_blank(), 
          axis.text.y = element_blank())
  print(cat_p2)
dev.off()


# RM coverage -------------------------------------------------------

rm_wide <- reshape2::dcast(
  data = rm_long, 
  formula = eval_t + rep_num + scenario ~ variable, 
  value.var = "value")

rm_30w <- rm_wide[rm_wide$eval_t == plot_t, ]
rm_30w <- rm_30w[order(rm_30w$rm), ]
rm_30w$sorted_rep <- 1:NROW(rm_30w)
rm_30w$no_cover <- as.factor(pmax((rm_30w$rm < rm_30w$rm_ci_lb),
                                  (rm_30w$rm > rm_30w$rm_ci_ub)))

pdf(here::here("intermediates", qrify(paste0( "finite_sample_rm_sace_t", plot_t, ".pdf"))), width = 4, height = 8)
  cat_p3 <- ggplot(data = rm_30w) + 
    geom_errorbar(mapping = aes(x = sorted_rep, ymin = rm_ci_lb, ymax = rm_ci_ub), width = 0) + 
    geom_point(mapping = aes(x = sorted_rep, y = rm, color = no_cover), size = 0.6, stroke = 0.6) + 
    geom_hline(yintercept = superpop_rm$rm_sace[superpop_rm == plot_t], color = "gray") + 
    scale_color_manual("", values = c("#000000", "#FF0000"), guide = FALSE) +
    scale_x_reverse() +
    labs(x = "Replicates (sorted by underlying finite sample truth)", 
         y = paste0("Finite sample RM-SACE(", plot_t, ", ", plot_t, ")"),
         title = "") +
    coord_flip() +
    theme(axis.ticks.y = element_blank(), 
          axis.text.y = element_blank())
  print(cat_p3)
dev.off()


# Operating characteristics 
oc_tab_raw <- as.data.frame(data.table::rbindlist(lapply(regs, FUN = get_registry_oc)))

saveRDS(oc_tab_raw, here::here("intermediates", qrify("oc_tab_raw.rds")))

oc_tab_raw <- readRDS(here::here("intermediates", qrify("oc_tab_raw.rds")))


# Triplet of all caterpillars together
cat_leg <- get_legend(cat_p3 + 
  scale_color_manual("Covers finite sample truth", 
                     labels = c("Yes", "No"),
                     values = c("black", "red")) +
  guides(color = guide_legend(override.aes = list(size = 1.8))) + 
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = "center"))


# Finite truth coverage probabilities
mean(aa_30w$no_cover == "0")
mean(tv_30w$no_cover == "0")
mean(rm_30w$no_cover == "0")

cov_fs <- c(mean(aa_30w$no_cover == "0"),
            mean(tv_30w$no_cover == "0"),
            mean(rm_30w$no_cover == "0")) * 100
cov_mat <- format(cbind(cov_fs, cov_p), digits = 3)
make_cat_lab <- function(plot_num, cov_mat) {
  paste0("95% interval coverage\n",
         "Finite sample: ", cov_mat[plot_num, 1], "%\n",
         "Population:      ", cov_mat[plot_num, 2], "%")
}
l1 <- make_cat_lab(plot_num = 1, cov_mat = cov_mat)
l2 <- make_cat_lab(plot_num = 2, cov_mat = cov_mat)
l3 <- make_cat_lab(plot_num = 3, cov_mat = cov_mat)

pdf(here::here("output", qrify(paste0("caterpillars", plot_t, ".pdf"))), width = 8, height = 8)
  cat3 <- plot_grid(
            cat_p1 + labs(caption = l1, y = "Fraction always-alive at 30") + 
              theme(plot.caption = element_text(hjust = 0.5, vjust = 0.5, 
                                                face = "italic")),
            cat_p2 + labs(caption = l2, y = "TV-SACE(30, 30)", x = "") + 
              theme(plot.caption = element_text(hjust = 0.5, vjust = 0.5,
                                                face = "italic")),
            cat_p3 + labs(caption = l3, y = "RM-SACE(30, 30)", x = "") +
              theme(plot.caption = element_text(hjust = 0.5, vjust = 0.5,
                                                face = "italic")),
            nrow = 1, 
            labels = c("A", "B", "C"),
            hjust = -0.3)
  plot_grid(cat3, cat_leg, nrow = 2, rel_heights = c(10, 0.5))
dev.off()


# Convergence checks ------------------------------------------------------


r1_conv <- sapply(1:replicates_per_scenario, FUN = function(x, reg = regs[[1]]) {
apply_simstudy_conv_criteria(loadResult(x, reg = reg))}
  ) 
any(r1_conv)

r2_conv <- sapply(1:replicates_per_scenario, FUN = function(x, reg = regs[[2]]) {
  apply_simstudy_conv_criteria(loadResult(x, reg = reg))}
) 
any(r2_conv)

r3_conv <- sapply(1:replicates_per_scenario, FUN = function(x, reg = regs[[3]]) {
  apply_simstudy_conv_criteria(loadResult(x, reg = reg))}
) 
any(r3_conv)



# Summary tables ----------------------------------------------------------

# Operating characteristics table
oc_tab <- make_oc_table(
  regs = regs,
  digits = 3,
  format = "latex",
  escape = FALSE,
  booktabs = TRUE)

writeLines(oc_tab, here::here("output", qrify( "sim_oc_table.tex")))

