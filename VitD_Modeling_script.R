# Vitamin D - Boxplots (Mean & Median) + Jitter Plot WITH p-values

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(svglite)
library(extrafont)

loadfonts(device = "all", quiet = TRUE)

data <- read_excel("VitDPAS_vitD_data.xlsx")

data_long <- data %>%
  dplyr::selectdplyr::select(Individual, d0, d1, d28, d29, d56, d57, d84) %>%
  pivot_longer(cols = starts_with("d"),
               names_to  = "day_label",
               values_to = "vitamin_d_level") %>%
  mutate(day = as.numeric(gsub("d", "", day_label))) %>%
  dplyr::selectdplyr::select(-day_label) %>%
  filter(!is.na(Individual))


day_mapping <- data.frame(
  day      = c(0, 1, 28, 29, 56, 57, 84),
  day_plot = c(0, 2.5, 12, 14.5, 24, 26.5, 36)
)

data_long_plot <- data_long %>% left_join(day_mapping, by = "day")

summary_stats <- data_long %>%
  group_by(day) %>%
  summarise(
    mean   = mean(vitamin_d_level, na.rm = TRUE),
    median = median(vitamin_d_level, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(day_mapping, by = "day")

# WILCOXON TESTS 

df_wide <- data %>% filter(!is.na(Individual))

run_wilcox <- function(col1, col2) {
  wilcox.test(df_wide[[col1]], df_wide[[col2]],
              paired = TRUE, exact = FALSE)
}

comparisons_list <- list(
  list(col1="d0",  col2="d1",  label="d0 vs d1"),
  list(col1="d28", col2="d29", label="d28 vs d29"),
  list(col1="d56", col2="d57", label="d56 vs d57"),
  list(col1="d0",  col2="d28", label="d0 vs d28"),
  list(col1="d0",  col2="d56", label="d0 vs d56"),
  list(col1="d0",  col2="d84", label="d0 vs d84")
)

tests  <- lapply(comparisons_list, function(x) run_wilcox(x$col1, x$col2))
names(tests) <- sapply(comparisons_list, function(x) x$label)
pvals <- sapply(tests, function(t) t$p.value)

fmt_p <- function(p) {
  if (p < 0.0001) "p < 0.0001"
  else paste0("p = ", formatC(p, format="f", digits=4))
}


cat("PAIRED WILCOXON SIGNED-RANK TEST RESULTS\n")
for (nm in names(pvals)) {
  cat(paste0("  ", nm, ": ", fmt_p(pvals[nm]), "\n"))
}

# BRACKET FUNCTIONS


BTEXT <- 8

sig_bracket <- function(x1, x2, y_top, label) {
  list(
    annotate("segment", x=x1,xend=x2,y=y_top,yend=y_top,linewidth=1.2),
    annotate("segment", x=x1,xend=x1,y=y_top,yend=y_top-1.2,linewidth=1.2),
    annotate("segment", x=x2,xend=x2,y=y_top,yend=y_top-1.2,linewidth=1.2),
    annotate("text", x=(x1+x2)/2,y=y_top+1,label=label,
             size=BTEXT,fontface="bold",family="Arial")
  )
}

sig_bracket_below <- function(x1, x2, y_bot, label) {
  list(
    annotate("segment", x=x1,xend=x2,y=y_bot,yend=y_bot,linewidth=1.2),
    annotate("segment", x=x1,xend=x1,y=y_bot,yend=y_bot+1.2,linewidth=1.2),
    annotate("segment", x=x2,xend=x2,y=y_bot,yend=y_bot+1.2,linewidth=1.2),
    annotate("text", x=(x1+x2)/2,y=y_bot-1,label=label,
             size=BTEXT,fontface="bold",family="Arial")
  )
}

# BRACKET POSITIONS


y_data_max <- max(data_long$vitamin_d_level, na.rm=TRUE)
y_data_min <- min(data_long$vitamin_d_level, na.rm=TRUE)

b_bolus_y  <- y_data_min - 3
b_cross1_y <- y_data_max + 4
b_cross2_y <- y_data_max + 8
b_cross3_y <- y_data_max + 12

all_brackets <- c(
  sig_bracket_below(0,2.5,b_bolus_y,fmt_p(pvals["d0 vs d1"])),
  sig_bracket_below(12,14.5,b_bolus_y,fmt_p(pvals["d28 vs d29"])),
  sig_bracket_below(24,26.5,b_bolus_y,fmt_p(pvals["d56 vs d57"])),
  sig_bracket(0,12,b_cross1_y,fmt_p(pvals["d0 vs d28"])),
  sig_bracket(0,24,b_cross2_y,fmt_p(pvals["d0 vs d56"])),
  sig_bracket(0,36,b_cross3_y,fmt_p(pvals["d0 vs d84"]))
)

y_axis_max <- y_data_max + 14
y_axis_min <- y_data_min - 7

x_scale <- scale_x_continuous(
  breaks=c(0,2.5,12,14.5,24,26.5,36),
  labels=c("0","1","28","29","56","57","84")
)

y_scale <- scale_y_continuous(
  limits=c(y_axis_min,y_axis_max),
  expand=c(0,0)
)

# BASE THEME

base_theme <- theme_classic() +
  theme(
    text=element_text(family="Arial",face="bold",size=22),
    axis.text=element_text(size=20),
    axis.title=element_text(size=24),
    axis.line=element_line(linewidth=2),
    axis.ticks=element_line(linewidth=2),
    axis.ticks.length=unit(0.3,"cm"),
    panel.grid=element_blank()
  )


# BOLUS SHADING LAYER


bolus_layer <- list(
  geom_rect(aes(xmin=-1, xmax=3.5, ymin=-Inf, ymax=Inf),
            fill="#cfe8ff", alpha=0.25, inherit.aes=FALSE),
  geom_rect(aes(xmin=11, xmax=15.5, ymin=-Inf, ymax=Inf),
            fill="#cfe8ff", alpha=0.25, inherit.aes=FALSE),
  geom_rect(aes(xmin=23, xmax=27.5, ymin=-Inf, ymax=Inf),
            fill="#cfe8ff", alpha=0.25, inherit.aes=FALSE)
)


# PLOT 1 — MEAN


p1 <- ggplot(data_long_plot,
             aes(x=day_plot,y=vitamin_d_level,group=day_plot)) +
  bolus_layer +
  all_brackets +
  geom_boxplot(width=2,linewidth=1.2,fatten=NULL) +
  stat_summary(fun=mean,geom="crossbar",
               width=2,linewidth=1.2,color="red") +
  geom_line(data=summary_stats,
            aes(x=day_plot,y=mean,group=1),
            color="red",linewidth=1.8,inherit.aes=FALSE) +
  geom_point(data=summary_stats,
             aes(x=day_plot,y=mean),
             color="red",size=4,inherit.aes=FALSE) +
  x_scale + y_scale +
  labs(x="Day", y="25(OH)D3 Level (ng/mL)") +
  base_theme

# PLOT 2 — MEDIAN


p2 <- ggplot(data_long_plot,
             aes(x=day_plot,y=vitamin_d_level,group=day_plot)) +
  bolus_layer +
  all_brackets +
  geom_boxplot(width=2,linewidth=1.2) +
  geom_line(data=summary_stats,
            aes(x=day_plot,y=median,group=1),
            color="darkblue",linewidth=1.8,inherit.aes=FALSE) +
  geom_point(data=summary_stats,
             aes(x=day_plot,y=median),
             color="darkblue",size=4,inherit.aes=FALSE) +
  x_scale + y_scale +
  labs(x="Day", y="25(OH)D3 Level (ng/mL)") +
  base_theme


# PLOT 3 — JITTER

p3 <- ggplot(data_long_plot,
             aes(x=day_plot,y=vitamin_d_level)) +
  bolus_layer +
  all_brackets +
  geom_jitter(width=0.8,size=3,alpha=0.5,color="gray40") +
  geom_line(data=summary_stats,
            aes(x=day_plot,y=mean,group=1),
            color="red",linewidth=2,inherit.aes=FALSE) +
  geom_point(data=summary_stats,
             aes(x=day_plot,y=mean),
             color="red",size=5,inherit.aes=FALSE) +
  x_scale + y_scale +
  labs(x="Day", y="25(OH)D3 Level (ng/mL)") +
  base_theme

# SAVE (PNG + SVG)

ggsave("vitd_boxplot_mean.png",plot=p1,width=16,height=12,dpi=300,bg="white")
ggsave("vitd_boxplot_median.png",plot=p2,width=16,height=12,dpi=300,bg="white")
ggsave("vitd_jitter.png",plot=p3,width=16,height=12,dpi=300,bg="white")

ggsave("vitd_boxplot_mean.svg",plot=p1,width=16,height=12,bg="white")
ggsave("vitd_boxplot_median.svg",plot=p2,width=16,height=12,bg="white")
ggsave("vitd_jitter.svg",plot=p3,width=16,height=12,bg="white")

print(p1)
print(p2)
print(p3)





## VITD MODeLLING SCRIPT
# 25(OH)D3 Pharmacokinetic Modeling — Decay and Supplementation Kinetics
# Phase 1  dC/dt = -d*C        =>  C(t) = C0 * exp(-d*t)
# Phase 2  dC/dt = k - d*C     =>  C(t) = (k/d) + (C0 - k/d)*exp(-d*t)


library(ggplot2)
library(cowplot)
library(svglite)

set.seed(as.integer(as.numeric(Sys.time()) * 1000) %% 10000)

DECAY_START  <- c(53.952, 45.694, 55.112)   # baseline C0  (ng/mL)
DECAY_END    <- c(38.103, 41.865, 42.774)   # follow-up Ct (ng/mL)
DECAY_DAYS   <- c(26.0,   30.0,   81.0)     # elapsed time (days)

INCREASE_DATA <- list(
  c(36.3, 42.2, 54.0),   # R1: C(0h), C(24h), C(48h)
  c(38.1, 45.9, 45.7),   # R2
  c(41.9, 48.2, 55.1)    # R3
)
TIME_POINTS <- c(0.0, 24.0, 48.0)   # hours

# Decay parameters from Phase 1 (per day), converted to per hour for Phase 2
DECAY_PER_DAY  <- c(0.013377, 0.002917, 0.003129)
DECAY_PER_HOUR <- DECAY_PER_DAY / 24.0

OUT_DIR <- "."


# PUBLICATION THEME

base_size <- 20

pub_theme <- function() {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      # Axes lines — only bottom and left, thick
      axis.line.x      = element_line(colour = "black", linewidth = 1.2),
      axis.line.y      = element_line(colour = "black", linewidth = 1.2),
      axis.ticks       = element_line(colour = "black", linewidth = 0.9),
      axis.ticks.length = unit(5, "pt"),
      
      # All text bold black 20 pt
      axis.text        = element_text(colour = "black", face = "bold",
                                      size = base_size),
      axis.title       = element_text(colour = "black", face = "bold",
                                      size = base_size),
      plot.title       = element_text(colour = "black", face = "bold",
                                      size = base_size, hjust = 0.5),
      
      # No grid, no panel border
      panel.grid       = element_blank(),
      panel.border     = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      
      # Legend OFF in all data plots (saved separately)
      legend.position  = "none",
      
      # Margins
      plot.margin = margin(10, 20, 10, 10)
    )
}

# Per-replicate colours
REP_COLS   <- c("royalblue", "firebrick", "forestgreen")
REP_LABELS <- c("Replicate 1", "Replicate 2", "Replicate 3")

# Helper: save both PNG and SVG
save_both <- function(p, stem, w, h) {
  ggsave(file.path(OUT_DIR, paste0(stem, ".png")),
         plot = p, width = w, height = h, dpi = 300, units = "in",
         bg = "white")
  ggsave(file.path(OUT_DIR, paste0(stem, ".svg")),
         plot = p, width = w, height = h, units = "in",
         bg = "white")
  cat("Saved:", stem, "\n")
}


# MATHEMATICAL MODELS 

# Exponential decay: C(t) = C0 * exp(-d*t)
VitD_decay <- function(C0, d, t) C0 * exp(-d * t)

# Analytical inversion: d = -ln(Ct/C0) / t
d_from_endpoints <- function(C0, Ct, t) -log(Ct / C0) / t

# ODE solution for supplementation
VitD_increase <- function(C0, k, d, t) {
  ss <- k / d
  ss + (C0 - ss) * exp(-d * t)
}

# Analytical k from a single follow-up point
k_from_point <- function(C0, C1, d, t1) {
  d * (C1 - C0 * exp(-d * t1)) / (1 - exp(-d * t1))
}

# Cost function for k optimisation (RSS over all time points)
cost_increase <- function(k, data, d, tpts) {
  C0  <- data[1]
  err <- 0.0
  for (i in 2:length(data)) {
    err <- err + (VitD_increase(C0, k, d, tpts[i]) - data[i])^2
  }
  err
}


# PHASE 1 — DECAY FITTING


cat("\n25(OH)D3 Decay Analysis\n")
cat(strrep("=", 60), "\n")

decay_results <- vector("list", 3)

for (ii in 1:3) {
  C0 <- DECAY_START[ii]
  Ct <- DECAY_END[ii]
  Dt <- DECAY_DAYS[ii]
  
  d_ana <- d_from_endpoints(C0, Ct, Dt)
  
  # Numerical cross-check via optim()
  obj   <- function(d) (VitD_decay(C0, d, Dt) - Ct)^2
  opt   <- tryCatch(
    optimise(obj, interval = c(max(1e-4, d_ana * 0.1),
                               min(10,   d_ana * 10))),
    error = function(e) list(minimum = d_ana, objective = 0)
  )
  d_num <- opt$minimum
  
  hl <- log(2) / d_ana
  
  decay_results[[ii]] <- list(
    rep          = ii,
    C0           = C0,
    Ct           = Ct,
    delta_t      = Dt,
    d_ana        = d_ana,
    d_num        = d_num,
    half_life    = hl
  )
  
  cat(sprintf("  R%d: d = %.6f day-1  |  t1/2 = %.1f days\n", ii, d_ana, hl))
}

d_day  <- sapply(decay_results, function(r) r$d_ana)
hl_day <- log(2) / d_day
cat(sprintf("\n  Mean d    = %.6f +/- %.6f day-1\n", mean(d_day), sd(d_day)))
cat(sprintf("  Mean t1/2 = %.1f +/- %.1f days\n\n",  mean(hl_day), sd(hl_day)))


# PHASE 2 — SUPPLEMENTATION FITTING


cat("25(OH)D3 Supplementation Kinetics\n")
cat(strrep("=", 60), "\n")

incr_results <- vector("list", 3)

for (ii in 1:3) {
  data <- INCREASE_DATA[[ii]]
  d    <- DECAY_PER_HOUR[ii]
  C0   <- data[1]
  
  k_24 <- k_from_point(C0, data[2], d, TIME_POINTS[2])
  k_48 <- k_from_point(C0, data[3], d, TIME_POINTS[3])
  
  # Numerical optimisation over k 
  obj   <- function(k) cost_increase(k, data, d, TIME_POINTS)
  opt   <- tryCatch(
    optimise(obj, interval = c(0, 10)),
    error = function(e) list(minimum = mean(c(k_24, k_48)), objective = NA)
  )
  k_opt <- opt$minimum
  ss    <- k_opt / d
  hl_h  <- log(2) / d
  
  incr_results[[ii]] <- list(
    rep         = ii,
    data        = data,
    d           = d,
    k_24        = k_24,
    k_48        = k_48,
    k_opt       = k_opt,
    ss          = ss,
    half_life_h = hl_h,
    half_life_d = hl_h / 24
  )
  
  cat(sprintf("  R%d: k = %.6f h-1  |  SS = %.1f ng/mL  |  t1/2 = %.2f days\n",
              ii, k_opt, ss, hl_h / 24))
}


# FIGURE 1 — DECAY CURVES  (3 stacked panels)

cat("\nBuilding Figure 1 — Decay curves...\n")

decay_panel_list <- list()

for (ii in 1:3) {
  r   <- decay_results[[ii]]
  col <- REP_COLS[ii]
  
  t_end <- r$delta_t * 1.15
  t_seq <- seq(0, t_end, length.out = 1000)
  C_fit <- VitD_decay(r$C0, r$d_ana, t_seq)
  
  df_line <- data.frame(t = t_seq, C = C_fit)
  df_pts  <- data.frame(t = c(0, r$delta_t), C = c(r$C0, r$Ct))
  
  y_lo <- r$Ct  - 0.12 * (r$C0 - r$Ct)
  y_hi <- r$C0  + 0.12 * (r$C0 - r$Ct)
  
  p <- ggplot() +
    geom_line(data = df_line, aes(x = t, y = C),
              colour = col, linewidth = 1.2) +
    geom_point(data = df_pts, aes(x = t, y = C),
               colour = col, fill = col, shape = 21,
               size = 4, stroke = 1.5) +
    scale_x_continuous(expand = c(0.02, 0)) +
    coord_cartesian(xlim = c(0, t_end), ylim = c(y_lo, y_hi)) +
    labs(
      x     = "Time (days)",
      y     = "25(OH)D3 (ng/mL)",
      title = sprintf("Replicate %d   d = %.5f day-1   t1/2 = %.1f days",
                      ii, r$d_ana, r$half_life)
    ) +
    pub_theme()
  
  decay_panel_list[[ii]] <- p
}

fig1 <- plot_grid(plotlist = decay_panel_list,
                  ncol = 1, align = "v", axis = "lr")

save_both(fig1, "fig1_25OHD3_decay_replicates", w = 8, h = 12)


# FIGURE 2 — SUPPLEMENTATION KINETICS  (3 stacked panels)
# Solid = k_opt fit.  Dashed = steady-state asymptote.
# Dotted lines = k from individual 24 h and 48 h points.
# Filled circles = observed data.

cat("Building Figure 2 — Supplementation kinetics...\n")

incr_panel_list <- list()

for (ii in 1:3) {
  r   <- incr_results[[ii]]
  col <- REP_COLS[ii]
  
  t_seq  <- seq(0, 60, length.out = 1000)
  C_kopt <- VitD_increase(r$data[1], r$k_opt, r$d, t_seq)
  C_k24  <- VitD_increase(r$data[1], r$k_24,  r$d, t_seq)
  C_k48  <- VitD_increase(r$data[1], r$k_48,  r$d, t_seq)
  
  df_kopt <- data.frame(t = t_seq, C = C_kopt, curve = "k optimised")
  df_k24  <- data.frame(t = t_seq, C = C_k24,  curve = "k from 24 h")
  df_k48  <- data.frame(t = t_seq, C = C_k48,  curve = "k from 48 h")
  df_line <- rbind(df_kopt, df_k24, df_k48)
  
  df_pts  <- data.frame(t = TIME_POINTS, C = r$data)
  
  dr   <- diff(range(r$data))
  y_lo <- min(r$data) - 0.15 * dr
  y_hi <- max(r$data) + 0.35 * dr
  
  p <- ggplot() +
    # k from 24 h (dotted)
    geom_line(data = df_k24, aes(x = t, y = C),
              colour = col, linewidth = 0.9, linetype = "dotted") +
    # k from 48 h (dashed)
    geom_line(data = df_k48, aes(x = t, y = C),
              colour = col, linewidth = 0.9, linetype = "dashed") +
    # k optimised (solid, thicker)
    geom_line(data = df_kopt, aes(x = t, y = C),
              colour = col, linewidth = 1.4, linetype = "solid") +
    # Steady-state asymptote
    geom_hline(yintercept = r$ss,
               colour = col, linewidth = 0.9, linetype = "longdash") +
    # Observed data points
    geom_point(data = df_pts, aes(x = t, y = C),
               colour = col, fill = col, shape = 21,
               size = 4, stroke = 1.5) +
    scale_x_continuous(expand = c(0.02, 0)) +
    coord_cartesian(xlim = c(0, 60), ylim = c(y_lo, y_hi)) +
    labs(
      x     = "Time (hours)",
      y     = "25(OH)D3 (ng/mL)",
      title = sprintf("Replicate %d   k = %.5f h-1   SS = %.1f ng/mL",
                      ii, r$k_opt, r$ss)
    ) +
    pub_theme()
  
  incr_panel_list[[ii]] <- p
}

fig2 <- plot_grid(plotlist = incr_panel_list,
                  ncol = 1, align = "v", axis = "lr")

save_both(fig2, "fig2_25OHD3_supplementation_replicates", w = 8, h = 12)

# FIGURE 3 — COMBINED: DECAY + SUPPLEMENTATION SIDE BY SIDE

cat("Building Figure 3 — Combined phases side-by-side...\n")

combined_panels <- list()
panel_idx <- 1

for (ii in 1:3) {
  dr  <- decay_results[[ii]]
  ir  <- incr_results[[ii]]
  col <- REP_COLS[ii]
  
  # --- Left: decay ---
  t_end <- dr$delta_t * 1.15
  t_seq <- seq(0, t_end, length.out = 1000)
  C_fit <- VitD_decay(dr$C0, dr$d_ana, t_seq)
  
  df_dl <- data.frame(t = t_seq, C = C_fit)
  df_dp <- data.frame(t = c(0, dr$delta_t), C = c(dr$C0, dr$Ct))
  
  y_lo_d <- dr$Ct - 0.12 * (dr$C0 - dr$Ct)
  y_hi_d <- dr$C0 + 0.12 * (dr$C0 - dr$Ct)
  
  pL <- ggplot() +
    geom_line(data = df_dl, aes(x = t, y = C),
              colour = col, linewidth = 1.2) +
    geom_point(data = df_dp, aes(x = t, y = C),
               colour = col, fill = col, shape = 21,
               size = 4, stroke = 1.5) +
    scale_x_continuous(expand = c(0.02, 0)) +
    coord_cartesian(xlim = c(0, t_end), ylim = c(y_lo_d, y_hi_d)) +
    labs(
      x     = "Time (days)",
      y     = "25(OH)D3 (ng/mL)",
      title = sprintf("R%d — Decay", ii)
    ) +
    pub_theme()
  
  # --- Right: supplementation ---
  t_seq_h <- seq(0, 60, length.out = 1000)
  C_i     <- VitD_increase(ir$data[1], ir$k_opt, ir$d, t_seq_h)
  
  df_il <- data.frame(t = t_seq_h, C = C_i)
  df_ip <- data.frame(t = TIME_POINTS, C = ir$data)
  
  dr_i  <- diff(range(ir$data))
  y_lo_i <- min(ir$data) - 0.15 * dr_i
  y_hi_i <- max(ir$data) + 0.35 * dr_i
  
  pR <- ggplot() +
    geom_line(data = df_il, aes(x = t, y = C),
              colour = col, linewidth = 1.2) +
    geom_hline(yintercept = ir$ss,
               colour = col, linewidth = 0.9, linetype = "dashed") +
    geom_point(data = df_ip, aes(x = t, y = C),
               colour = col, fill = col, shape = 21,
               size = 4, stroke = 1.5) +
    scale_x_continuous(expand = c(0.02, 0)) +
    coord_cartesian(xlim = c(0, 60), ylim = c(y_lo_i, y_hi_i)) +
    labs(
      x     = "Time (hours)",
      y     = "25(OH)D3 (ng/mL)",
      title = sprintf("R%d — Supplementation", ii)
    ) +
    pub_theme()
  
  combined_panels[[panel_idx]]     <- pL
  combined_panels[[panel_idx + 1]] <- pR
  panel_idx <- panel_idx + 2
}

fig3 <- plot_grid(plotlist = combined_panels,
                  ncol = 2, align = "hv", axis = "tblr")

save_both(fig3, "fig3_25OHD3_combined_phases", w = 14, h = 12)


# FIGURE 4 — PARAMETER SUMMARY BARS  (2x2 grid)

cat("Building Figure 4 — Parameter summary bars...\n")

k_hr    <- sapply(incr_results, function(r) r$k_opt)
ss_vals <- sapply(incr_results, function(r) r$ss)

df_params <- data.frame(
  rep   = factor(1:3, labels = c("R1","R2","R3")),
  col   = REP_COLS,
  d     = d_day,
  hl    = hl_day,
  k     = k_hr,
  ss    = ss_vals
)

make_bar <- function(df, y_col, ylab, ttl) {
  y_max <- max(df[[y_col]]) * 1.25
  ggplot(df, aes(x = rep, y = .data[[y_col]], fill = rep)) +
    geom_col(colour = "black", linewidth = 0.8, width = 0.6) +
    geom_point(shape = 21, size = 4, fill = "black",
               colour = "white", stroke = 1.5) +
    scale_fill_manual(values = REP_COLS) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, y_max)) +
    labs(x = "Replicate", y = ylab, title = ttl) +
    pub_theme() +
    theme(legend.position = "none")
}

pA <- make_bar(df_params, "d",  "d (day-1)",   "Decay rate constant")
pB <- make_bar(df_params, "hl", "t1/2 (days)", "25(OH)D3 half-life")
pC <- make_bar(df_params, "k",  "k (h-1)",     "Production rate constant")
pD <- make_bar(df_params, "ss", "SS (ng/mL)",  "Predicted steady state")

fig4 <- plot_grid(pA, pB, pC, pD, ncol = 2, align = "hv", axis = "tblr")

save_both(fig4, "fig4_25OHD3_parameter_summary", w = 12, h = 10)


# LEGEND FIGURE — saved separately, not attached to any data figure


cat("Building legend figure...\n")

# Replicate colour legend
df_leg_rep <- data.frame(
  rep  = factor(1:3, labels = REP_LABELS),
  x    = 1:3,
  y    = rep(1, 3),
  col  = REP_COLS
)

leg_rep <- ggplot(df_leg_rep, aes(x = x, y = y, colour = rep)) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 4, stroke = 1.5) +
  scale_colour_manual(
    values = setNames(REP_COLS, REP_LABELS),
    name   = "Subject"
  ) +
  theme_void(base_family = "Arial") +
  theme(
    legend.position  = "center",
    legend.title     = element_text(face = "bold", size = base_size,
                                    colour = "black"),
    legend.text      = element_text(face = "bold", size = base_size,
                                    colour = "black"),
    legend.key.width = unit(2, "cm"),
    legend.key.size  = unit(1, "cm"),
    plot.background  = element_rect(fill = "white", colour = NA)
  )

# Curve-type legend
df_leg_curve <- data.frame(
  ltype = factor(
    c("Fitted model (k opt)",
      "k from 24 h point",
      "k from 48 h point",
      "Steady state",
      "Observed data"),
    levels = c("Fitted model (k opt)",
               "k from 24 h point",
               "k from 48 h point",
               "Steady state",
               "Observed data")
  ),
  x = 1:5, y = rep(1, 5)
)

leg_curve <- ggplot(df_leg_curve, aes(x = x, y = y,
                                      linetype = ltype,
                                      shape    = ltype)) +
  geom_line(colour = "black", linewidth = 1.1, na.rm = TRUE) +
  geom_point(colour = "black", fill = "black", size = 4,
             data = df_leg_curve[df_leg_curve$ltype == "Observed data", ],
             na.rm = TRUE) +
  scale_linetype_manual(
    name   = "Curve type",
    values = c(
      "Fitted model (k opt)" = "solid",
      "k from 24 h point"    = "dotted",
      "k from 48 h point"    = "dashed",
      "Steady state"         = "longdash",
      "Observed data"        = "blank"
    )
  ) +
  scale_shape_manual(
    name   = "Curve type",
    values = c(
      "Fitted model (k opt)" = NA,
      "k from 24 h point"    = NA,
      "k from 48 h point"    = NA,
      "Steady state"         = NA,
      "Observed data"        = 21
    )
  ) +
  theme_void(base_family = "Arial") +
  theme(
    legend.position  = "center",
    legend.title     = element_text(face = "bold", size = base_size,
                                    colour = "black"),
    legend.text      = element_text(face = "bold", size = base_size,
                                    colour = "black"),
    legend.key.width = unit(2, "cm"),
    legend.key.size  = unit(1, "cm"),
    plot.background  = element_rect(fill = "white", colour = NA)
  )

leg_rep_grob   <- get_legend(leg_rep)
leg_curve_grob <- get_legend(leg_curve)

fig_legend <- plot_grid(leg_rep_grob, leg_curve_grob,
                        ncol = 1, rel_heights = c(1, 1.6))

save_both(fig_legend, "fig_legend", w = 8, h = 5)


# TXT OUTPUT — DECAY EQUATIONS

sink(file.path(OUT_DIR, "25OHD3_decay_equations.txt"))
cat("25(OH)D3 Decay Model — Fitted Equations\n")
cat(strrep("=", 65), "\n")
cat("ODE:      dC/dt = -d*C\n")
cat("Solution: C(t)  = C0 * exp(-d * t)\n")
cat("Units:    C ng/mL,  t days,  d day^-1\n\n")

for (r in decay_results) {
  cat(strrep("-", 65), "\n")
  cat(sprintf("Replicate %d\n", r$rep))
  cat(strrep("-", 65), "\n")
  cat(sprintf("  C0            = %8.3f ng/mL\n",  r$C0))
  cat(sprintf("  Ct            = %8.3f ng/mL\n",  r$Ct))
  cat(sprintf("  dt            = %8.1f days\n",   r$delta_t))
  cat(sprintf("  d (analytical)= %12.6f day^-1\n", r$d_ana))
  cat(sprintf("  d (numerical) = %12.6f day^-1\n", r$d_num))
  cat(sprintf("  t1/2          = %8.2f days\n",   r$half_life))
  cat(sprintf("  Equation: C(t) = %.3f * exp( -%.6f * t )\n\n",
              r$C0, r$d_ana))
}

cat(strrep("-", 65), "\n")
cat("Summary Statistics\n")
cat(strrep("-", 65), "\n")
cat(sprintf("  Mean d    = %.6f +/- %.6f day^-1\n", mean(d_day), sd(d_day)))
cat(sprintf("  Mean t1/2 = %.2f  +/- %.2f  days\n",  mean(hl_day), sd(hl_day)))
cat(sprintf("\nGenerated: %s\n", Sys.time()))
sink()
cat("Saved: 25OHD3_decay_equations.txt\n")


# TXT OUTPUT — SUPPLEMENTATION EQUATIONS

sink(file.path(OUT_DIR, "25OHD3_supplementation_equations.txt"))
cat("25(OH)D3 Supplementation Kinetics — Fitted Equations\n")
cat(strrep("=", 65), "\n")
cat("ODE:      dC/dt = k - d*C\n")
cat("Solution: C(t)  = (k/d) + [C0 - k/d] * exp(-d * t)\n")
cat("Units:    C ng/mL,  t hours,  d h^-1,  k ng/mL/h\n\n")

for (r in incr_results) {
  cat(strrep("-", 65), "\n")
  cat(sprintf("Replicate %d\n", r$rep))
  cat(strrep("-", 65), "\n")
  cat(sprintf("  C(0h)              = %8.3f ng/mL\n",  r$data[1]))
  cat(sprintf("  C(24h)             = %8.3f ng/mL\n",  r$data[2]))
  cat(sprintf("  C(48h)             = %8.3f ng/mL\n",  r$data[3]))
  cat(sprintf("  d (h^-1)           = %12.8f\n",       r$d))
  cat(sprintf("  k from 24h         = %12.6f h^-1\n",  r$k_24))
  cat(sprintf("  k from 48h         = %12.6f h^-1\n",  r$k_48))
  cat(sprintf("  k optimised        = %12.6f h^-1\n",  r$k_opt))
  cat(sprintf("  Steady state (k/d) = %8.2f ng/mL\n",  r$ss))
  cat(sprintf("  t1/2               = %8.2f days\n",   r$half_life_d))
  cat(sprintf("  Equation: C(t) = %.2f + (%.3f - %.2f) * exp(-%.8f * t)\n\n",
              r$ss, r$data[1], r$ss, r$d))
  
  p24 <- VitD_increase(r$data[1], r$k_opt, r$d, 24.0)
  p48 <- VitD_increase(r$data[1], r$k_opt, r$d, 48.0)
  cat("  Validation:\n")
  cat(sprintf("    C(24h): pred=%.3f  obs=%.3f  err=%.3f ng/mL\n",
              p24, r$data[2], abs(p24 - r$data[2])))
  cat(sprintf("    C(48h): pred=%.3f  obs=%.3f  err=%.3f ng/mL\n\n",
              p48, r$data[3], abs(p48 - r$data[3])))
}

k_vec  <- sapply(incr_results, function(r) r$k_opt)
ss_vec <- sapply(incr_results, function(r) r$ss)
cat(strrep("-", 65), "\n")
cat("Summary Statistics\n")
cat(strrep("-", 65), "\n")
cat(sprintf("  Mean k  (h^-1)   = %.6f +/- %.6f\n",    mean(k_vec),    sd(k_vec)))
cat(sprintf("  Mean k  (day^-1) = %.4f  +/- %.4f\n",   mean(k_vec)*24, sd(k_vec)*24))
cat(sprintf("  Mean SS (ng/mL)  = %.2f   +/- %.2f\n",  mean(ss_vec),   sd(ss_vec)))
cat(sprintf("\nGenerated: %s\n", Sys.time()))
sink()
cat("Saved: 25OHD3_supplementation_equations.txt\n")

# CONSOLE SUMMARY TABLE

cat("\n", strrep("=", 72), "\n")
cat("FINAL PARAMETER SUMMARY — 25(OH)D3 MODELING\n")
cat(strrep("=", 72), "\n")
cat(sprintf("%-5s  %-14s %-12s %-14s %-12s %-12s\n",
            "Rep", "d (day^-1)", "t1/2 (days)", "k (h^-1)",
            "k (day^-1)", "SS (ng/mL)"))
cat(strrep("-", 72), "\n")
for (ii in 1:3) {
  dr <- decay_results[[ii]]
  ir <- incr_results[[ii]]
  cat(sprintf("R%-4d  %-14.6f %-12.2f %-14.6f %-12.4f %-12.2f\n",
              ii, dr$d_ana, dr$half_life, ir$k_opt, ir$k_opt*24, ir$ss))
}
cat(strrep("-", 72), "\n")
cat(sprintf("Mean   %-14.6f %-12.2f %-14.6f %-12.4f %-12.2f\n",
            mean(d_day), mean(hl_day), mean(k_hr), mean(k_hr)*24, mean(ss_vals)))
cat(sprintf("+/-SD  %-14.6f %-12.2f %-14.6f %-12.4f %-12.2f\n",
            sd(d_day), sd(hl_day), sd(k_hr), sd(k_hr)*24, sd(ss_vals)))
cat(strrep("=", 72), "\n")
cat(sprintf("Completed: %s\n", Sys.time()))

