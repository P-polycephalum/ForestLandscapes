library(pacman)
p_load(dclone, MASS, ggplot2, snow, tidyverse, parallel, gridExtra)

diag_test <- function(diag_table, tol = 0.10) {
  diag_table$standard_eigenvalue <- diag_table$lambda.max / diag_table$lambda.max[1]
  diag_table$expected             <- 1 / diag_table$n.clones
  diag_table$ratio                <- diag_table$standard_eigenvalue / diag_table$expected
  diag_table$signed_deviation     <- diag_table$ratio - 1
  diag_table$abs_deviation        <- abs(diag_table$signed_deviation)
  diag_table$within_expected_band <- diag_table$abs_deviation <= tol
  
  ratios <- diag_table$ratio[diag_table$n.clones > 1]
  signed_deviation <- ratios - 1
  
  rate_metric_abs    <- mean(abs(signed_deviation))
  rate_metric_signed <- mean(signed_deviation)
  max_deviation      <- max(abs(signed_deviation))
  max_positive_dev   <- max(signed_deviation)
  
  estimable <- max_deviation <= tol
  
  warning <- if (estimable) {
    "expected_1_over_K_decline"
  } else if (max_positive_dev > tol) {
    "warning_slower_than_expected_possible_nonestimability"
  } else {
    "warning_irregular_or_faster_than_expected_decline"
  }
  
  list(
    table = diag_table,
    rate_metric_abs = rate_metric_abs,
    rate_metric_signed = rate_metric_signed,
    max_deviation = max_deviation,
    max_positive_deviation = max_positive_dev,
    estimable = estimable,
    warning = warning
  )
}
day_to_DOY <- function(day, pheno_year, start_month = pheno_start_month) {
  if (is.na(start_month)) start_month <- 1
  
  format(
    as.Date(
      paste0(
        as.numeric(as.character(pheno_year)),
        "-",
        sprintf("%02d", start_month),
        "-01"
      )
    ) + day - 1,
    "%b %d"
  )
}
DOY_to_day <- function(month, day, pheno_year, start_month = pheno_start_month) {
  if (is.na(start_month)) start_month <- 1
  cal_year <- if (month >= start_month) pheno_year else pheno_year + 1
  as.numeric(
    as.Date(paste0(cal_year, "-", sprintf("%02d", month), "-", sprintf("%02d", day))) -
      as.Date(paste0(pheno_year, "-", sprintf("%02d", start_month), "-01"))
  )
}
dur_trans <- function(model,type="90_10"){
  kd <- summary(model)$statistics["kd", "Mean"]
  transition_width_90_10 <- 2 * log(9) / kd
  transition_width_99_01 <- 2 * log(90) / kd
  print(paste0("Transition width (90-10): ", round(transition_width_90_10, 2), " days"))
  print(paste0("Transition width (99-01): ", round(transition_width_99_01, 2), " days"))
}
mid_trans <- function(model){
  Td <- summary(model)$statistics["Td", "Mean"]
  print(Td)
}
dur_trans <- function(model, type="90_10"){
  kd <- summary(model)$statistics["kd", "Mean"]
  transition_width_90_10 <- 2 * log(9) / kd
  transition_width_99_01 <- 2 * log(90) / kd
  if (type == "90_10") {
    return(transition_width_90_10)
  } else if (type == "99_01") {
    return(transition_width_99_01)
  }
}
mid_trans <- function(model){
  Td <- summary(model)$statistics["Td", "Mean"]
  return(Td)
}
# Read and prepare data
# 8 for drop in cava, 11 for drop in dipteryx
pheno_start_month <- 10
data <- read.csv("timeseriesv2//data//dipteryx_oleifera_leafing_valid.csv")
data <- data %>%
  mutate(
    y_norm = pmin(pmax(leafing / 100, 1e-4), 1 - 1e-4),
    date = as.Date(time),
    DOY = yday(date),
    year = year(date),
    month = month(date),
    pheno_year = if (is.na(pheno_start_month)) year else if_else(month >= pheno_start_month, year, year - 1),
    day = if (is.na(pheno_start_month)) DOY else as.numeric(difftime(
      date,
      as.Date(paste0(pheno_year, "-", sprintf("%02d", pheno_start_month), "-01")),
      units = "days"
    )) + 1L,
    tree = as.factor(tag),
    pheno_year = as.factor(pheno_year),
    tree_year = as.factor(paste0(tree, "_", pheno_year)))

trees1 <- unique(data$tree)
years1 <- unique(data$pheno_year)
all_before_threshold <- data.frame()
all_after_threshold <- data.frame()
for (i in 1:length(trees1)) {
  for (j in 1:length(years1)) {
    subset_data <- data %>%
      filter(tree == trees1[i], pheno_year == years1[j]) %>%
      arrange(day)
    if (nrow(subset_data) == 0 || all(is.na(subset_data$y_norm))) {
      next
    }

    min_indices <- which(subset_data$y_norm == min(subset_data$y_norm, na.rm = TRUE))
    cut_idx <- min_indices[1]
    cut_idx_end <- min_indices[length(min_indices)]
    seg <- subset_data[1:cut_idx_end, ]
    seg_after <- subset_data[cut_idx:nrow(subset_data), ]
    
    #loop for drop
    n_obs_before   <- nrow(seg)
    n_obs_after    <- nrow(seg_after)
    y_range_before <- max(seg$y_norm, na.rm = TRUE) - min(seg$y_norm, na.rm = TRUE)
    y_range_after  <- max(seg_after$y_norm, na.rm = TRUE) - min(seg_after$y_norm, na.rm = TRUE)

    if (n_obs_before < 2 || y_range_before < 0.5) {
      print(paste0("SKIPPED tree=", trees1[i], " year=", years1[j],
                   " n=", n_obs_before, " range=", round(y_range_before, 3)))
      next
    }
    if (n_obs_after < 2 || y_range_after < 0.5) {
      print(paste0("SKIPPED tree=", trees1[i], " year=", years1[j],
                   " n=", n_obs_after, " range=", round(y_range_after, 3)))
      next
    }
    print(paste0("tree=", trees1[i], " year=", years1[j],
                 " n=", n_obs_before, " range=", round(y_range_before, 3)))
    print(paste0("tree=", trees1[i], " year=", years1[j],
                 " n=", n_obs_after, " range=", round(y_range_after, 3)))
    all_before_threshold <- bind_rows(all_before_threshold, seg)
    all_after_threshold <- bind_rows(all_after_threshold, seg_after)
  }
}

drop_seg <- all_before_threshold %>%
  filter(!tree %in% c("3811", "4250")) %>%
  mutate(
    pheno_year = droplevels(as.factor(pheno_year)),
    tree = droplevels(as.factor(tree))
  )

flush_seg <- all_after_threshold %>%
  filter(!pheno_year %in% c("2019")) %>%
  mutate(
    pheno_year = droplevels(as.factor(pheno_year)),
    tree = droplevels(as.factor(tree))
  )
  
drop_year_id <- as.numeric(as.factor(drop_seg$pheno_year))
drop_tree_id <- as.numeric(as.factor(drop_seg$tree))
flush_year_id <- as.numeric(as.factor(flush_seg$pheno_year))
flush_tree_id <- as.numeric(as.factor(flush_seg$tree))

#cava models
model.drop<- readRDS("timeseriesv2//data//cava_drop.rds") #8
model.flush <- readRDS("timeseriesv2//data//cava_flush.rds") # 10
model.drop_year_fe <- readRDS("timeseriesv2//data//cava_drop_year_fe.rds") #i think we skip 2017 #8
model.flush_year_fe <- readRDS("timeseriesv2//data//cava_flush_year_fe.rds") #10 skipped 2019
model.drop_indv_fe <- readRDS("timeseriesv2//data//cava_drop_indv_fe.rds") # processing now, 8
model.flush_indv_fe <- readRDS("timeseriesv2//data//cava_flush_indv_fe.rds") # processing now, 10
model.drop_yearfe_indvre<- readRDS("timeseriesv2//data//cava_drop_yearfe_indvre.rds") # processing now, 8

 #dipterix models
dipt.drop <- readRDS("timeseriesv2//data//dipt_drop.rds") #10
dipt.flush <- readRDS("timeseriesv2//data//dipt_flush.rds") #10 also, seems to fit fine
dipt.drop_year_fe <- readRDS("timeseriesv2//data//dipt_drop_year_fe.rds") # we skip 2017
dipt.flush_year_fe <- readRDS("timeseriesv2//data//dipt_flush_year_fe.rds") #we skip 2019
dipt.drop_indv_fe <- readRDS("timeseriesv2//data//dipt_drop_indv_fe.rds") # we skip trees 3811 and 4250
dipt.flush_indv_fe <- readRDS("timeseriesv2//data//dipt_flush_indv_fe.rds") # we skip trees 3811 and 4250


model.drop_indv_fe <- readRDS("timeseriesv2//data//dipt_drop_indv_fe.rds") # we skip trees 3811 and 4250
base<- summary(model.drop_indv_fe)$statistics["base", "Mean"]
amp<- summary(model.drop_indv_fe)$statistics["amp", "Mean"]
kd<- summary(model.drop_indv_fe)$statistics["kd", "Mean"]
uniqueTrees<- unique(drop_seg$tree)
iTd<- c()
for (i in 1:length(uniqueTrees)){
  iTd[i]<- summary(model.drop_indv_fe)$statistics[paste0("iTd[", i, "]"), "Mean"]
}
df_drop_indv <- data.frame()
for (i in 1:length(uniqueTrees)) {
  days <- seq(min(drop_seg$day), max(drop_seg$day), by = 1)
  y_norm <- base + amp / (1 + exp(kd * (days - iTd[i])))
  date <- as.Date(day_to_DOY(days, drop_seg$pheno_year[drop_seg$tree == uniqueTrees[i]][1], start_month = 11), format="%B %d")
  df_drop_indv <- rbind(df_drop_indv, data.frame(day = days, tree = uniqueTrees[i], y_norm = y_norm, date = date, transition="drop"))
}

ggplot(drop_seg, aes(x = day, y = y_norm, color=tree_year)) +
  geom_point()+
  geom_line(data = df_drop_indv, aes(x = day, y = y_norm, color=tree))+
  facet_wrap(~tree)+
  theme(legend.position = "none")


#Code to create interannual shifts plot
#define the drop and flush models as
drop<- dipt.drop_year_fe
flush<- dipt.flush_year_fe
dur_trans(drop)
summary(drop)
dur_trans(flush)
amp_flush <- summary(flush)$statistics["amp", "Mean"]
base_flush <- summary(flush)$statistics["base", "Mean"]
kd_flush <- summary(flush)$statistics["kd", "Mean"]
amp <- summary(drop)$statistics["amp", "Mean"]
base <- summary(drop)$statistics["base", "Mean"]
kd <- summary(drop)$statistics["kd", "Mean"]
amp<- amp_flush
base<- base_flush

years_drop <- sort(unique(drop_seg$pheno_year))
years_flush <- sort(unique(flush_seg$pheno_year))
yTd_drop <- summary(drop)$statistics[
  paste0("yTd[", seq_along(years_drop), "]"),
  "Mean"
]
yTd_flush <- summary(flush)$statistics[
  paste0("yTd[", seq_along(years_flush), "]"),
  "Mean"
]

#one curve per year 
df_drop <- data.frame()
for (year in years_drop) {
  yTd <- yTd_drop[which(years_drop == year)]
  days <- seq(min(drop_seg$day), max(drop_seg$day), by = 1)
  y_norm <- base + amp / (1 + exp(kd * (days - yTd)))
  date <- as.Date(day_to_DOY(days, year, start_month = 8), format="%B %d")
  df_drop <- rbind(df_drop, data.frame(day = days, pheno_year = year, y_norm = y_norm, date = date, transition="drop"))
}

df_flush <- data.frame()
for (year in years_flush) {
  yTd <- yTd_flush[which(years_flush == year)]
  days <- seq(min(flush_seg$day), max(flush_seg$day), by = 1)
  y_norm <- base_flush + amp_flush / (1 + exp(-kd_flush * (days - yTd)))
  date <- as.Date(day_to_DOY(days, year, start_month = 10), format="%B %d")
  df_flush <- rbind(df_flush, data.frame(day = days, pheno_year = year, y_norm = y_norm, date = date, transition="flush"))
}


limit_drop_left <- min(yTd_drop) - dur_trans(drop) -10
limit_drop_right <- max(yTd_drop) + dur_trans(drop) +10
limit_flush_left <- min(yTd_flush) - dur_trans(flush) -10
limit_flush_right <- max(yTd_flush) + dur_trans(flush) +10

tmp <- (limit_drop_right+limit_flush_left)/2
limit_drop_right <- tmp
limit_flush_left <- tmp

merged_drop<- left_join(df_drop, drop_seg[, c("day", "pheno_year", "y_norm")], by=c("day", "pheno_year")) %>%
  filter(day>=limit_drop_left & day<limit_drop_right)

merged_flush<- left_join(df_flush, flush_seg[, c("day", "pheno_year", "y_norm")], by=c("day", "pheno_year")) %>%
  filter(day >= limit_flush_left & day < limit_flush_right)
combined_df <- rbind(merged_drop, merged_flush)


years <- sort(unique(as.character(combined_df$pheno_year)))

year_labels <- setNames(
  paste0(years, "-", as.numeric(years) + 1),
  years
)
transition_labels <- c(
  "drop" = "Leaf Drop",
  "flush" = "Leaf Flush"
)
windows()
ggplot(combined_df, aes(x = day, y = y_norm.x, color = pheno_year)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(y = y_norm.y), alpha = 0.3) +
  facet_grid(
    pheno_year ~ transition,
    scales = "free_x",
    drop = FALSE,
    labeller = labeller(pheno_year = year_labels, transition = transition_labels)
  )+
  scale_x_continuous(
    breaks = c(seq(0, 365, by = 30)),
    labels = function(x) day_to_DOY(x, 2021, start_month = ifelse(combined_df$transition[1] == "drop", 8, 10))
  ) + theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(y = "Leaf Cover",
       title = "Interannual variation in Leaf Drop and Flush of Dipteryx oleifera",
       subtitle = "Data points (faded) and fitted curves (solid lines) for each phenological year") 

ggsave("timeseriesv2//figures//dipt_drop_fixed_effect.jpg", width = 10, height = 8, dpi = 300)
#to create a plot of dipteryx drop and flush we will
#pull the data points of the drop cycle
#parameters needed
amp <- summary(dipt.drop)$statistics["amp", "Mean"]
amp_flush <- summary(dipt.flush)$statistics["amp", "Mean"]
base <- summary(dipt.drop)$statistics["base", "Mean"]
base_flush <- summary(dipt.flush)$statistics["base", "Mean"]
Td <- summary(dipt.drop)$statistics["Td", "Mean"]
Td_flush <- summary(dipt.flush)$statistics["Td", "Mean"]
kd <- summary(dipt.drop)$statistics["kd", "Mean"]
kd_flush <- summary(dipt.flush)$statistics["kd", "Mean"]

midpoint<- (Td_flush+ Td)/2

df_drop<- data.frame(day=seq(min(drop_seg$day), round(midpoint, 0), by=1))
for (i in 1:length(df_drop$day)){
  df_drop$y_norm[i] <- (amp/(1+exp(kd*(df_drop$day[i]-Td))))+base
  df_drop$date[i]<-as.Date(paste(day_to_DOY(df_drop$day[i], 2019, start_month = pheno_start_month), "2019"), format="%B %d %Y")
}

df_flush<- data.frame(day=seq(round(midpoint, 0), max(flush_seg$day), by=1))
for (i in 1:length(df_flush$day)){
  df_flush$y_norm[i] <- (amp_flush/(1+exp(-kd_flush*(df_flush$day[i]-Td_flush))))+base_flush
  df_flush$date[i] <- as.Date(paste(day_to_DOY(df_flush$day[i], 2019, start_month = pheno_start_month), "2019"), format="%B %d %Y")
}

x_left <- Td - round(dur_trans(dipt.drop)/2, 0)
x_right <- Td + round(dur_trans(dipt.drop)/2, 0)
x_left_flush <- Td_flush - round(dur_trans(dipt.flush)/2, 0)
x_right_flush <- Td_flush + round(dur_trans(dipt.flush)/2, 0)
drop_panel <- ggplot(drop_seg[drop_seg$day <= round(midpoint, 0), ], aes(x = day, y = y_norm, color = tree_year)) +
  geom_point(alpha = 0.3) +
  geom_line(data = df_drop, aes(x = day, y = y_norm), color = "red", size = 1) +
  geom_vline(xintercept = Td, color = "black", linetype = "dashed") +
  geom_vline(xintercept = x_left, color = "black", linetype = "dashed") +
  geom_vline(xintercept = x_right, color = "black", linetype = "dashed") +
  geom_segment(aes(x=x_left+1, xend=x_right-1, y=0.5, yend=0.5), inherit.aes = FALSE, color="black", linewidth=0.5)+
  geom_segment(aes(x=x_left+1, xend=x_left+1, y=0.48, yend=0.52), inherit.aes = FALSE, color="black", linewidth=0.5)+
  geom_segment(aes(x=x_right-1, xend=x_right-1, y=0.48, yend=0.52), inherit.aes = FALSE, color="black", linewidth=0.5)+
  geom_text(aes(x=x_left-15, y=0.50, label=paste0(round(dur_trans(dipt.drop), 1), " days")), inherit.aes = FALSE, color="black", size=4)+
  scale_x_continuous(breaks = c(seq(0, 365, by = 30)), labels = function(x) day_to_DOY(x, 2019, start_month = pheno_start_month)) +
  ylab("Leaf Cover") +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
  )

flush_panel<-ggplot(flush_seg[flush_seg$day >= round(midpoint, 0), ], aes(x = day, y = y_norm, color=tree_year)) +
  geom_point(alpha=0.3)+
  geom_line(data=df_flush, aes(x=day, y=y_norm), color="blue", size=1)+
  geom_vline(xintercept=Td_flush, color="black", linetype="dashed")+
  geom_vline(xintercept=x_left_flush, color="black", linetype="dashed")+
  geom_vline(xintercept=x_right_flush, color="black", linetype="dashed")+
  geom_segment(aes(x=x_left_flush+1, xend=x_right_flush-1, y=0.5, yend=0.5), inherit.aes = FALSE, color="black", linewidth=0.5)+
  geom_segment(aes(x=x_left_flush+1, xend=x_left_flush+1, y=0.48, yend=0.52), inherit.aes = FALSE, color="black", linewidth=0.5)+
  geom_segment(aes(x=x_right_flush-1, xend=x_right_flush-1, y=0.48, yend=0.52), inherit.aes = FALSE, color="black", linewidth=0.5)+
  geom_text(aes(x=x_left_flush+32, y=0.50, label=paste0(round(dur_trans(dipt.flush), 1), " days")), inherit.aes = FALSE, color="black", size=4)+
  scale_x_continuous(breaks = seq(0, 365, by = 30), labels = function(x) day_to_DOY(x, 2019, start_month = pheno_start_month)) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank()
        )
windows()
p<- arrangeGrob(
  drop_panel,
  flush_panel,
  ncol = 2,
  top="Dipteryx oleifera leafing phenology",
  bottom="Day of Year (DOY)"
)
ggsave("timeseriesv2//figures//dipteryx_leafing_phenology.jpg", p, width = 12, height = 6, dpi = 300)



#to create a plot of dipteryx drop and flush we will
#pull the data points of the drop cycle
#parameters needed
amp <- summary(model.drop)$statistics["amp", "Mean"]
amp_flush <- summary(model.flush)$statistics["amp", "Mean"]
base <- summary(model.drop)$statistics["base", "Mean"]
base_flush <- summary(model.flush)$statistics["base", "Mean"]
Td <- summary(model.drop)$statistics["Td", "Mean"]
Td_flush <- summary(model.flush)$statistics["Td", "Mean"]
kd <- summary(model.drop)$statistics["kd", "Mean"]
kd_flush <- summary(model.flush)$statistics["kd", "Mean"]
midpoint<- (Td_flush+ Td)/2

df_drop<- data.frame(day=seq(min(drop_seg$day), round(midpoint, 0), by=1))
for (i in 1:length(df_drop$day)){
  df_drop$y_norm[i] <- (amp/(1+exp(kd*(df_drop$day[i]-Td))))+base
  df_drop$date[i]<-as.Date(paste(day_to_DOY(df_drop$day[i], 2019, start_month = 8), "2019"), format="%B %d %Y")
}

df_flush<- data.frame(day=seq(round(midpoint, 0), max(flush_seg$day), by=1))
for (i in 1:length(df_flush$day)){
  df_flush$y_norm[i] <- (amp_flush/(1+exp(-kd_flush*(df_flush$day[i]-Td_flush))))+base_flush
  df_flush$date[i] <- as.Date(paste(day_to_DOY(df_flush$day[i], 2019, start_month = 10), "2019"), format="%B %d %Y")
}


x_left <- Td - round(dur_trans(model.drop)/2, 0)
x_right <- Td + round(dur_trans(model.drop)/2, 0)
x_left_flush <- Td_flush - round(dur_trans(model.flush)/2, 0)
x_right_flush <- Td_flush + round(dur_trans(model.flush)/2, 0)
drop_panel <- ggplot(drop_seg[drop_seg$day <= round(midpoint, 0), ], aes(x = day, y = y_norm, color = tree_year)) +
  geom_point(alpha = 0.3) +
  geom_line(data = df_drop, aes(x = day, y = y_norm), color = "red", size = 1) +
  geom_vline(xintercept = Td, color = "black", linetype = "dashed") +
  geom_vline(xintercept = x_left, color = "black", linetype = "dashed") +
  geom_vline(xintercept = x_right, color = "black", linetype = "dashed") +
  geom_segment(aes(x=x_left+1, xend=x_right-1, y=0.5, yend=0.5), inherit.aes = FALSE, color="black", linewidth=0.5)+
  geom_segment(aes(x=x_left+1, xend=x_left+1, y=0.48, yend=0.52), inherit.aes = FALSE, color="black", linewidth=0.5)+
  geom_segment(aes(x=x_right-1, xend=x_right-1, y=0.48, yend=0.52), inherit.aes = FALSE, color="black", linewidth=0.5)+
  geom_text(aes(x=x_left-15, y=0.50, label=paste0(round(dur_trans(model.drop), 1), " days")), inherit.aes = FALSE, color="black", size=4)+
  scale_x_continuous(breaks = c(seq(0, 365, by = 30)), labels = function(x) day_to_DOY(x, 2019, start_month = 8)) +
  ylab("Leaf Cover") +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
  )

flush_panel<-ggplot(flush_seg[flush_seg$day >= round(midpoint, 0), ], aes(x = day, y = y_norm, color=tree_year)) +
  geom_point(alpha=0.3)+
  geom_line(data=df_flush, aes(x=day, y=y_norm), color="blue", size=1)+
  geom_vline(xintercept=Td_flush, color="black", linetype="dashed")+
  geom_vline(xintercept=x_left_flush, color="black", linetype="dashed")+
  geom_vline(xintercept=x_right_flush, color="black", linetype="dashed")+
  geom_segment(aes(x=x_left_flush+1, xend=x_right_flush-1, y=0.5, yend=0.5), inherit.aes = FALSE, color="black", linewidth=0.5)+
  geom_segment(aes(x=x_left_flush+1, xend=x_left_flush+1, y=0.48, yend=0.52), inherit.aes = FALSE, color="black", linewidth=0.5)+
  geom_segment(aes(x=x_right_flush-1, xend=x_right_flush-1, y=0.48, yend=0.52), inherit.aes = FALSE, color="black", linewidth=0.5)+
  geom_text(aes(x=x_left_flush+32, y=0.50, label=paste0(round(dur_trans(model.flush), 1), " days")), inherit.aes = FALSE, color="black", size=4)+
  scale_x_continuous(breaks = seq(0, 365, by = 30), labels = function(x) day_to_DOY(x, 2019, start_month = 10)) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank()
        )

p<- arrangeGrob(
  drop_panel,
  flush_panel,
  ncol = 2,
  top="Cavallinesia platanifolia leafing phenology",
  bottom="Day of Year (DOY)"
)
ggsave("timeseriesv2//figures//cavallinesia_leafing_phenology.jpg", p, width = 12, height = 6, dpi = 300)
