library('RBi.helpers')
library('data.table')
library('ggthemr')
library('devtools')
library('scales')
library('reshape2')

ggthemr('fresh')

code_dir <- path.expand("~/code/ebola_lofa/")

## 1. case numbers
inc_filename <- "lofa_incidence.rds"
incidence <- readRDS(paste(code_dir, "data", inc_filename, sep = "/"))

weekly_admissions <- incidence[["admissions"]]
weekly_admissions <- weekly_admissions[date >= "2014-06-01" & date < as.Date("2014-11-01") - 7]
weekly_admissions[, type := factor("hospital admissions", levels = c("hospital admissions", "community deaths", "admission delay", "People reached through HP activities", "Time line", "R0 without isolation", "Proportion seeking healthcare"))]
weekly_admissions[is.na(classification), classification := "unknown"]

## 2. community burials

community_deaths <- incidence[["deaths"]]
community_deaths <- community_deaths[date >= "2014-06-01" & date < as.Date("2014-11-01") - 7]
community_deaths[, type := factor("community deaths", levels = c("hospital admissions", "community deaths", "admission delay", "People reached through HP activities", "Time line", "R0 without isolation", "Proportion seeking healthcare"))]
community_deaths[is.na(classification), classification := "unknown"]

## 3. admission delay

adm_filename <- "lofa_admission_delays.rds"
admission_delay <- readRDS(paste(code_dir, "data", adm_filename, sep = "/"))
admission_delay <- admission_delay[date >= as.Date("2014-06-01") - 7]
admission_delay <- admission_delay[date < as.Date("2014-11-01") - 7]
admission_delay[, type := factor("admission delay", levels = c("hospital admissions", "community deaths", "admission delay", "People reached through HP activities", "Time line", "R0 without isolation", "Proportion seeking healthcare"))]
admission_delay[, date_next := date + 7]
admission_delay[nrow(admission_delay), min.50 := NA]
admission_delay[nrow(admission_delay), max.50 := NA]
admission_delay[nrow(admission_delay), min.90 := NA]
admission_delay[nrow(admission_delay), max.90 := NA]

## 4. HP activities

hp_filename <- "hp_activities.rds"
reached <- readRDS(paste(code_dir, "data", hp_filename, sep = "/"))
reached <- reached[date < as.Date("2014-11-01") - 7]
reached[, type := factor("People reached through HP activities", levels = c("hospital admissions", "community deaths", "admission delay", "People reached through HP activities", "Time line", "R0 without isolation", "Proportion seeking healthcare"))]

## 5. Timeline

timeline <- data.table(rbind(
    c(date = "2014-04-12", event = "ETC established with 10 beds"),
    c(date = "2014-08-06", event = "State of emergency declared by the government"),
    c(date = "2014-08-10", event = "Expanded ETC (100 beds) opened"), 
    c(date = "2014-09-15", event = "Additional support by MOH"), 
    c(date = "2014-07-15", event = "WHO starts training community volunteers"), 
    c(date = "2014-04-15", event = "Initial GCHV training"), 
    c(date = "2014-09-15", event = "Health structure improved by MOH"), 
    c(date = "2014-08-15", event = "Social mobilisation initiated"), 
    c(date = "2014-09-15", event = "Social mobilisation expanded"), 
    c(date = "2014-07-01", event = "Hand-washing buckets distributed")
))
timeline <- timeline[date >= "2014-06-01"]
timeline[, date := as.Date(date)]
setkey(timeline, date)
timeline[, type := factor("Time line", levels = c("hospital admissions", "community deaths", "admission delay", "People reached through HP activities", "Time line", "R0 without isolation", "Proportion seeking healthcare"))]
timeline[, y := 0]
timeline[, yend := (nrow(timeline):1) + 1]

## 6. Community reproduction number
traces <- readRDS("lofa_parallel_traces.rds")
cR0 <- plot_libbi(traces, model = "~/code/ebola_lofa/ebola_lofa_fit.bi",
                  date.origin = as.Date("2014-06-02") - 7, date.unit = "week",
                  states = "R0", params = NULL, noises = NULL,
                  trend = "median", step = TRUE, quantile = c(0.5, 0.9))
cR0$data$states[, type := factor("R0 without isolation", levels = c("hospital admissions", "community deaths", "admission delay", "People reached through HP activities", "Time line", "R0 without isolation", "Proportion seeking healthcare"))]
cR0$data$states <- cR0$data$states[-nrow(cR0$data$states)]
cR0$data$states <- rbind(cR0$data$states, cR0$data$states[nrow(cR0$data$states)])
cR0$data$states[nrow(cR0$data$states), time_next := NA]
cR0$data$states[nrow(cR0$data$states), time := time + 7]

## 7. Healthcare-seeking behaviour
H <- plot_libbi(traces, model = "~/code/ebola_lofa/ebola_lofa_fit.bi",
                date.origin = as.Date("2014-06-02") - 7, date.unit = "week",
                states = "H", params = NULL, noises = NULL,
                trend = "median", step = TRUE, quantile = c(0.5, 0.75))
H$data$states[, type := factor("Proportion seeking healthcare", levels = c("hospital admissions", "community deaths", "admission delay", "People reached through HP activities", "Time line", "R0 without isolation", "Proportion seeking healthcare"))]
H$data$states <- H$data$states[-nrow(H$data$states)]
H$data$states <- rbind(H$data$states, H$data$states[nrow(H$data$states)])
H$data$states[nrow(H$data$states), time_next := NA]
H$data$states[nrow(H$data$states), time := time + 7]
H$data$states[, value := value * 100]
H$data$states[, min.1 := min.1 * 100]
H$data$states[, min.2 := min.2 * 100]
H$data$states[, max.1 := max.1 * 100]
H$data$states[, max.2 := max.2 * 100]

R0_line <- data.table(median = 1, type = factor("R0 without isolation", levels = c("hospital admissions", "community deaths", "admission delay", "People reached through HP activities", "Time line", "R0 without isolation", "Proportion seeking healthcare")))

p <- ggplot()
p <- p + geom_bar(stat = "identity", data = weekly_admissions, mapping = aes(x = date, y = admissions, fill = factor(classification, c("confirmed", "probable", "not a case", "unknown")), order = factor(classification, c("confirmed", "probable", "not a case", "unknown"))))
p <- p + geom_bar(stat = "identity", data = community_deaths, mapping = aes(x = date, y = deaths, fill = factor(classification, c("confirmed", "probable", "not a case", "unknown")), order = factor(classification, c("confirmed", "probable", "not a case", "unknown"))))
p <- p + geom_step(data = admission_delay, mapping = aes(x = date, y = median), color = "black", lwd = 1.5)
p <- p + geom_rect(data = admission_delay, mapping = aes(xmin = date, xmax = date_next, ymin = min.50, ymax = max.50), alpha = 0.3)
p <- p + geom_bar(stat = "identity", data = reached, mapping = aes(date, y = total), color = "black")
p <- p + geom_segment(data = timeline, mapping = aes(x = date, xend = date, y = y, yend = yend))
p <- p + geom_text(data = timeline, mapping = aes(x = date + 2, y = yend, label = event), hjust = 0, size = 3, color = "black")
p <- p + geom_step(data = cR0$data$states, mapping = aes(x = time, y = value), lwd = 1.5)
p <- p + geom_rect(data = cR0$data$states, mapping = aes(xmin = time, xmax = time_next, ymin = min.1, ymax = max.1), alpha = 0.3)
p <- p + geom_hline(data = R0_line, mapping = aes(yintercept = median), color = "black", linetype = "dashed")
p <- p + geom_step(data = H$data$states, mapping = aes(x = time, y = value), lwd = 1.5)
p <- p + geom_rect(data = H$data$states, mapping = aes(xmin = time, xmax = time_next, ymin = min.1, ymax = max.1), alpha = 0.3)
p <- p + scale_fill_brewer("", palette = "Set1")
p <- p + facet_grid(type~., scales = "free")
p <- p + theme(legend.position = "top", strip.text.y = element_blank())
p <- p + scale_y_continuous(breaks = pretty_breaks(),
                            expression(atop("Healthcare            Community              Timeline of             # of people                 Days to             Community              Hospital      ",
                                        "    seeking (in %)   reproduction number         events              reached by HP          admission               deaths               admissions         ")))
p <- p + scale_x_date("", date_breaks = "1 month", date_minor_breaks = "1 week", date_labels = "%e %b")
p <- p + theme(text = element_text(size = 14))
p <- p + expand_limits(y = 0)
ggsave("~/timeline.pdf", p, width = 8.6, height = 13)
ggsave("~/timeline.svg", p, width = 8.6, height = 13)
ggsave("~/timeline.png", p, width = 8.6, height = 13)

## write data files
weekly_admissions[, date := NULL]
community_deaths[, date := NULL]
admission_delay[, date := NULL]
reached[, date := NULL]
write.table(weekly_admissions, quote = FALSE, sep = ",", row.names = FALSE, file = "etc_admissions.csv")
write.table(community_deaths, quote = FALSE, sep = ",", row.names = FALSE, file = "community_deaths.csv")
write.table(admission_delay, quote = FALSE, sep = ",", row.names = FALSE, file = "admission_delay.csv")
write.table(reached, quote = FALSE, sep = ",", row.names = FALSE, file = "reached.csv")

## plot just hospital admission
p <- ggplot(weekly_admissions, aes(x=date, y=admissions))+geom_bar(stat="identity", fill="black")+ylab("New hospital admissions")+scale_x_date(date_labels="%b %Y")+xlab("")+theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot("lofa_admission_data.pdf", p)
save_plot("lofa_admission_data.png", p)
