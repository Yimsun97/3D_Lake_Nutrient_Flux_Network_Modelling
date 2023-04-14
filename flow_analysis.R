library(tidyverse)
library(network)
library(enaR)
library(plyr)
library(readxl)
library(ggpubr)
library(patchwork)
library(extrafont)

loadfonts(device = "win")

path <- 'D:/Files/Projects/NonLinearity/'

correct_ben_flux <- function(flow.mat, varname){
  flow_ <- flow.mat
  if (varname == 'N'){
    if (flow_['NH4', 'SED_N'] > 0){
      flow_['NH4', 'SED_N'] <- flow_['SED_N', 'NH4'] + flow_['NH4', 'SED_N']
    }
    if (flow.mat['NO3', 'SED_N'] > 0){
      flow_['NO3', 'SED_N'] <- flow_['SED_N', 'NO3'] + flow_['NO3', 'SED_N']
    }
  } else {
    if (flow_['PO4', 'SED_P'] > 0){
      flow_['PO4', 'SED_P'] <- flow_['SED_P', 'PO4'] + flow_['PO4', 'SED_P']
    }
  }
}

extract_sc_stats <- function(scenario, varname, time, freq = "m"){
  path_ <- paste(path, "network tables/", sep = "")
  flow.mat <- read.csv(paste(path_, scenario, "/flux_", varname, 
                             "_", freq, time-1, ".csv", 
                             sep = ""), header=TRUE, row.names=1)
  input.mat <- read.csv(paste(path_, scenario, "/input_", varname,
                              "_", freq, time-1, ".csv", 
                              sep = ""), header=TRUE, row.names=1)
  output.mat <- read.csv(paste(path_, scenario, "/output_", varname, 
                               "_", freq, time-1, ".csv",  
                               sep = ""), header=TRUE, row.names=1)
  storage.mat <- read.csv(paste(path_, scenario, "/mass_", varname, 
                                "_", freq, time-1, ".csv", 
                                sep = ""), header=TRUE, row.names=1)
  flow <- as.matrix(flow.mat)
  input <- c(input.mat)[[1]]
  export <- c(output.mat)[[1]]
  storage <- c(storage.mat)[[1]]

  lake.model <- enaR::pack(
    flow = abs(flow),
    input = abs(input),
    export = abs(export),
    storage = abs(storage),
  )
  
  lake.flow <- enaFlow(lake.model, balance.override = TRUE)
  lake.storage <- enaStorage(lake.model, balance.override = TRUE)
  
  model.res <- list(model = lake.model,
                    TST = lake.flow$ns[,'TST'][[1]],
                    APL = lake.flow$ns[,'APL'][[1]],
                    FCI = lake.flow$ns[,'FCI'][[1]],
                    TCF = lake.flow$ns[,'mode2.F'][[1]],
                    SCI = lake.storage$ns[,'CIS'][[1]],
                    TSS = lake.storage$ns[,'TSS'][[1]])
  return(model.res)
}


extract_all_stats <- function(scs, varnames, times, freq = "m"){
  all.res <- list()
  i <- 1
  for (sc in scs) {
    for (varname in varnames) {
      for (time in times) {
        lake.res <- extract_sc_stats(sc, varname, time, freq = freq)
        lake.stats <- as.data.frame(lake.res[2:length(lake.res)])
        lake.stats <- lake.stats %>% mutate(
          sc = sc, 
          varname = varname,
          time = time,
          input = sum(lake.res[[1]]%v%"input"),
          output = sum(lake.res[[1]]%v%"output")
        )
        
        all.res[[i]] <- lake.stats
        print(i)
        print(sc)
        print(varname)
        print(time)
        i = i + 1
      }
    }
  }
  all.res.df <- ldply(all.res, rbind)
  all.res.tb <- as_tibble(all.res.df)
  return(all.res.tb)
}

lake.res <- extract_sc_stats(scenario = "S06", 
                             varname = "P", 
                             time = 4, 
                             freq = "m")

wqs <- read_excel(paste(path, '/WQ tables/wq_max+min_TN.xlsx', sep = ""))
colnames(wqs)[2:12] <- c("S00", "S01", "S02", "S03", "S04", "S05",
                         "S06", "S07", "S08", "S09", "S10")
wqs <- wqs %>% mutate(time = ceiling(TIME - 2193))

wqs.long <- wqs %>% pivot_longer(S00:S10, names_to = "sc", values_to = "wq")

scs <- c("S00", "S01", "S02", "S03", "S04", "S05",
         "S06", "S07", "S08", "S09", "S10")
scs.vis <- c("S00", "S01", "S03", "S06")
varnames <- c("N", "P")
freq = 'd'
if (freq == 'm') {
  times <- 1:12
} else {
  times <- 1:364
}

all.res <- extract_all_stats(scs, varnames, times, freq = freq)

merged.res <- left_join(all.res, wqs.long, by = c("sc", "time"))
month.map <- data.frame(month = 1:12, season = c(4, 1, 1, 1, 2, 2, 2, 
                                                 3, 3, 3, 4, 4))
month.map <- data.frame(month = 1:12, season = c(0, 0, 0, 0, 1, 1,
                                                 1, 1, 1, 1, 0, 0))
seasons.vis <- c(1, 3)

merged.res <- left_join(merged.res, month.map, by = "month")

all.res.N <- merged.res %>% filter(varname == 'N')
all.res.P <- merged.res %>% filter(varname == 'P')

all.res.N <- all.res.N %>% mutate(rTST = TST / max(all.res.N[all.res.N$sc == "S00", "TST"]))
all.res.P <- all.res.P %>% mutate(rTST = TST / max(all.res.P[all.res.P$sc == "S00", "TST"]))

all.res.N <- all.res.N %>% mutate(rTSS = TSS / max(all.res.N[all.res.N$sc == "S00", "TSS"]))
all.res.P <- all.res.P %>% mutate(rTSS = TSS / max(all.res.P[all.res.P$sc == "S00", "TSS"]))


N.res.long <- all.res.N %>% select(FCI, rTST, rTSS, season, sc) %>% 
  pivot_longer(FCI:rTSS, names_to = "index", values_to = "value")
P.res.long <- all.res.P %>% select(FCI, rTST, rTSS, season, sc) %>% 
  pivot_longer(FCI:rTSS, names_to = "index", values_to = "value")

N.res.long.month <- all.res.N %>% select(FCI, rTST, rTSS, month, sc) %>% 
  pivot_longer(FCI:rTSS, names_to = "index", values_to = "value")
P.res.long.month <- all.res.P %>% select(FCI, rTST, rTSS, month, sc) %>% 
  pivot_longer(FCI:rTSS, names_to = "index", values_to = "value")

# for box plots
np.res.month <- rbind(N.res.long.month %>% mutate(varname = "N"),
                      P.res.long.month %>% mutate(varname = "P"))
write.csv(np.res.month, file = "NetworkIndicators.csv", row.names = FALSE)

scs.vis <- c("S00", "S01", "S03", "S06", "S07", "S09")
labels <- c("S0", "S1", "S3", "S6", "S7", "S9")
colors <- c('#000000', '#1f77b4',  '#ff7f0e', '#2ca02c', '#d62728', 
            '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', 
            '#17becf')
colors.vis <- colors[c(1, 2, 4, 7, 8, 10)]
month.vis <- c(5, 6, 7, 8, 9, 10)
linewidth <- 1 / (ggplot2::.pt * 72.27/96)
varnames <- c(
  `N` = "N",
  `P` = "P"
)
font.size <- 7
(p.var.fci <- ggboxplot(np.res.month %>% filter(sc %in% scs.vis,
                                                month %in% month.vis,
                                                index == 'FCI'), 
                        x = "sc", y = "value", 
                        size = 0.3,
                        width = 0.3, 
                        notch = FALSE, 
                        color = "sc", 
                        # palette = "jco",
                        add = "mean_se",
                        error.plot = "errorbar",
                        outlier.shape = NA,
                        # facet.by = c("varname"),
                        # scales = "free",
                        add.params = list(size = 0.3, color = "black"),
                        short.panel.labs = FALSE) +
    stat_compare_means(method = "wilcox.test", 
                       label = "p.signif",
                       ref.group = "S00",
                       family = "Arial",
                       color = 'red',
                       size = font.size * 0.3528
                       ) +
    facet_grid(. ~ varname, labeller = as_labeller(varnames))+
    scale_y_continuous(name = "FCI")+
    scale_color_manual(values = colors.vis)+
    scale_x_discrete(labels = labels)+
    theme_bw()+
    theme(strip.background = element_rect(fill = NA, color = "black",
                                          linewidth = linewidth),
          # strip.text.x = element_blank(),
          legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          text=element_text(size=font.size, family="Arial", color = "black"),
          axis.text = element_text(size=font.size, family="Arial", color = "black"),
          axis.title.x = element_blank(),
          panel.background = element_rect(fill=NA),
          plot.background = element_rect(fill=NA, color=NA),
          panel.border = element_rect(fill=NA, colour = "black", 
                                      linewidth = linewidth))
    )
(p.var.tss <- ggboxplot(np.res.month %>% filter(sc %in% scs.vis,
                                                month %in% month.vis,
                                                index == 'rTSS'), 
                        x = "sc", y = "value", 
                        size = 0.3,
                        width = 0.3, 
                        notch = FALSE, 
                        color = "sc", 
                        # palette = "jco",
                        add = "mean_se", 
                        error.plot = "errorbar",
                        outlier.shape = NA,
                        facet.by = "varname",
                        scales = "free",
                        add.params = list(size = 0.3, color = "black"),
                        short.panel.labs = FALSE) +
    stat_compare_means(method = "wilcox.test", 
                       label = "p.signif",
                       ref.group = "S00",
                       family = "Arial",
                       color = 'red',
                       size = font.size * 0.3528) +
    scale_y_continuous(name = "rTSS")+
    scale_x_discrete(labels = labels)+
    scale_color_manual(values = colors.vis)+
    theme_bw()+
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          text=element_text(size=font.size, family="Arial", color = "black"),
          axis.title.x = element_blank(),
          axis.text = element_text(size=font.size, family="Arial", color = "black"),
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          panel.border = element_rect(fill=NA, colour = "black", 
                                      linewidth = linewidth))
)

(p.var.tst <- ggboxplot(np.res.month %>% filter(sc %in% scs.vis,
                                                month %in% month.vis,
                                                index == 'rTST'), 
                        x = "sc", y = "value", 
                        size = 0.3,
                        width = 0.3, 
                        notch = FALSE, 
                        color = "sc", 
                        # palette = "jco",
                        add = "mean_se",
                        error.plot = "errorbar",
                        outlier.shape = NA,
                        # facet.by = c("varname"),
                        # scales = "free",
                        add.params = list(size = 0.3, color = "black"),
                        short.panel.labs = FALSE) +
    stat_compare_means(method = "wilcox.test", 
                       label = "p.signif",
                       ref.group = "S00",
                       family = "Arial",
                       color = 'red',
                       size = font.size * 0.3528
    ) +
    facet_grid(. ~ varname, labeller = as_labeller(varnames))+
    scale_y_continuous(name = "rTST")+
    scale_color_manual(values = colors.vis)+
    scale_x_discrete(name = "Scenario", 
                     labels = labels)+
    theme_bw()+
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          text=element_text(size=font.size, family="Arial", color = "black"),
          axis.text = element_text(size=font.size, family="Arial", color = "black"),
          # axis.title.x = element_blank(),
          panel.background = element_rect(fill=NA),
          plot.background = element_rect(fill=NA, color=NA),
          panel.border = element_rect(fill=NA, colour = "black", 
                                      linewidth = linewidth))
)

(p.var.pub <- p.var.fci / p.var.tss / p.var.tst)

ggsave("figures/fs_seasons_NP.svg", 
       width = 14, height = 18, units = "cm",
       dpi = 500, bg='transparent', 
       plot = p.var.pub)
ggsave("figures/fs_seasons_NP.png", 
       width = 14, height = 18, units = "cm",
       dpi = 500, bg='transparent', 
       plot = p.var.pub)


# for box plots in the SI
np.res.month <- rbind(N.res.long.month %>% mutate(varname = "N"),
                      P.res.long.month %>% mutate(varname = "P"))
scs.vis <- c("S00", "S02", "S04", "S05", "S08", "S10")
labels <- c("S0", "S2", "S4", "S5", "S8", "S10")
colors.vis <- colors[c(1, 3, 5, 6, 9, 11)]
month.vis <- c(5, 6, 7, 8, 9, 10)
linewidth <- 1 / (ggplot2::.pt * 72.27/96)
varnames <- c(
  `N` = "N",
  `P` = "P"
)
font.size <- 7
(p.var.fci <- ggboxplot(np.res.month %>% filter(sc %in% scs.vis,
                                                month %in% month.vis,
                                                index == 'FCI'), 
                        x = "sc", y = "value", 
                        size = 0.3,
                        width = 0.3, 
                        notch = FALSE, 
                        color = "sc", 
                        # palette = "jco",
                        add = "mean_se",
                        error.plot = "errorbar",
                        outlier.shape = NA,
                        # facet.by = c("varname"),
                        # scales = "free",
                        add.params = list(size = 0.3, color = "black"),
                        short.panel.labs = FALSE) +
    stat_compare_means(method = "wilcox.test", 
                       label = "p.signif",
                       ref.group = "S00",
                       family = "Arial",
                       color = 'red',
                       size = font.size * 0.3528
    ) +
    facet_grid(. ~ varname, labeller = as_labeller(varnames))+
    scale_y_continuous(name = "FCI")+
    scale_color_manual(values = colors.vis)+
    scale_x_discrete(labels = labels)+
    theme_bw()+
    theme(strip.background = element_rect(fill = NA, color = "black",
                                          linewidth = linewidth),
          # strip.text.x = element_blank(),
          legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          text=element_text(size=font.size, family="Arial", color = "black"),
          axis.text = element_text(size=font.size, family="Arial", color = "black"),
          axis.title.x = element_blank(),
          panel.background = element_rect(fill=NA),
          plot.background = element_rect(fill=NA, color=NA),
          panel.border = element_rect(fill=NA, colour = "black", 
                                      linewidth = linewidth))
)
(p.var.tss <- ggboxplot(np.res.month %>% filter(sc %in% scs.vis,
                                                month %in% month.vis,
                                                index == 'rTSS'), 
                        x = "sc", y = "value", 
                        size = 0.3,
                        width = 0.3, 
                        notch = FALSE, 
                        color = "sc", 
                        # palette = "jco",
                        add = "mean_se", 
                        error.plot = "errorbar",
                        outlier.shape = NA,
                        facet.by = "varname",
                        scales = "free",
                        add.params = list(size = 0.3, color = "black"),
                        short.panel.labs = FALSE) +
    stat_compare_means(method = "wilcox.test", 
                       label = "p.signif",
                       ref.group = "S00",
                       family = "Arial",
                       color = 'red',
                       size = font.size * 0.3528) +
    scale_y_continuous(name = "rTSS")+
    scale_x_discrete(labels = labels)+
    scale_color_manual(values = colors.vis)+
    theme_bw()+
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          text=element_text(size=font.size, family="Arial", color = "black"),
          axis.title.x = element_blank(),
          axis.text = element_text(size=font.size, family="Arial", color = "black"),
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          panel.border = element_rect(fill=NA, colour = "black", 
                                      linewidth = linewidth))
)

(p.var.tst <- ggboxplot(np.res.month %>% filter(sc %in% scs.vis,
                                                month %in% month.vis,
                                                index == 'rTST'), 
                        x = "sc", y = "value", 
                        size = 0.3,
                        width = 0.3, 
                        notch = FALSE, 
                        color = "sc", 
                        # palette = "jco",
                        add = "mean_se",
                        error.plot = "errorbar",
                        outlier.shape = NA,
                        # facet.by = c("varname"),
                        add.params = list(size = 0.3, color = "black"),
                        short.panel.labs = FALSE) +
    stat_compare_means(method = "wilcox.test", 
                       label = "p.signif",
                       ref.group = "S00",
                       family = "Arial",
                       color = 'red',
                       size = font.size * 0.3528
    ) +
    facet_grid(. ~ varname, labeller = as_labeller(varnames))+
    scale_y_continuous(name = "rTST")+
    scale_color_manual(values = colors.vis)+
    scale_x_discrete(name = "Scenario", 
                     labels = labels)+
    theme_bw()+
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          text=element_text(size=font.size, family="Arial", color = "black"),
          axis.text = element_text(size=font.size, family="Arial", color = "black"),
          # axis.title.x = element_blank(),
          panel.background = element_rect(fill=NA),
          plot.background = element_rect(fill=NA, color=NA),
          panel.border = element_rect(fill=NA, colour = "black", 
                                      linewidth = linewidth))
)

(p.var.supp <- p.var.fci / p.var.tss / p.var.tst)


ggsave("figures/fs_seasons_NP_supp.png", 
       width = 14, height = 18, units = "cm",
       dpi = 500, bg='transparent', 
       plot = p.var.supp)

# some statistics
mean((np.res.month %>% filter(sc == "S06",
                              month %in% month.vis, varname == 'N',
                              index == 'FCI'))$value)

sd((np.res.month %>% filter(sc == "S06",
                            month %in% month.vis, varname == 'N',
                            index == 'FCI'))$value)

mean((np.res.month %>% filter(sc == "S00",
                              month %in% month.vis, varname == 'N',
                              index == 'rTST'))$value)

sd((np.res.month %>% filter(sc == "S00",
                            month %in% month.vis, varname == 'N',
                            index == 'rTST'))$value)

mean((np.res.month %>% filter(sc == "S10",
                              month %in% month.vis, varname == 'N',
                              index == 'rTSS'))$value)

mean((np.res.month %>% filter(sc == "S10",
                              month %in% month.vis, varname == 'P',
                              index == 'rTSS'))$value)
