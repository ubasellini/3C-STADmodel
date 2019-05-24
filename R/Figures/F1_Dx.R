## clean workspace
rm(list=ls(all=TRUE))

## load packages
library(msm)
library(RColorBrewer)
library(ggplot2)
library(ggridges)
library(demography)
library(tidyverse)
library(viridis)

## load country-specific data
setwd("~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Data")
cou <- "CHE"      ## "CHE" or "SWE"
sex <- "Males"    ## "Females" or "Males"

## load country-specific data
name <- paste(cou,"_HMD.Rdata",sep="")
load(name)

## select country, ages and years
years <- 1950:2016
x <- ages <- 0:110
n <- length(years)
m <- length(ages)
y <- round(seq(years[1],years[n],by = 10))
y <- c(round(seq(years[1],years[n],by = 10)),2016)
FittingData <- extract.years(FullData, years=years)
FittingData <- extract.ages(FittingData, ages=ages)
if (sex=="Males"){
  E <- FittingData$pop$male
  MX <- FittingData$rate$male
  Z <- E*MX
  LHAZact <- log(Z/E)
}else if(sex=="Females"){
  E <- FittingData$pop$female
  MX <- FittingData$rate$female
  Z <- E*MX
  LHAZact <- log(Z/E)
}

## set seed for reproducibility
set.seed(2019)
sim.tot <- 100000

## simulate life-times from the observed mortality pattern
## (useful to consruct age-at-death distributions)
df <- data.frame(x=NA,y=rep(y,each=sim.tot))
i <- 1
for (i in 1:length(y)){
  j <- which(years==y[i])
  z <- as.vector(Z[,j])
  e <- as.vector(E[,j])
  sim.lifes <- rpexp(n=sim.tot,rate=z/e,x)
  df[df$y==y[i],]$x <- sim.lifes
}
df$y <- as.factor(df$y)

## Z theme
z_theme <- function() {
  library(RColorBrewer)
  # Generate the colors for the chart procedurally with RColorBrewer
  palette <- brewer.pal("Greys", n=9)
  color.background = "#FFFFFF"
  color.grid.major = palette[3]
  color.axis.text = palette[7]
  color.axis.title = palette[7]
  color.title = palette[8]
  # Begin construction of chart
  theme_bw(base_size=9) +
    # Set the entire chart region to a light gray color
    theme(panel.background=element_rect(fill=color.background, color=color.background)) +
    theme(plot.background=element_rect(fill=color.background, color=color.background)) +
    theme(panel.border=element_rect(color=color.background)) +
    # Format the grid
    theme(panel.grid.major=element_line(color=color.grid.major,size=.25)) +
    theme(panel.grid.minor=element_blank()) +
    theme(axis.ticks=element_blank()) +
    # Format the legend, but hide by default
    theme(legend.position="none") +
    theme(legend.background = element_rect(fill=color.background)) +
    theme(legend.text = element_text(size=7,color=color.axis.title)) +
    # Set title and axis labels, and format these and tick marks
    theme(plot.title=element_text(color=color.title, size=20, vjust=1.25)) +
    theme(axis.text.x=element_text(size=14,color=color.axis.text)) +
    theme(axis.text.y=element_text(size=14,color=color.axis.text)) +
    theme(axis.title.x=element_text(size=16,color=color.axis.title, vjust=0)) +
    theme(axis.title.y=element_text(size=16,color=color.axis.title, vjust=1.25))
}

## Choose colors
palette1 <- brewer.pal("YlGnBu", n=9)
palette2 <- brewer.pal("YlOrBr", n=9)
my.cols <- c(palette1[6],palette2[5])
my.cols <- adjustcolor(my.cols, alpha=0.95)
# display.brewer.pal(9,"YlOrBr")
scale.fig <- 0.98

ggplot(df, aes(x = x, y = y, fill = ..x..)) +
  geom_density_ridges_gradient(scale = scale.fig,quantile_lines = F,from=0,to=110,
                               calc_ecdf = T,bandwidth = 3,quantiles = c(0.25, 0.75)) +
  scale_fill_viridis(option = "C", guide = "none") +
  labs(title = 'Age-at-death distributions') + z_theme()
F1 <- ggplot(df, aes(x = x, y = y, fill = ..x..)) +
  geom_density_ridges_gradient(scale = scale.fig,quantile_lines = F,from=0,to=110,
                               calc_ecdf = T,bandwidth = 3,quantiles = c(0.25, 0.75)) +
  scale_fill_viridis(option = "C", guide = "none") +
  labs(title = 'Age-at-death distributions') + z_theme()
ingredients <- ggplot_build(F1) %>% purrr::pluck("data", 1)
density_lines <- ingredients %>%
  group_by(group) %>% filter(density == max(density)) %>% ungroup()

## compute IQR
bla <- ingredients %>% group_by(group,quantile) %>%
  mutate(Qmin=min(x),Qmax=max(x)) 
bla2 <- bla %>% select(Qmin,Qmax,y) %>% distinct(Qmin,Qmax,y)
bla3 <- bla2 %>% ungroup() %>% mutate(Qmin2=lag(Qmax),Qmax2=lead(Qmin)) %>% 
  filter(quantile==2) %>%
  mutate(QminFIN = (Qmin+Qmin2)/2,QmaxFIN = (Qmax+Qmax2)/2,
         IQR = round(QmaxFIN-QminFIN,1),z=group+0.35)
bla3$x <- density_lines$x-4

F1 + geom_segment(data = density_lines, 
      aes(x = x, y = ymin, xend = x, yend = ymin+density*scale*iscale)) +
    geom_text(data=bla3,inherit.aes = F, 
            aes(x=x,y=z,label=sprintf("%1.1f", IQR)), 
            position=position_nudge(y=-0.1), colour="black", size=3.5)

ggplot(df, aes(x=x, y=y, fill=factor(..quantile..))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = T,bandwidth = 3,
                      scale = scale.fig,quantiles = c(0.25, 0.75),from=0,to=110,
                      quantile_lines = F) +
  scale_fill_manual(
    name = "", values = c(my.cols[1], my.cols[2], my.cols[1])) + 
  z_theme() + scale_x_continuous(breaks=seq(0,120,20)) +
  labs(x="Age",
       y="Year") +
  scale_y_discrete(expand = expand_scale(add = c(0.1, 1.03))) + 
  geom_segment(data = density_lines,inherit.aes = F, 
                 aes(x = x, y = ymin, xend = x, yend = ymin+density*iscale*scale),
                 linetype="dashed") +
  geom_text(data=bla3,inherit.aes = F,
            aes(x=x,y=z,label=sprintf("%1.1f", IQR)),
            position=position_nudge(y=-0.1), colour="black", size=3.5)

setwd("~/Documents/Demography/Work/STADall/Github/3C-STADmodel/Chapter/Figures")
ggsave("F1.pdf",width = 10, height = 6)
