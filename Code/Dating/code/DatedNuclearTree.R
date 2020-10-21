
library(devtools)
library(ggplot2)
library(treeio)
library(ggtree)
library(ape)
library(dplyr)
library(deeptime)
library(ggthemes)

phy <- read.beast("mcc.trees")
phy
class(phy)

#in ggtree, x axis starts where your tree ends (here, left side near the root). So if the oldest node is 137 mya, the 0 of coordinate system will start from there#
#you would want to leave the white space before tree starts 
#so your tree starts from 137 (max egde length) but you would want the space till say 150, set lower limit of x axis (xlim) as -13 (~150-137)
#you want the upper limit of x axis (xlim) till some point after 137 so that your labels fit#
#so we will put the upper limit of x-axis as 167 (~137+52)#

#for y axis, the 0 is probably where the lowest branch is#
#so we will start y axis a little below the tree (-3.5 here) and
#end it a little above where there tree end#
#so u give the upperlimit of y axis as 59+1.5 (59 is the n_nodes or no. of nodes#

#epoch names and details are coming from deeptime
epochs2 <- epochs

#u cannot simply plot the absolute min and max age of epochs because of where the x-axis of ggtree starts
#so u need to substract from the max length of edge#
epochs2$max_age <- max(phy@phylo$edge.length) - epochs2$max_age
epochs2$min_age <- max(phy@phylo$edge.length) - epochs2$min_age

#make vector of tips you want to colour in green#
taxa_to_color<-c("Sixalix_farinosa", "Lomelosia_albocinta", "Pycnocomon_canus",   
                 "Pteracephalodes_hookeri", "Bassecoia_siamensis", 
                 "Dipsacus_fullomum", "Cephalaria_gigantea", 
                 "Succisa_pratensis", "Knautia_macedonica", 
                 "Triplostegia_grandulifera", "Valeriana_microphylla", 
                 "Fedia_cornucopiae","Plectritis_congesta", "Valeriana_lobata",
                 "Centranthus_ruber", "Nardostachys_jatamansi", 
                 "Patrina_triloba","Morina_longifolia", 
                 "Cryptothladia_kokonorica", "Acanthocalyx_nepalensis_ssp_delevaya",                  "Triosteum_himalayanum" )


# following are calibration nodes
nodes.to.label <- c(87,69, 103, 64, 97, 94, 93, 109, 117, 110, 62)

# This is replcating NA as many times as sum of nodes (59) and no. of tip labels(60)
#so it is making a vector of 119 NA
node.label.cols <- rep(NA, phy@phylo$Nnode + length(phy@phylo$tip.label))

node.label.cols[nodes.to.label] <- 'Black' 

nodes.text <- rep(NA, phy@phylo$Nnode + length(phy@phylo$tip.label))
nodes.text[nodes.to.label] <- c(1:length(nodes.to.label))

#for making a legend for calibration points
legend.df <- data.frame(
  nodelab = c(1:11),
  text = c("Patrinia hokiana fruit (5.34-7.24)", "Dipelta bovayensis fruit (48-49)", "Lonicera krassilovii leaf (4.5-15.3)", "Diervilla sp. pollen (40-42)", "Lonicera protojaponica leaf (5.34-7.24)", "Lonicera sp. pollen (34.07 ? 0.10)", "Heptacodium spinosa pollen (20.4-23)", "Weigela weichangensis seed (22.1)", "Viburnum pollen type 1b (45)", "Viburnum pollen type 1a (47)", "Secondary calibration (49-130)"),
  x = rep(-8.5, 11),#to get suitable position for points in the legend on x axis; repeat 126-180 11 times
  y = seq(58, (58 - (10 * 1.5)), -1.5))#to get suitable position for points in the legend on y axis-the y axis will go from 58 to 28, every 2 units#

#using groupOTU function to colour selected tip labels in different colours
tree <- groupOTU(phy, taxa_to_color)
p <- ggtree(tree, size  = 0)#size refers to edge width which is set to 0 here because q ill be overwritten by p later and we don't want to have two replaictes of edge lengths

# ggtree incorrectly rotates 95%HPD
#see this description which has been taken from Michael Landis's github page from 
#vib_div repository- https://github.com/mlandis/vib_div/blob/master/code/plot/plot_fig3_mcc.R
## Encountered problems with using geom_range to plot age HPDs in ggtree. It
# appears that geom_range incorrectly rotates the HPD relative to the height
# of the node unnecessarily. My guess for this would be because older version
# of ggtree primarily supported length measurements, and not height measurements
# so the new capability to handle height might contain a "reflection" bug.
#
# For example, suppose a node has height 3 with HPD [2, 7]. You can think of
# this offset as h + [2-h, 7-h]. ggtree seems to "rotate" this
# causing the HPD to appear as [-1, 4]. Figtree displays this correctly.
# 
# See this excellent trick by Tauana: https://groups.google.com/forum/#!msg/bioc-ggtree/wuAlY9phL9Q/L7efezPgDAAJ

#To solve this problem, we will manually extract and plot HPD data from p ggplot object
hpd <-  p +
  geom_tree() +
  geom_rootedge(rootedge = 9) +
  coord_cartesian(xlim = c(-13, 137+52),
                  ylim = c(-3.5, 59+1.5),
                  expand = FALSE) +
  geom_range("height_0.95_HPD", color='blue', size=0.5, alpha=.3)

hpd <- hpd$data[, c("height", "height_0.95_HPD", "x", "y")]
hpd$min <- NA
hpd$max <- NA
for(i in 1:nrow(hpd)){
  hpd$min[i] <- hpd$height_0.95_HPD[i][[1]][1]
  hpd$max[i] <- hpd$height_0.95_HPD[i][[1]][2]
}

hpd$min <- max(phy@phylo$edge.length) - hpd$min
hpd$max <- max(phy@phylo$edge.length) - hpd$max

#here we are plotting names and alternating bands of grey for different epochs
#make sure that all the grey scales are plotted first and then the texts otherwise the text will hide behind the grey rectangles#
q <- p + annotate(geom = "rect",
                  xmin = epochs2$max_age[3], xmax = epochs2$min_age[1],
                  ymin = -4, ymax = 59+1.5,
                  fill = "gray95") +
  annotate(geom = "rect",
           xmin = epochs2$max_age[5], xmax = epochs2$min_age[5],
           ymin = -7, ymax = 59+1.5,
           fill = "gray95") +
  annotate(geom = "rect",
           xmin = epochs2$max_age[7], xmax = epochs2$min_age[7],
           ymin = -7, ymax = 59+1.5,
           fill = "gray95") +
  annotate(geom = "rect",
           xmin = epochs2$max_age[9], xmax = epochs2$min_age[9],
           ymin = -7, ymax = 59+1.5,
           fill = "gray95") +
  annotate(geom = "text",
           x = mean(c(epochs2$max_age[3], epochs2$min_age[1])),
           y = -2,
           label = "Rec",
           hjust = 0.5, size = 2.5) +
  annotate(geom = "text",
           x = mean(c(epochs2$max_age[4], epochs2$min_age[4])),
           y = -2,
           label = "Mio",
           hjust = 0.5, size = 2.5) +
  annotate(geom = "text",
           x = mean(c(epochs2$max_age[5], epochs2$min_age[5])),
           y = -2,
           label = "Oli",
           hjust = 0.5, size = 2.5) +
  annotate(geom = "text",
           x = mean(c(epochs2$max_age[6], epochs2$min_age[6])),
           y = -2,
           label = "Eoc",
           hjust = 0.5, size = 2.5) +
  annotate(geom = "text",
           x = mean(c(epochs2$max_age[7], epochs2$min_age[7])),
           y = -2,
           label = "Pal",
           hjust = 0.5, size = 2.5) +
  annotate(geom = "text",
           x = mean(c(epochs2$max_age[8], epochs2$min_age[8])),
           y = -2,
           label = "Late Cret",
           hjust = 0.5, size = 2.5) +
  annotate(geom = "text",
           x = mean(c(epochs2$max_age[9], epochs2$min_age[9])),
           y = -2,
           label = "Early Cret",
           hjust = 0.5, size = 2.5) 

# to get minor ticks on x-axis
breaks <- (max(phy@phylo$edge.length) - c(seq(140, 0, -5)))
labels <- seq(140, 0, -5)
labs <- NA
for(i in 1:length(labels)) {
  if (labels[i] %in% seq(140, 0, -20)) {
    labs[i] <- labels[i]
  }
}
labs[is.na(labs)] <- ""

q <- q +
  geom_tree(size = 0.35) +
  geom_text2(aes(label = label, 
                 subset = isTip, 
                 color = group), 
             size = 2.5, offset=0, hjust = 0, angle = 0, offset.text = 1.1,
             family = "Helvetica",
             fontface = 4) +
  scale_color_manual(values = c("#e67e00", "#009E73")) +
  geom_rootedge(rootedge = 9.8) +
  coord_cartesian(xlim = c(-13, 137+52),
                  ylim = c(-3.5, 59+1.5),
                  expand = FALSE) +
  scale_x_continuous(
    breaks = breaks,
    labels = labs) +
  labs(x = "Age (Ma)") +
  geom_segment(data = hpd, 
               aes(x = min, y = y, xend = max, yend = y),
               size = 1.5, alpha = 0.3, color = "blue") +
  geom_point(data = ggtree(phy@phylo)$data, 
             aes(x, y), 
             color = node.label.cols,
             size = 5, shape = 21, fill = "#ffffc6", stroke = 0.5) +
  geom_text(data = ggtree(phy@phylo)$data,
            aes(x, y),
            label = nodes.text, size = 3,
            hjust = 0.5, vjust = 0.5) +
#  theme_tree2() +
  theme(
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, "in"),
    axis.line.x = element_line(color = "black", size = 0.25),
    axis.ticks.x = element_line(color = "black", size = 0.25),
    axis.text.x = element_text(color = "black", size = 7),
    axis.title.x = element_text(color = "black", size = 8),
    legend.position = "none"
  )

  
#to plot legend for calibration points
for(i in 1:nrow(legend.df)) {
  
  q <- q +
    annotate(geom = "point", 
             x = legend.df$x[i], y = legend.df$y[i],
             shape = 21, fill = "#ffffc6", size = 5, color = "Black", stroke = 0.3) +
    annotate(geom = "text", 
             x = legend.df$x[i], y = legend.df$y[i], 
             label = legend.df$nodelab[i], 
             size = 3,
             hjust = 0.5, vjust = 0.5) +
    annotate(geom = "text",
             x = legend.df$x[i] + 4, y = legend.df$y[i],
             label = legend.df$text[i],
             size = 2.5,family = "Helvetica", fontface = 3,
             hjust = 0, vjust = 0.5)
}

q <- q + annotate(geom = "text",
                  x = max(phy@phylo$edge.length)-130, y = 60,
                  label = "Calibration (Ma)",
                  size = 3, family = "Helvetica", fontface = 1,
                  hjust = 0, vjust = 0.5) # this needs to be outside of the loop, otherwise it gets plotted multiple times.

#to save the plot#
ggsave(q, filename = "DatedNuclearTree.pdf", width = 7.25 , height = 9, units = "in")




