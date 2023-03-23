setwd("C:/Users/nadja/Documents/R/trees_for_quaranja_new/MAPS")
#install.packages("sf")
#install.packages("mapdata")
library(dplyr)
library(ggplot2)
library(sf)
library(maps)
library(ggrepel)
library(patchwork)
############################################
# reading in table
City <- read.csv("screened_samples_coord_new.txt", header = TRUE, sep = "\t") 

# reading in world data
world_map <- map_data("world")
world_map <- mutate(world_map, fill = ifelse(region %in% c("China", "Brazil"), "dark gray", "light grey"))
# plot map data
p <- ggplot(label=cities$geo_loc_name_country) + coord_fixed() +
  xlab("") + ylab("") + scale_y_continuous(limits=c(-60,90)) + scale_x_continuous(limits=c(-150,170))

# making base world
base_world_messy <- p + geom_polygon(data=world_map, aes(x=long, y=lat, group=group, fill=fill), 
                                     colour="white") + scale_fill_identity()
cleanup <- 
  theme(text = element_text(size = 20), 
        legend.key = element_blank(), 
        legend.position = c(0.1, 0.25), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white', colour = 'white'), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        legend.box = "horizontal",
        legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size=14), #change legend text font size
        legend.background=element_blank())

base_world <- base_world_messy + cleanup
## transcriptomic data
map_data_sized <- 
  base_world +
  geom_point(na.rm = TRUE, data= City[City$Source=="RNA-Seq",], aes(x=long, y=lat, size= Freq, color= Tribes),alpha= 0.4, position=position_jitter(h=1,w=2)) + 
  scale_size(range = c(4, 10), breaks = c(1,5,10,20,30)) +
  scale_color_manual(values= c("red3", "dodgerblue3", "gray5")) +
  labs(size="sample size \n per location") +
  annotate("text", x=-50, y=-43, label="screened SRA datasets: 538", color="black", size= 6, hjust =0) +
  annotate("text", x=-50, y=-50, label="unknown location: 10", color="black", size= 6, hjust = 0) +
  annotate("text", x=100, y=36, label="n=117", color="black", size= 5) +
  annotate("text", x=-50, y=-12, label="n=7", color="black", size= 5) +
  annotate("text", x=-150, y=80, label="a)", color="black", size= 6)

map_data_sized

## Orthomyxovirus hits

# reading in table
cities <- read.csv("screened_hit_samples.txt", header = TRUE, sep = "\t") 

# reading in world data
world_map <- map_data("world")
# plot map data
p <- ggplot(label=cities$geo_loc_name_country) + coord_fixed() +
  xlab("") + ylab("") + scale_y_continuous(limits=c(-60,90)) + scale_x_continuous(limits=c(-150,170))

# making base world
base_world_messy <- p + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
                                     colour="white", fill="light grey")
cleanup <- 
  theme(panel.grid.major = element_blank(), text = element_text(size = 20), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white', colour = 'white'), 
        legend.position = c(0.1, 0.25), legend.box = "horizontal", 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size=14),
        legend.background=element_blank())

base_world <- base_world_messy + cleanup

## putting in points data

map_data_sized2 <- 
  base_world +
  geom_point(na.rm = TRUE, data= cities[cities$Source=="SRR",], aes(x=lng, y=lat,size = Freq, color= Tribes),alpha= 0.4, position=position_jitter(h=1,w=2)) +
  scale_color_manual(values= c("red3", "dodgerblue3", "gray5")) +
  scale_size(range = c(4, 7), breaks = c(1,5,10)) +
  labs(size="sample size \n per location") +
  scale_shape_manual(values=c(19,8)) +
  scale_alpha_manual(values = c(0.3, 1)) + 
  geom_text_repel(na.rm = TRUE, data = cities, size =5, aes(x=lng, y=lat, label= geo_loc_name),hjust=0, vjust=0, max.overlaps = 20) +
  annotate("text", x=-50, y=-43, label="positive contig hits: 95", color="black", size= 6, hjust =0) +
  annotate("text", x=-50, y=-50, label="unknown location: 3", color="black", size= 6, hjust = 0) + 
  annotate("text", x=-150, y=80, label="b)", color="black", size= 6)

map_data_sized2

map_data_sized / map_data_sized2
