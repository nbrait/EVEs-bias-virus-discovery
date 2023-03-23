library(ggplot2)
setwd("C:/Users/nadja/Documents/R/trees_for_quaranja_new/THRESHOLD")
matches_table <- read.csv("project_thresholds_new", header = TRUE, sep = "\t")
ggplot(matches_table[metadata_threshold$PROJECT=="PRJNA564787",]) + 
  geom_jitter(aes(DIFF_BITSCORE, pmax(EXO_BITSCORE, -EVE_BITSCORE), color= EVE_SIMILARITY/(EVE_SIMILARITY-EXO_SIMILARITY)), size=4) +
  scale_color_viridis_c(name = "identity to known EVE", 
                        limits = c(0,1)) +
  scale_y_log10(limits = c(10,1800)) +
  theme(legend.position = "none",
        legend.direction = "vertical",panel.background = element_rect(fill = "white"), 
        plot.title = element_text(hjust = 0.5),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1)) +
  xlab("Δ BLAST bit score") +# for the x axis label
  ylab("maximum BLAST bit score per contig (log10)") + xlim(-400, 1200)+
  annotate("text", x = -400, y = 1800, label = "PRJNA564787", size = 15, hjust = 0)

ggplot(metadata_threshold[metadata_threshold$PROJECT=="PRJNA605178",]) + 
  geom_jitter(aes(DIFF_BITSCORE, pmax(EXO_BITSCORE, -EVE_BITSCORE), color= EVE_SIMILARITY/(EVE_SIMILARITY-EXO_SIMILARITY)), size=4) +
  scale_color_viridis_c(name = "identity to known EVE", 
                        limits = c(0,1)) +
  scale_y_log10(limits = c(10,1800)) +
  theme(legend.position = "none",
        legend.direction = "vertical",panel.background = element_rect(fill = "white"), 
        plot.title = element_text(hjust = 0.5),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1)) +
  xlab("Δ BLAST bit score") +# for the x axis label
  ylab("maximum BLAST bit score per contig (log10)") + xlim(-400, 1200)+
  annotate("text", x = -400, y = 1800, label = "PRJNA605178", size = 15, hjust = 0)

ggplot(metadata_threshold[metadata_threshold$PROJECT=="PRJNA233429",]) + 
  geom_jitter(aes(DIFF_BITSCORE, pmax(EXO_BITSCORE, -EVE_BITSCORE), color= EVE_SIMILARITY/(EVE_SIMILARITY-EXO_SIMILARITY)), size=4) +
  scale_color_viridis_c(name = "identity to known EVE", 
                        limits = c(0,1)) +
  scale_y_log10(limits = c(10,1800)) +
  theme(legend.position = "none",
        legend.direction = "vertical",panel.background = element_rect(fill = "white"), 
        plot.title = element_text(hjust = 0.5),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1)) +
  xlab("Δ BLAST bit score") +# for the x axis label
  ylab("maximum BLAST bit score per contig (log10)") + xlim(-400, 1200)+
  annotate("text", x = -400, y = 1800, label = "PRJNA233429", size = 15, hjust = 0)

ggplot(metadata_threshold[metadata_threshold$PROJECT=="PRJNA276135",]) + 
  geom_jitter(aes(DIFF_BITSCORE, pmax(EXO_BITSCORE, -EVE_BITSCORE), color= EVE_SIMILARITY/(EVE_SIMILARITY-EXO_SIMILARITY)), size=4) +
  scale_color_viridis_c(name = "identity to known EVE", 
                        limits = c(0,1)) +
  scale_y_log10(limits = c(10,1800)) +
  theme(legend.position = "none",
        legend.direction = "vertical",panel.background = element_rect(fill = "white"), 
        plot.title = element_text(hjust = 0.5),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1)) +
  xlab("Δ BLAST bit score") +# for the x axis label
  ylab("maximum BLAST bit score per contig (log10)") + xlim(-400, 1200)+
  annotate("text", x = -400, y = 1800, label = "PRJNA276135", size = 15, hjust = 0)

ggplot(metadata_threshold[metadata_threshold$PROJECT=="PRJNA310013",]) + 
  geom_jitter(aes(DIFF_BITSCORE, pmax(EXO_BITSCORE, -EVE_BITSCORE), color= EVE_SIMILARITY/(EVE_SIMILARITY-EXO_SIMILARITY)), size=4) +
  scale_color_viridis_c(name = "identity to known EVE", 
                        limits = c(0,1)) +
  scale_y_log10(limits = c(10,1800)) +
  theme(legend.position = "none",
        legend.direction = "vertical",panel.background = element_rect(fill = "white"), 
        plot.title = element_text(hjust = 0.5),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1)) +
  xlab("Δ BLAST bit score") +# for the x axis label
  ylab("maximum BLAST bit score per contig (log10)") + xlim(-400, 1200)+
  annotate("text", x = -400, y = 1800, label = "PRJNA310013", size = 15, hjust = 0)

ggplot(metadata_threshold[metadata_threshold$PROJECT=="PRJNA356542",]) + 
  geom_jitter(aes(DIFF_BITSCORE, pmax(EXO_BITSCORE, -EVE_BITSCORE), color= EVE_SIMILARITY/(EVE_SIMILARITY-EXO_SIMILARITY)), size=4) +
  scale_color_viridis_c(name = "identity to known EVE", 
                        limits = c(0,1)) +
  scale_y_log10(limits = c(10,1800)) +
  theme(legend.position = "none",
        legend.direction = "vertical",panel.background = element_rect(fill = "white"), 
        plot.title = element_text(hjust = 0.5),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1)) +
  xlab("Δ BLAST bit score") +# for the x axis label
  ylab("maximum BLAST bit score per contig (log10)") + xlim(-400, 1200)+
  annotate("text", x = -400, y = 1800, label = "PRJNA356542", size = 15, hjust = 0)

ggplot(metadata_threshold[metadata_threshold$PROJECT=="PRJNA388696",]) + 
  geom_jitter(aes(DIFF_BITSCORE, pmax(EXO_BITSCORE, -EVE_BITSCORE), color= EVE_SIMILARITY/(EVE_SIMILARITY-EXO_SIMILARITY)), size=4) +
  scale_color_viridis_c(name = "identity to known EVE", 
                        limits = c(0,1)) +
  scale_y_log10(limits = c(10,1800)) +
  theme(legend.position = "none",
        legend.direction = "vertical",panel.background = element_rect(fill = "white"), 
        plot.title = element_text(hjust = 0.5),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1)) +
  xlab("Δ BLAST bit score") +# for the x axis label
  ylab("maximum BLAST bit score per contig (log10)") + xlim(-400, 1200)+
  annotate("text", x = -400, y = 1800, label = "PRJNA388696", size = 15, hjust = 0)

ggplot(metadata_threshold[metadata_threshold$PROJECT=="PRJNA413709",]) + 
  geom_jitter(aes(DIFF_BITSCORE, pmax(EXO_BITSCORE, -EVE_BITSCORE), color= EVE_SIMILARITY/(EVE_SIMILARITY-EXO_SIMILARITY)), size=4) +
  scale_color_viridis_c(name = "identity to known EVE", 
                        limits = c(0,1)) +
  scale_y_log10(limits = c(10,1800)) +
  theme(legend.position = "none",
        legend.direction = "vertical",panel.background = element_rect(fill = "white"), 
        plot.title = element_text(hjust = 0.5),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1)) +
  xlab("Δ BLAST bit score") +# for the x axis label
  ylab("maximum BLAST bit score per contig (log10)") + xlim(-400, 1200)+
  annotate("text", x = -400, y = 1800, label = "PRJNA413709", size = 15, hjust = 0)

ggplot(metadata_threshold[metadata_threshold$PROJECT=="PRJNA515586",]) + 
  geom_jitter(aes(DIFF_BITSCORE, pmax(EXO_BITSCORE, -EVE_BITSCORE), color= EVE_SIMILARITY/(EVE_SIMILARITY-EXO_SIMILARITY)), size=4) +
  scale_color_viridis_c(name = "identity to known EVE", 
                        limits = c(0,1)) +
  scale_y_log10(limits = c(10,1800)) +
  theme(legend.position = "none",
        legend.direction = "vertical",panel.background = element_rect(fill = "white"), 
        plot.title = element_text(hjust = 0.5),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1)) +
  xlab("Δ BLAST bit score") +# for the x axis label
  ylab("maximum BLAST bit score per contig (log10)") + xlim(-400, 1200)+
  annotate("text", x = -400, y = 1800, label = "PRJNA515586", size = 15, hjust = 0)

ggplot(metadata_threshold[metadata_threshold$PROJECT=="PRJNA547758",]) + 
  geom_jitter(aes(DIFF_BITSCORE, pmax(EXO_BITSCORE, -EVE_BITSCORE), color= EVE_SIMILARITY/(EVE_SIMILARITY-EXO_SIMILARITY)), size=4) +
  scale_color_viridis_c(name = "identity to known EVE", 
                        limits = c(0,1)) +
  scale_y_log10(limits = c(10,1800)) +
  theme(legend.position = "none",
        legend.direction = "vertical",panel.background = element_rect(fill = "white"), 
        plot.title = element_text(hjust = 0.5),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1)) +
  xlab("Δ BLAST bit score") +# for the x axis label
  ylab("maximum BLAST bit score per contig (log10)") + xlim(-400, 1200)+
  annotate("text", x = -400, y = 1800, label = "PRJNA547758", size = 15, hjust = 0)

ggplot(metadata_threshold[metadata_threshold$PROJECT=="PRJNA814001",]) + 
  geom_jitter(aes(DIFF_BITSCORE, pmax(EXO_BITSCORE, -EVE_BITSCORE), color= EVE_SIMILARITY/(EVE_SIMILARITY-EXO_SIMILARITY)), size=4) +
  scale_color_viridis_c(name = "identity to known EVE", 
                        limits = c(0,1)) +
  scale_y_log10(limits = c(10,1800)) +
  theme(legend.position = "none",
        legend.direction = "vertical",panel.background = element_rect(fill = "white"), 
        plot.title = element_text(hjust = 0.5),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1)) +
  xlab("Δ BLAST bit score") +# for the x axis label
  ylab("maximum BLAST bit score per contig (log10)") + xlim(-400, 1200)+
  annotate("text", x = -400, y = 1800, label = "PRJNA814001", size = 15, hjust = 0)

ggplot(metadata_threshold[metadata_threshold$PROJECT=="PRJNA413709",]) + 
  geom_jitter(aes(DIFF_BITSCORE, pmax(EXO_BITSCORE, -EVE_BITSCORE), color= EVE_SIMILARITY/(EVE_SIMILARITY-EXO_SIMILARITY)), size=4) +
  scale_color_viridis_c(name = "identity to known EVE", 
                        limits = c(0,1)) +
  scale_y_log10(limits = c(10,1800)) +
  theme(legend.position = "none",
        legend.direction = "vertical",panel.background = element_rect(fill = "white"), 
        plot.title = element_text(hjust = 0.5),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1)) +
  xlab("Δ BLAST bit score") +# for the x axis label
  ylab("maximum BLAST bit score per contig (log10)") + xlim(-400, 1200)+
  annotate("text", x = -400, y = 1800, label = "PRJNA413709", size = 15, hjust = 0)

ggplot(metadata_threshold[metadata_threshold$PROJECT=="PRJNA431084",]) + 
  geom_jitter(aes(DIFF_BITSCORE, pmax(EXO_BITSCORE, -EVE_BITSCORE), color= EVE_SIMILARITY/(EVE_SIMILARITY-EXO_SIMILARITY)), size=4) +
  scale_color_viridis_c(name = "identity to known EVE", 
                        limits = c(0,1)) +
  scale_y_log10(limits = c(10,1800)) +
  theme(legend.position = "none",
        legend.direction = "vertical",panel.background = element_rect(fill = "white"), 
        plot.title = element_text(hjust = 0.5),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1)) +
  xlab("Δ BLAST bit score") +# for the x axis label
  ylab("maximum BLAST bit score per contig (log10)") + xlim(-400, 1200)+
  annotate("text", x = -400, y = 1800, label = "PRJNA431084", size = 15, hjust = 0)
  
ggplot(metadata_threshold[metadata_threshold$PROJECT=="PRJNA208524",]) + 
  geom_jitter(aes(DIFF_BITSCORE, pmax(EXO_BITSCORE, -EVE_BITSCORE), color= EVE_SIMILARITY/(EVE_SIMILARITY-EXO_SIMILARITY)), size=4) +
  scale_color_viridis_c(name = "identity to known EVE", 
                        limits = c(0,1)) +
  scale_y_log10(limits = c(10,1800)) +
  theme(legend.position = "none",
        legend.direction = "vertical",panel.background = element_rect(fill = "white"), 
        plot.title = element_text(hjust = 0.5),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1)) +
  xlab("Δ BLAST bit score") +# for the x axis label
  ylab("maximum BLAST bit score per contig (log10)") + xlim(-400, 1200)+
  annotate("text", x = -400, y = 1800, label = "PRJNA208524", size = 15, hjust = 0)

ggplot(metadata_threshold[metadata_threshold$PROJECT=="PRJNA483089",]) + 
  geom_jitter(aes(DIFF_BITSCORE, pmax(EXO_BITSCORE, -EVE_BITSCORE), color= EVE_SIMILARITY/(EVE_SIMILARITY-EXO_SIMILARITY)), size=4) +
  scale_color_viridis_c(name = "identity to known EVE", 
                        limits = c(0,1)) +
  scale_y_log10(limits = c(10,1800)) +
  theme(legend.position = "none",
        legend.direction = "vertical",panel.background = element_rect(fill = "white"), 
        plot.title = element_text(hjust = 0.5),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1)) +
  xlab("Δ BLAST bit score") +# for the x axis label
  ylab("maximum BLAST bit score per contig (log10)") + xlim(-400, 1200)+
  annotate("text", x = -400, y = 1800, label = "PRJNA483089", size = 15, hjust = 0)

