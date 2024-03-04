library(readxl)
library(openxlsx)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

data_transformed <- read.csv("../Figures/HeatMap/HeatMap_transformed.csv")
### DEFINE DATA ORDER BY CONSEQUENCE 
data_transformed$order <- ifelse(data_transformed$CONSEQUENCE == "SPLICE", 1,
                                 ifelse(data_transformed$CONSEQUENCE == "MISSENSE", 2,
                                        ifelse(data_transformed$CONSEQUENCE == "MISSENSE, SPLICE REGION" , 3,
                                        ifelse(data_transformed$CONSEQUENCE == "NONSENSE",4,
                                        ifelse(data_transformed$CONSEQUENCE == "NONSENSE, SPLICE REGION", 5, 6)))))

data_transformed <- data_transformed[order(data_transformed$order, decreasing = TRUE),]
a <- unique(data_transformed[,c("CONSEQUENCE", "VARIANT")])
variant_order <- a$VARIANT
consequence_order <- a$CONSEQUENCE

values = c("#D53E4F","#FC8D59", "#B2DF8A","#F8766D00")
data_subset <- read.csv("../Figures/HeatMap/HeatMap_transformed2.csv")
data_subset <- data_subset[data_subset$EFFECT == "CONSEQUENCE",c("VARIANT", "EFFECT", "RANGE")]
data_subset$order <- ifelse(data_subset$RANGE == "SPLICE", 1,
                            ifelse(data_subset$RANGE == "MISSENSE", 2,
                                   ifelse(data_subset$RANGE == "MISSENSE, SPLICE REGION" , 3,
                                          ifelse(data_subset$RANGE == "NONSENSE",4,
                                                 ifelse(data_subset$RANGE == "NONSENSE, SPLICE REGION", 5, 6)))))



g_main <- ggplot(data_transformed, aes(x=VARIANT, y=EFFECT, fill=RANGE)) +
  geom_tile() +
  coord_flip() +
  scale_x_discrete(limits = variant_order) +
  scale_y_discrete(limits = c("SpliceAI_AG" ,"SpliceAI_AL","SpliceAI_DG","SpliceAI_DL",
                              "MaxEntScan_AG", "MaxEntScan_AL","MaxEntScan_DG", "MaxEntScan_DL"),
                   position = "left") +
  scale_fill_manual(breaks = c("No impact", "Low", "Intermediate", "High"),
                    values = c("#F8766D00","#B2DF8A","#FC8D59","#D53E4F"),
                    name = "Variant Impact",
                    labels = c("No impact", "Low", "Intermediate", "High")) + 
  theme(axis.text.x = element_text(angle = 90, hjust=1, size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        aspect.ratio = 14/4,
        axis.text.y = element_text(size = rel(1)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(30,15,30,15), "pt"),
        legend.position = "bottom",
        legend.text=element_text(size=rel(1))) 
g_main

g_left <- ggplot() +
  theme_void() +
  coord_flip() +
  geom_tile(data = data_subset,
            aes(x = VARIANT, y = 1, fill = RANGE, width = 1)) +
  scale_fill_manual(breaks = c("SPLICE", "MISSENSE", "MISSENSE, SPLICE REGION", "NONSENSE","NONSENSE, SPLICE REGION"),
    values = c("#FB9A99","#A6CEE3","#1F78B4","#CAB2D6","#6A3D9A"),
    name = "Consequence",
    labels = c("SPLICE", "MISSENSE", "MISSENSE, SPLICE REGION", "NONSENSE","NONSENSE, SPLICE REGION")) +
  scale_x_discrete(limits = variant_order) +
  scale_y_discrete(limits = c("CONSEQUENCE")) +
  theme(aspect.ratio = 16/0.685,
        legend.position = "right",
        legend.text=element_text(size=rel(1)),
        plot.margin=unit(c(30,15,30,15), "pt"))

g_left

g_main + g_left

ggsave("FINAL_FIGURES_AMR_legendPaper.png",
       plot = last_plot(),
       device = "png",
       width = 12,
       height = 30)


## GROUP SPLICE:
data_transformed_sub <- data_transformed[data_transformed$CONSEQUENCE == "SPLICE",]
a_sub <- unique(data_transformed_sub[,c("CONSEQUENCE", "VARIANT")])
variant_order <- a_sub$VARIANT
consequence_order <- a_sub$CONSEQUENCE

data_subset <- data_subset[data_subset$EFFECT == "CONSEQUENCE",c("VARIANT", "EFFECT", "RANGE")]
data_subset_sub <- data_subset[data_subset$VARIANT %in% variant_order,]

g_SPLICE<- ggplot(data_transformed_sub, aes(x=VARIANT, y=EFFECT, fill=RANGE)) +
  geom_tile() +
  coord_flip() +
  scale_x_discrete(limits = variant_order) +
  scale_y_discrete(limits = c("SpliceAI_AG" ,"SpliceAI_AL","SpliceAI_DG","SpliceAI_DL",
                              "MaxEntScan_AG", "MaxEntScan_AL","MaxEntScan_DG", "MaxEntScan_DL"),
                   position = "left") +
  scale_fill_manual(breaks = c("No impact", "Low", "Intermediate", "High"),
                    values = c("#F8766D00","#B2DF8A","#FC8D59","#D53E4F"),
                    name = "Variant Impact",
                    labels = c("No impact", "Low", "Intermediate", "High")) + 
  theme(axis.text.x = element_text(angle = 90, hjust=1, size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        aspect.ratio = 14/4,
        axis.text.y = element_text(size = rel(1)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(30,30,30,30), "pt"),
        legend.position  = c(0.4, -0.145),
        legend.direction = 'horizontal',
        legend.text=element_text(size=rel(1))) 
g_SPLICE
ggsave("FINAL_FIGURES_Splice.png",
       plot = last_plot(),
       device = "png",
       width = 7,
       height = 12)

g_left <- ggplot() +
  theme_void() +
  coord_flip() +
  geom_tile(data = data_subset_sub,
            aes(x = VARIANT, y = 1, fill = RANGE, width = 1)) +
  scale_fill_manual(breaks = c("SPLICE", "MISSENSE", "MISSENSE, SPLICE REGION", "NONSENSE","NONSENSE, SPLICE REGION"),
                    values = c("#FB9A99","#A6CEE3","#1F78B4","#CAB2D6","#6A3D9A"),
                    name = "Consequence",
                    labels = c("SPLICE", "MISSENSE", "MISSENSE, SPLICE REGION", "NONSENSE","NONSENSE, SPLICE REGION")) +
  scale_x_discrete(limits = variant_order) +
  scale_y_discrete(limits = c("CONSEQUENCE")) +
  theme(aspect.ratio = 16/0.685,
        legend.position = "right",
        legend.text=element_text(size=rel(1)),
        plot.margin=unit(c(30,15,30,15), "pt"))

g_left

g_SPLICE + g_left
ggsave("/media/adminiis/HematoLaFe/variantes_Papaemmanuil/SplicingVariants/FIGURES/HeatMap/FINAL_FIGURES_Splice.png",
       plot = last_plot(),
       device = "png",
       width = 7,
       height = 12)

### GROUP MISSENSE
data_transformed_sub <- data_transformed[data_transformed$CONSEQUENCE == "MISSENSE" | data_transformed$CONSEQUENCE == "MISSENSE, SPLICE REGION",]
a_sub <- unique(data_transformed_sub[,c("CONSEQUENCE", "VARIANT")])
variant_order <- a_sub$VARIANT
consequence_order <- a_sub$CONSEQUENCE


data_subset <- data_subset[data_subset$EFFECT == "CONSEQUENCE",c("VARIANT", "EFFECT", "RANGE")]
data_subset_sub <- data_subset[data_subset$VARIANT %in% variant_order,]

g_MISSENSE <- ggplot(data_transformed_sub, aes(x=VARIANT, y=EFFECT, fill=RANGE)) +
  geom_tile() +
  coord_flip() +
  scale_x_discrete(limits = variant_order) +
  scale_y_discrete(limits = c("SpliceAI_AG" ,"SpliceAI_AL","SpliceAI_DG","SpliceAI_DL",
                              "MaxEntScan_AG", "MaxEntScan_AL","MaxEntScan_DG", "MaxEntScan_DL"),
                   position = "left") +
  scale_fill_manual(breaks = c("No impact", "Low", "Intermediate", "High"),
                    values = c("#F8766D00","#B2DF8A","#FC8D59","#D53E4F"),
                    name = "Variant Impact",
                    labels = c("No impact", "Low", "Intermediate", "High")) + 
  theme(axis.text.x = element_text(angle = 90, hjust=1, size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        aspect.ratio = 14/4,
        axis.text.y = element_text(size = rel(1)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(30,30,30,30), "pt"),
        legend.position  = "bottom",
        legend.text=element_text(size=rel(1))) 
g_MISSENSE

g_left <- ggplot() +
  theme_void() +
  coord_flip() +
  geom_tile(data = data_subset_sub,
            aes(x = VARIANT, y = 1, fill = RANGE, width = 1)) +
  scale_fill_manual(breaks = c("SPLICE", "MISSENSE", "MISSENSE, SPLICE REGION", "NONSENSE","NONSENSE, SPLICE REGION"),
                    values = c("#FB9A99","#A6CEE3","#1F78B4","#CAB2D6","#6A3D9A"),
                    name = "Consequence",
                    labels = c("SPLICE", "MISSENSE", "MISSENSE, SPLICE REGION", "NONSENSE","NONSENSE, SPLICE REGION")) +
  scale_x_discrete(limits = variant_order) +
  scale_y_discrete(limits = c("CONSEQUENCE")) +
  theme(aspect.ratio = 16/0.685,
        legend.position = "right",
        legend.text=element_text(size=rel(1)),
        plot.margin=unit(c(30,15,30,15), "pt"))

g_left

g_MISSENSE + g_left

ggsave("FINAL_FIGURES_Missense&SpliceRegion.png",
       plot = last_plot(),
       device = "png",
       width = 9,
       height = 10)

### GROUP NONSENSE
data_transformed_sub <- data_transformed[data_transformed$CONSEQUENCE == "NONSENSE" | data_transformed$CONSEQUENCE == "NONSENSE, SPLICE REGION",]
a_sub <- unique(data_transformed_sub[,c("CONSEQUENCE", "VARIANT")])
variant_order <- a_sub$VARIANT
consequence_order <- a_sub$CONSEQUENCE


data_subset <- data_subset[data_subset$EFFECT == "CONSEQUENCE",c("VARIANT", "EFFECT", "RANGE")]
data_subset_sub <- data_subset[data_subset$VARIANT %in% variant_order,]

g_NONSENSE <- ggplot(data_transformed_sub, aes(x=VARIANT, y=EFFECT, fill=RANGE)) +
  geom_tile() +
  coord_flip() +
  scale_x_discrete(limits = variant_order) +
  scale_y_discrete(limits = c("SpliceAI_AG" ,"SpliceAI_AL","SpliceAI_DG","SpliceAI_DL",
                              "MaxEntScan_AG", "MaxEntScan_AL","MaxEntScan_DG", "MaxEntScan_DL"),
                   position = "left") +
  scale_fill_manual(breaks = c("No impact", "Low", "Intermediate", "High"),
                    values = c("#F8766D00","#B2DF8A","#FC8D59","#D53E4F"),
                    name = "Variant Impact",
                    labels = c("No impact", "Low", "Intermediate", "High")) + 
  theme(axis.text.x = element_text(angle = 90, hjust=1, size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        aspect.ratio = 14/4,
        axis.text.y = element_text(size = rel(1)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(30,30,30,30), "pt"),
        legend.position  = "bottom",
        legend.text=element_text(size=rel(1))) 
g_MISSENSE

g_left <- ggplot() +
  theme_void() +
  coord_flip() +
  geom_tile(data = data_subset_sub,
            aes(x = VARIANT, y = 1, fill = RANGE, width = 1)) +
  scale_fill_manual(breaks = c("SPLICE", "MISSENSE", "MISSENSE, SPLICE REGION", "NONSENSE","NONSENSE, SPLICE REGION"),
                    values = c("#FB9A99","#A6CEE3","#1F78B4","#CAB2D6","#6A3D9A"),
                    name = "Consequence",
                    labels = c("SPLICE", "MISSENSE", "MISSENSE, SPLICE REGION", "NONSENSE","NONSENSE, SPLICE REGION")) +
  scale_x_discrete(limits = variant_order) +
  scale_y_discrete(limits = c("CONSEQUENCE")) +
  theme(aspect.ratio = 16/0.685,
        legend.position = "right",
        legend.text=element_text(size=rel(1)),
        plot.margin=unit(c(30,15,30,15), "pt"))

g_left

g_NONSENSE + g_left

ggsave("FINAL_FIGURES_Nonsense&SpliceRegion.png",
       plot = last_plot(),
       device = "png",
       width = 9,
       height = 10)
