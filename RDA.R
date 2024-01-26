#load packages
library(vegan)
library(EnvStats)
library(ggplot2)
library(ggfortify)

#set working directory
setwd("~/Desktop/TXST/Haley/")

#import data (row names are sample IDs and columns are immune assay metrics)
data <- read.csv("assay_rawData.csv", row.names = 1)

#test for outliers (rosner test) and normality (shapiro-wilk p>0.05 data are normal)
#--------SYMBIONT DENSITY-----------
hist(data$symbiontDensity, breaks = 10, main = "symbiont density")
rosnerTest(data$symbiontDensity, k=3)
#no outliers in symbiont density
shapiro.test(data$symbiontDensity)
#data are normal with no transformation

#--------PPO-----------
hist(data$PPO_mgProtein, breaks = 10, main = "PPO")
rosnerTest(data$PPO_mgProtein, k=15)
#12 outliers in PPO filtered out
PPO_data <- subset(data, PPO_mgProtein > 0)
PPO_data <- subset(PPO_data, PPO_mgProtein < 0.0081)
hist(PPO_data$PPO_mgProtein, breaks = 10, main = "PPO")
#normality testing
shapiro.test(PPO_data$PPO_mgProtein)
#p<0.05 data are not normal
PPO_data$PPO_transformed <- (PPO_data$PPO_mgProtein)^0.3
hist(PPO_data$PPO_transformed, breaks = 10, main = "PPO")
shapiro.test(PPO_data$PPO_transformed)
#normal following transformation above

#--------CATALASE-----------
hist(data$CAT_mgProtein, breaks = 10, main = "catalase")
rosnerTest(data$CAT_mgProtein, k=15)
#11 outliers in catalase filtered out
CAT_data <- subset(data, CAT_mgProtein > 0)
CAT_data <- subset(CAT_data, CAT_mgProtein < 2001)
hist(CAT_data$CAT_mgProtein, breaks = 10, main = "catalase")
shapiro.test(CAT_data$CAT_mgProtein)
#p<0.05, data are not normal
CAT_data$CAT_transformed <- (CAT_data$CAT_mgProtein)^0.01
hist(CAT_data$CAT_transformed, breaks = 10, main = "catalase")
shapiro.test(CAT_data$CAT_transformed)
#data normal following transformation above

#--------MELANIN-----------
hist(data$mgMelanin_mgTissue, breaks = 10, main = "melanin")
rosnerTest(data$mgMelanin_mgTissue, k=5)
#11 outliers in melanin filtered out
Mel_data <- subset(data, mgMelanin_mgTissue < 15)
hist(Mel_data$mgMelanin_mgTissue, breaks = 10, main = "melanin")
shapiro.test(Mel_data$mgMelanin_mgTissue)
#p<0.05, data are not normal
Mel_data$Melanin_transformed <- sqrt(Mel_data$mgMelanin_mgTissue)
hist(Mel_data$Melanin_transformed, breaks = 10, main = "melanin")
shapiro.test(Mel_data$Melanin_transformed)
#Data are normal following transformation above

#Subset data frames to only normal immune parameter columns and then use merge command to assemble data frame for principle components
PPO_data <- PPO_data[-c(1:2,4:7)]
CAT_data <- CAT_data[-c(1:2,4:7)]
Mel_data <- Mel_data[-c(1:2,4:7)]
merged_data <- merge(PPO_data, CAT_data, by = "row.names")
row.names(merged_data) <- merged_data$Row.names
merged_data <- merged_data[-c(1,5:6)]
merged_data <- merge(merged_data, Mel_data, by = "row.names")
row.names(merged_data) <- merged_data$Row.names
merged_data <- merged_data[-c(1:3)]
merged_data$Melanin_transformed <- replace_na(merged_data$Melanin_transformed, 0)

pca_data <- merged_data[c(1,2,5)]
#run principle components analysis
pca <- prcomp(pca_data, center = TRUE, scale. = TRUE)
summary(pca)

#merge metadata with PCA data
metadata_merged <- merge(merged_data, pca_data, by = "row.names")

#plot PCA
pca_plot <- autoplot(pca, data = metadata_merged, color = "dayCat.x")
pca_plot + theme_bw() + theme(panel.grid = element_blank())

#set up data for RDA (PCA data will be used from above, metadata will be filtered to variables of interest)
row.names(metadata_merged) <- metadata_merged$Row.names
metadata_merged <- metadata_merged[-c(1)]
RDA_meta <- metadata_merged[c(3,9:12)]
#remove NA values from designated columns
RDA_meta <- na.omit(RDA_meta[c("day.x","symbiontDensity","PPO_mgProtein","CAT_mgProtein","mgMelanin_mgTissue")])
pca_data <- na.omit(pca_data[c("CAT_transformed","PPO_transformed","Melanin_transformed")])

#run RDA analysis
rda <- rda(rda_data ~ ., data = RDA_meta)
summary(rda)
#selecting statistically significant variables
forwardSel <- ordiR2step(rda(rda_data ~ 1, data = RDA_meta), scope = formula(rda), direction = "forward",
                         R2scope = TRUE, pstep = 1000, trace = FALSE)
#see what variables are retained by forward call
forwardSel$call
#write new model
rda_sig <- rda(rda_data ~ mgMelanin_mgTissue + day.x, data = RDA_meta)
RsquareAdj(rda_sig)
#significance testing
anova.cca(rda_sig, step = 1000)
anova.cca(rda_sig, step = 1000, by = "term")
anova.cca(rda_sig, step = 1000, by = "axis")

#rda plotting
rda_plot <- ordiplot(rda_sig, scaling = 2, type = "text")

rda_coords <- scores(rda_sig, display = "sites", choices = c(1,2), scaling = 1)
rda_arrows <- scores(rda_sig, display = "bp", choices = c(1,2), scaling = 1)
rda_plotting <- merge(rda_coords, metadata_merged, by = "row.names")

rda_plot <-  ggplot() + geom_segment(data = rda_arrows, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), arrow = arrow(length = unit(0.1, "inches"))) +
  geom_point(data = rda_plotting, size = 3, pch = 21, aes(RDA1, RDA2, fill = treatment)) +
  theme_bw() + theme(panel.grid = element_blank()) + scale_fill_manual(values = c('control' = "#2252aa",
                                                                                  'heat' = "#ea4d33")) +
  ylim(-0.3, 0.3) + xlab("RDA1 (93.49%)") + ylab("RDA2 (<1%)")
rda_plot

