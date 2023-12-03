# Generating RF based ML models for predicting total ARG abundance, ARG type/sybtype count and distribution pattern 

#1.Data preprocessing 
dataset = read.csv('phylum_ML.csv')
dataset = dataset[,2:39]



#2.Splitting the data into training set and test set
library(caTools)
set.seed(123)
split = sample.split(dataset$total_ARGs, SplitRatio = 2/3)
training_set = subset(dataset, split == TRUE)
test_set = subset(dataset, split == FALSE)

#3. fitting regression to the dataset
#install.packages('randomForest')
library(randomForest)
set.seed(1234)
regressor = randomForest(x= training_set[38],
                         y= training_set$total_ARGs,
                         ntree = 500)


#4. predicting a new result with random forest regression 
y_pred = predict(regressor, test_set[38])
y_pred
y_pred_train  = predict(regressor, training_set[38])
y_pred_train
RF = randomForest(total_ARGs ~., data = training_set, ntrees = 1000, keep.forest=FALSE, importance = TRUE)

#5. Important variables, R squared 
library(Metrics)

mae(test_set$total_ARGs,y_pred)

library(caret)

postResample(y_pred, test_set$total_ARGs)['RMSE']^2 

postResample(y_pred , test_set$total_ARGs)['Rsquared']

varImpPlot(RF, type =2)

importance(RF)

#5. Predicted value Vs. Actual value plot 

library(ggplot2)
library(ggpubr)
sc_plot <- xle1 %>%
  ggplot(aes(x=Actual_Value, y=Value, color=Model))+
  geom_point(size= 2.5)+
  theme_minimal()+theme(axis.text.x = element_text(angle=0, size=12, hjust= 0.5, vjust= 0.5, face="bold"), 
                        panel.background = element_rect(fill="white", colour="gray", linewidth =0.2),
                        axis.title.y = element_text(face="bold", size = 14), 
                        axis.title.x = element_text(face="bold", size = 14),
                        axis.text.y = element_text(size=12, face ="bold"), panel.grid = element_blank(), 
                        legend.text = element_text(size =12, face ="bold"), 
                        legend.title = element_text(size =12, colour = "white"), 
                        plot.title = element_text(size = 17, face = "bold"))+
  labs(x="Actual Value")+ 
  labs(y="Predicted Value")+  scale_color_discrete(guide = guide_legend(override.aes = list(shape = 19)))+
  labs(title= "Total ARG abundance")

sc_plot +
  geom_smooth(method="lm", se = FALSE, formula = y ~ x)+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.x.npc = 0.01, label.y.npc =1,
           size = 4.2, fontface = "bold")
		   
# Generating heatmaps to show ARG abundance and distribution pattern
library(tidyverse)
library(pheatmap)
library(matrixStats)
library(readxl)
library(dplyr)
library(xlsx)
data <- read_excel("D:/Desktop/machine learning/Bracken_genuses/nmds.xlsx", sheet = "arg") 

data %>%
  select(1:18) %>%
  column_to_rownames("site") -> heatmap
heatmap - rowMeans((heatmap)) -> heatmap_meanSubtract

heatmap_meanSubtract %>%
  pheatmap()

heatmap_meanSubtract/rowSds(as.matrix(heatmap)) ->
  heatmap_zscores

read_excel("D:/Desktop/machine learning/Bracken_genuses/nmds.xlsx", sheet = "annotation_col") ->
  col_ann
col_ann %>% column_to_rownames("site")->col_ann
annotation_col <- data.frame(Type = col_ann$type)
rownames(annotation_col)<-colnames(heatmap_zscores)

# Generaring Bray-Curtis based NMDS plots based on bacterial composition and ARG profiles

library(readxl)
library(ggplot2)
library(wesanderson)
library(vegan)
pc <- read_excel("D:/Desktop/machine learning/Bracken_genuses/nmds.xlsx", sheet = "tax")
com = pc[, 3:ncol(pc)]
m_com =as.matrix(com)
set.seed(123)
nmds = metaMDS(m_com, distance = "bray")
nmds
plot(nmds)
data.scores = as.data.frame(scores(nmds)$sites)
data.scores$type = pc$type
data.scores$site = pc$site
head(data.scores)


ggplot() + 
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=type,colour=site),size=4) + # add the point markers
  scale_colour_manual(values = wes_palette(n=37, name="Darjeeling1", type = "continuous")) +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 9, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "type", y = "NMDS2", shape = "site")
  
