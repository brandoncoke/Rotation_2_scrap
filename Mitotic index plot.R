################################################################################
#Package install
################################################################################
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
library(ggplot2)
################################################################################
#Data import
################################################################################
mitotic_index= 
  read.csv("~//Rotation 2//Raw data//Microscopy//Mitotic index.csv",
           header = T)
colnames(mitotic_index)[1]= "Mitotic_stage"
mitotic_index$Mitotic_stage=  factor(mitotic_index$Mitotic_stage, levels=unique(mitotic_index$Mitotic_stage))
total=mitotic_index[1,3]
mitotic_index= mitotic_index[mitotic_index$Mitotic_stage != "Total",]

mitotic_index$Count= ((mitotic_index$Count/500)*100)
################################################################################
#Graph function creation
################################################################################
brandontheme=theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"),
  plot.title = element_text(color="black", size=14,
                                             face="bold",hjust = 0.5),
                   axis.text = element_text(color="black", size=12, 
                                            face="bold",hjust = 0.5),
                   axis.title = element_text(color="black", size=12, face="bold",
                                             hjust = 0.5))
ggplot(mitotic_index, aes(fill=Group, y=Count, x=Mitotic_stage)) + 
  geom_bar(position="dodge", stat="identity") +
  brandontheme +
  labs(y="Percentage of observed cells in stage",x="Stage of mitosis",title =
        "Breakdown of mitotic cell stage in different cell cultures") +
  scale_y_continuous(expand = c(0, 0.01), breaks =  seq(from=0,
                                                        to= max(mitotic_index$Count)+30,
                                                        by = 1)) +
  coord_flip() +
  scale_fill_manual(values=c("#990000", "#008080"))