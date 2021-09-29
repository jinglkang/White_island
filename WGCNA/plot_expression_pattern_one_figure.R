library(ggpubr)

ggarrange(Blenny_pH_pos, Blenny_pH_neg, 
          Common_pH_pos, Common_pH_neg, 
          Blue_eyed_pH_pos, Blue_eyed_pH_neg,
          nrow=3, ncol=2,common.legend = T,align = "h")

ggarrange(Blenny_Length_pos, Blenny_Length_neg, 
          Common_Length_pos, Common_Length_neg, 
          Blue_eyed_Length_pos, Blue_eyed_Length_neg,
          Yaldwyn_Length_pos, Yaldwyn_Length_neg, 
          nrow=4, ncol=2,common.legend = T,align = "h")

ggarrange(Blenny_Salinity_pos, Blenny_Salinity_neg, 
          Common_Salinity_pos, Common_Salinity_neg, 
          Blue_eyed_Salinity_pos, Blue_eyed_Salinity_neg,
          Yaldwyn_Salinity_pos, Yaldwyn_Salinity_neg, 
          nrow=4, ncol=2,common.legend = T,align = "h")
