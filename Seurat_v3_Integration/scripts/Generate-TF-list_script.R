




TF_list_2 <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/GO_lists/GO_DNA_binding.csv", head = F)
TF_list_3 <- read.csv("~/Downloads/GO_transcription_DNA-template.csv", head = F)
TF_list_4 <- read.csv("~/Downloads/GO_transcription_factor_activity.csv", head = F)

TF_list_5 <- read.csv("~/Downloads/GO_transcription_factor.csv", head = F)
TF_list_6 <- read.csv("~/Downloads/Danio_rerio_TF.txt", sep = "\t")

# This gives a df with 3141 unique entries for V2 (gene symbol)
TF_list <- Reduce(rbind, list(TF_list_3[,1:6], TF_list_4[,1:6], TF_list_5[,1:6], TF_list_6))

write.table(TF_list, file = "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/GO_lists/GO_TF_list.csv",col.names = F)
