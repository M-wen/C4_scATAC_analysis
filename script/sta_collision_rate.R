### Get the parameters
parser = argparse::ArgumentParser(description="Script to plot collision rate")
parser$add_argument('-I','--input', help='input file')
parser$add_argument('-O','--out', help='out directory')
args = parser$parse_args()

###get library
library(ggplot2)
library(cowplot)
library(data.table)
data=as.data.frame(fread(paste0(args$input),header = T))

###plot Number of beads per droplet
name=data[1]
name$Num=sapply(strsplit(as.character(name$CellBarcode),"N0"),"[",2)
table=as.data.frame(table(name$Num))
colnames(table)=c("Num","Count")

out_name_1=paste0(args$out,"/03.d2cfile/Plot4_Number_of_beads_per_droplet.svg")
#out_name_1="Plot4_Number_of_beads_per_droplet.svg"
svg(out_name_1,width = 5,height = 4)
#out_name_1=paste0(args$out,"/Plot4_Number_of_beads_per_droplet.svg")

p1=ggplot(data = table, mapping = aes(x = factor(Num), y = Count, fill = Num)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_brewer(palette = 'Set3',name = paste("Total",nrow(data),sep = " "),labels = paste(levels(table$Num)," ",table$Count))+
  theme_bw()+ theme_cowplot() + xlab("Number of beads per droplet")+
  theme(text = element_text(size=10),axis.text=element_text(size=8),plot.margin = unit(c(0.5, 0.2, 0.3, 0.2), "inches"))
#ggsave(file=out_name_1,plot=p1,width=5, height=4)
print(p1)
dev.off()

###plot Pllution collision rate after merge
#result=data[,c(1,3,7,9,14)]
result=data
result$hs_rate=result$uniqueHumanFrags/result$uniqueFrags
result$mm_rate=result$uniqueMouseFrags/result$uniqueFrags
result$map_rate_max=apply(result[,c("hs_rate","mm_rate")],1,max)
result$species <- as.factor(ifelse(result$map_rate_max > 0.9,
                                     ifelse(result$hs_rate > result$mm_rate ,'Human','Mouse'),'Mixed'))

pollution_rate=paste0(round(nrow(result[result$species=="Mixed",])/nrow(result),4)*100,"%")
a=paste("collision_rate",pollution_rate,sep=":")
write.table(a,paste0(args$out,"/03.d2cfile/1.collision_rate.tsv"),quote=F,col.names=F,row.names=F)

num=c(result$uniqueHumanFrags*2,result$uniqueMouseFrags*2)
num_sort=sort(num,decreasing = T)
num_sort=num_sort/10000

out_name_2=paste0(args$out,"/03.d2cfile/Plot3_Merged_Species_barcode_pollution.svg")
svg(out_name_2,width = 5,height = 4)

#col=c("Human"="#E41A1C","Mouse"="#377EB8","Mixied"="grey")
col=c("#E41A1C","#377EB8","grey")
result$species <- factor(result$species,levels= c("Human","Mouse","Mixed"))
Lables_y <- pretty(range(result$uniqueMouseFrags*2/10000),5)
Lables_x <- pretty(range(result$uniqueHumanFrags*2/10000),5)

p2=ggplot(data=result,aes(x=uniqueHumanFrags*2/10000, y =uniqueMouseFrags*2/10000,
                         color=species))+
  geom_point(alpha=0.9, size=0.5) +
  xlab(expression(Human~usable~reads~"(" %*% 10^4~")"))+
  ylab(expression(Mouse~usable~reads~"(" %*% 10^4~")")) +
  ggtitle(paste0("collision rate ",pollution_rate))+
  theme_gray() +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=10),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  scale_x_continuous(breaks = Lables_x,limits = c(0,num_sort[30]),labels=Lables_x)+
  scale_y_continuous(breaks = Lables_y,limits = c(0,num_sort[30]),labels=Lables_y)+
  scale_color_manual(values = col,labels = paste(table(result$species),"  ",levels(result$species)))+
  theme(panel.grid =element_blank()) + coord_fixed()
#ggsave(file=out_name_2,plot=p2,width=5, height=4)
print(p2)
dev.off()

###out put barcode list

sub_human=subset(result,result$species=="Human")[1]
sub_mm=subset(result,result$species=="Mouse")[1]
out_name_3=paste0(args$out,"/03.d2cfile/Human_barcode_list.txt")
out_name_4=paste0(args$out,"/03.d2cfile/Mouse_barcode_list.txt")
write.table(sub_human,out_name_3,sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(sub_mm,out_name_4,sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)

###split human and muouse QC stat
hg_data=data[match(sub_human$CellBarcode,data$CellBarcode),
                   c("CellBarcode","totalHumanFrags","uniqueHumanFrags","tssProportion")]
colnames(hg_data) <- c("CellBarcode","totalNuclearFrags","uniqueNuclearFrags",
                          "tssProportion")

mm_data=data[match(sub_mm$CellBarcode,data$CellBarcode),
                   c("CellBarcode","totalMouseFrags","uniqueMouseFrags","tssProportion")]
colnames(mm_data) <- c("CellBarcode","totalNuclearFrags","uniqueNuclearFrags",
                          "tssProportion")
write.csv(hg_data,file=paste(args$out,"/03.d2cfile/Human_QC_stat.csv",sep=""),quote=F,row.names=F)
write.csv(mm_data,file=paste(args$out,"/03.d2cfile/Mouse_QC_stat.csv",sep=""),quote=F,row.names=F)
###TSS QC
#result$QC <- as.factor(ifelse(result$tssProportion > 0.1 & log10(result$uniqueHumanFrags*2) >3,
#                              'TRUE','FALSE'))


#out_name_5=paste0(args$out,"/Summary.txt")
#write.table(result,out_name_5,sep = "\t",quote = FALSE)


sum_Human=subset(result,result$species=="Human")
sum_Mouse=subset(result,result$species=="Mouse")
sum_Human$QC <- as.factor(ifelse(sum_Human$tssProportion > 0.1 & log10(sum_Human$uniqueHumanFrags*2) >3,'TRUE','FALSE'))
sum_Mouse$QC <- as.factor(ifelse(sum_Mouse$tssProportion > 0.1 & log10(sum_Mouse$uniqueMouseFrags*2) >3,'TRUE','FALSE'))
sum_Human$QC <- factor(sum_Human$QC,levels=c("TRUE","FALSE"))
sum_Mouse$QC <- factor(sum_Mouse$QC,levels=c("TRUE","FALSE"))
plot=list()

plot[[1]]=ggplot(data=sum_Human,aes(x=log10(sum_Human$uniqueHumanFrags*2), y =tssProportion,
                                    color=QC))+
  geom_point(alpha=1, size=1.0) +
  xlab("log10 Human usable reads") +
  theme_gray() +
  theme_bw()+ geom_hline(aes(yintercept=0.1),color="black",linetype="dashed")+
  geom_vline(aes(xintercept=3.0),color="black",linetype="dashed")+
  scale_y_continuous(breaks =seq(0,0.5,0.1))+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),panel.grid = element_blank(),axis.line = element_line(colour = "black"))+
  scale_color_manual(values = c("#E41A1C","grey"),labels = paste(table(sum_Human$QC),"  ",levels(sum_Human$QC)))


plot[[2]]=ggplot(data=sum_Mouse,aes(x=log10(sum_Mouse$uniqueMouseFrags*2), y =tssProportion,
                                    color=QC))+
  geom_point(alpha=1, size=1.0) +
  xlab("log10 Mouse usable reads") +
  theme_gray() +
  theme_bw()+ geom_hline(aes(yintercept=0.1),color="black",linetype="dashed")+
  geom_vline(aes(xintercept=3.0),color="black",linetype="dashed")+
  scale_y_continuous(breaks =seq(0,0.5,0.1))+
  theme(plot.title = element_text(hjust = 0.5),panel.grid = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  scale_color_manual(values = c("#377EB8","grey"),labels = paste(table(sum_Mouse$QC),"  ",levels(sum_Mouse$QC)))

#out_name_4=paste0(args$out,"/TSS_QC.png")
out_name_4=paste0(args$out,"/03.d2cfile/Plot5_human_TSS_QC.svg")
svg(out_name_4,height = 4,width = 5)
#do.call(grid.arrange,c(plot,ncol=2))
print(plot[[1]])
dev.off()

 out_name_5=paste0(args$out,"/03.d2cfile/Plot6_mouse_TSS_QC.svg")
 svg(out_name_5,height = 4,width = 5)
 print(plot[[2]])
 dev.off()

