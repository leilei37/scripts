#twas整理
tg<-read.csv('~/Documents/twasgene.csv',header = F)
udg<-as.character(unique(tg[which(duplicated(tg[,1])),1]))
op<-data.frame(matrix(nrow = 1,ncol = 2))

i<-1
tg[intersect(which(duplicated(as.character(tg[,1]))),which(duplicated(tg[,2]))),]
for (i in 1:length(udg)) {
  op<-rbind(op,tg[which(as.character(tg[,1])==udg[i]),])
}
unique(tg[which(duplicated(tg[,1])),1])
#eqtl整理
load('~/Downloads/LEI/data/7p.RData')
#upset
library(UpSetR)
load('~/Documents/7ti/sum.RData')
listInput <- list(RT = grootmaf, SH = gshootmaf, LB =l3basemaf,LT =l3tipmaf,LD=ldmaf,LN=lnmaf,KN=kernmaf)
colorPalette<-c("#e41a1c","#377eb8","#4daf4a","9ecae1","#6baed6","#4292c6")
upset(fromList(listInput), nsets = 7,order.by = "degree",
      mb.ratio=c(0.7, 0.3),
#      set_size.scale_max = 500,
#      set_size.show = TRUE,
      nintersects = 150,color.pal = 1,
      keep.order = TRUE,
#      matrix.color ="#b35806", main.bar.color = colorPalette,
#      sets.bar.color = c("#e41a1c","#377eb8","#4daf4a"), 
      )
load('~/Documents/图数据/allpart2.RData')
length(grootcount[,1])
length(unique(grootcount[,1]))
listInput <- list(RT = grootcount[,1], SH = gshootcount[,1], LB =l3basecount[,1],LT =l3tipcount[,1],LD=ldcount[,1],LN=lncount[,1])
upset(fromList(listInput), nsets = 6,order.by = "freq",
      mb.ratio=c(0.7, 0.3),
      #      set_size.scale_max = 500,
      #      set_size.show = TRUE,
      nintersects = 110,color.pal = 1,
      keep.order = TRUE,
      #      matrix.color ="#b35806", main.bar.color = colorPalette,
      #      sets.bar.color = c("#e41a1c","#377eb8","#4daf4a"), 
)
#"degree"
listInput <- list(RT = grootsnpsum1[,7], SH = gshootsnpsum1[,7], LB =l3basesnpsum1[,7],LT =l3tipsnpsum1[,7],LD=ldsnpsum1[,7],LN=lnsnpsum1[,7])
upset(fromList(listInput), nsets = 6,order.by = "degree",
      mb.ratio=c(0.7, 0.3),
      #      set_size.scale_max = 500,
      #      set_size.show = TRUE,
      nintersects = 100,color.pal = 1,
      keep.order = TRUE,
      #      matrix.color ="#b35806", main.bar.color = colorPalette,
      #      sets.bar.color = c("#e41a1c","#377eb8","#4daf4a"), 
)
#柱状图重制 百分比
library(ggplot2)
eg<-read.csv('~/Documents/图数据/表达基因.csv')
library(forcats)
library(dplyr)
eg %>%
  mutate(GeneType = fct_relevel(GeneType,'Genes expressed in 7 tissues','Tissue specific genes','Genes expressed in 2-6 tissues')) %>%
#gene.type <-factor(gene.type,levels=c("common gene","other gene","tissue specific gene")) 
ggplot(aes(x=Tissues, y=count)) +
geom_col(aes(fill=GeneType),position = position_fill(reverse = TRUE))+
coord_flip()+theme_minimal()+
  geom_text(aes(label=count)
            , color="white", size=2,position=position_fill(0.5))+ylab("")+xlab("")+
  scale_fill_manual(breaks = c('Genes expressed in 7 tissues','Tissue specific genes','Genes expressed in 2-6 tissues'), 
#                    values = c('#c2a5cf','#9970ab','#762a83'))+
values = c('#0575B4','#58B4E9','#D8A4C0'))+
  scale_colour_pander()+
  theme(legend.position="left",legend.key.size = unit(10, "pt"))
###########
egenebai<-read.csv('~/Documents/图数据/egene.csv')
#write.csv(egenebai,file = '~/Documents/图数据/egene.csv')
egenebai %>%
  mutate(eGeneType = fct_relevel(eGeneType,"cis-eGene","trans-eGene","both-eGene")) %>%
ggplot(aes(x=tissue, y=count)) +
  geom_col(aes(fill=eGeneType),position = position_fill(reverse = T))+
  scale_fill_grey(start=0.3, end=0.6) +
  coord_flip()+
  theme_minimal()+
  geom_text(aes(label=count)
            , color="white", size=2.5,position=position_fill(0.5))+ylab("")+xlab('')+
  scale_fill_manual(breaks = c("cis-eGene","trans-eGene","both-eGene"), 
                    values = c('#35978f','#D2B48C','#A27333'))+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  theme(axis.text.x = element_text(size = 5))
#80cdc1
#35978f
#01665e
  scale_y_continuous(breaks=seq(0,1,0.1))+
  theme(axis.text.x = element_text(size = 5))
#  theme(axis.text.x = element_text(size = 5, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))
###########
eqtlbai<-read.csv('~/Documents/图数据/eQTL.csv')
eqtlbai %>%
  mutate(eQTLType = fct_relevel(eQTLType,"cis-eQTL","trans-eQTL")) %>%
  ggplot(aes(x=tissue, y=count)) +
  geom_col(aes(fill=eQTLType),position = position_fill(reverse = T))+
  #scale_fill_grey(start=0.3, end=0.6) +
  coord_flip()+theme_minimal()+
  geom_text(aes(label=count)
            , color="white", size=2,position=position_fill(0.5))+ylab("")+xlab('')+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  theme(axis.text.x = element_text(size = 5))+
    scale_fill_manual(breaks = c("cis-eQTL","trans-eQTL"), 
                      values = c('#D2B48C','#A27333'))
eqtlbai %>%
  mutate(eQTLType = fct_relevel(eQTLType,"cis-eQTL","trans-eQTL")) %>%
  ggplot(aes(x=tissue, y=count,fill=eQTLType)) +
  theme_classic()+
  theme(legend.position = c(0.83,0.5))+
  geom_bar(stat = 'identity', position = 'dodge')+coord_flip()+ylab("")+xlab('')+theme(legend.title=element_blank())+
  scale_fill_manual(breaks = c("cis-eQTL","trans-eQTL"), 
                    values = c('#D2B48C','#A27333'))
####################################¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥
###
eff<-read.csv('~/Documents/图数据/effect.csv')
colnames(eff)[2]<-'effect.type'
#eff %>%
#  mutate(eQTL.type = fct_relevel(eQTL.type,"cis-eQTL","trans-eQTL")) %>%
  ggplot(eff,aes(x=tissue, y=count)) +
  geom_col(aes(fill=effect.type),position = position_fill(reverse = T))+
  scale_fill_grey(start=0.3, end=0.6) +coord_flip()+theme_minimal()+
  geom_text(aes(label=count)
            , color="white", size=2,position=position_fill(0.5))+ylab("")+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  theme(axis.text.x = element_text(size = 5))
########效应值
  library(ggplot2)
  eff<-read.csv('~/Documents/图数据/effect百分比.csv')
  library(forcats)
  library(dplyr)
  eff %>%
    mutate(type = fct_relevel(type,"all positive effect","all negative effect","opposite effect")) %>%
    #gene.type <-factor(gene.type,levels=c("common gene","other gene","tissue specific gene")) 
    ggplot(aes(x=tissue, y=count)) +
    geom_col(aes(fill=type),position = position_fill(reverse = TRUE))+
    scale_fill_manual(
      breaks = c("all positive effect","all negative effect","opposite effect"), 
      values = c("#745550","#6A979A", "#B4D8EE")
    )+
    coord_flip()+theme_minimal()+
    geom_text(aes(label=count)
              , color="white", size=2,position=position_fill(0.5))+ylab("")+theme(legend.position="bottom",legend.title=element_blank(),legend.key.size = unit(10, "pt"))
#5ti4ti
  library(ggplot2)
  eff<-read.csv('~/Documents/图数据/4ti5tieffect.csv')
  eff<-read.csv('~/Documents/图数据/5tieffect.csv')
  eff<-read.csv('~/Documents/图数据/4tieffect.csv')
  library(forcats)
  library(dplyr)
  install.packages('ggrepel')
  library(ggrepel)
  eff %>%
    mutate(type = fct_relevel(Type,"All Positive Effect","All Negative Effect","Opposite Effect")) %>%
    #gene.type <-factor(gene.type,levels=c("common gene","other gene","tissue specific gene")) 
    ggplot(aes(x=Tissues, y=Counts,fill=Type)) +
    geom_bar(stat="identity") +
#    geom_col(aes(fill=Type))+
    scale_fill_manual(
      breaks = c("All Positive Effect","All Negative Effect","Opposite Effect"), 
      values = c("#745550","#6A979A", "#B4D8EE")
    )+
    coord_flip()+theme_minimal()+
    geom_text(aes( x=eff$Tissues,y=eff$label_y,label=eff$Counts)
              , color="white", size=2,hjust=1.5)+ylab("")+theme(legend.position="bottom",legend.title=element_blank(),legend.key.size = unit(10, "pt"))
#    geom_label_repel(aes(label=Counts),data=eff)
#108effect
  library(ggplot2)
  eff<-read.csv('~/Documents/图数据/108effect.csv')
  library(forcats)
  library(dplyr)
  eff %>%
    mutate(type = fct_relevel(type,"Positive Effect","Negative Effect")) %>%
    #gene.type <-factor(gene.type,levels=c("common gene","other gene","tissue specific gene")) 
    ggplot(aes(x=tissue, y=count)) +
    geom_col(aes(fill=type),position = position_fill(reverse = TRUE))+
    scale_fill_manual(
      breaks = c("Positive Effect","Negative Effect"), 
      values = c("#6A979A", "#B4D8EE")
    )+
    coord_flip()+theme_minimal()+
    geom_text(aes(label=count)
              , color="white", size=2,position=position_fill(0.5))+ylab("")+xlab('')+theme(legend.position="bottom",legend.title=element_blank(),legend.key.size = unit(10, "pt"))
  
#饼状图

library(scales)
  class(egp[,2])
  egp[,2]<-as.numeric(egp[,2])
egp<-read.csv('~/Documents/图数据/表达基因饼.csv')
egp<-egp[-3,]
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
attach(egp)
egp %>%
#  mutate(GeneType=fct_reorder('Genes expressed in 1-6 tissues','Conserved genes')) %>%
ggplot(aes(x="", y=egp$count, fill=egp$GeneType))+
geom_col(width = 1)+
coord_polar("y", start=0)+
  scale_fill_manual(breaks = c('Conserved genes','Genes expressed in 1-6 tissues'), 
                    values = c('#a6dba0','#9970ab'))+
  #scale_fill_grey(start=0.8, end=0.9)+
  blank_theme +guides(fill=guide_legend(title='GeneType'))+
  theme(axis.text.x=element_blank()) +
  geom_text(aes(y = egp$count/2 + c(0, cumsum(egp$count)[-length(egp$count)])),x = sum(egp$count)/24399,
                label = paste0(label_percent()(egp$count/24399),',',egp$count), size=4)
########nouse
class(egp$count)
label_value <- paste(round(egp$count/sum(egp$count) * 100, 0), '%', sep = '')
label <- paste(label_value,egp$count,sep = ',')
p <- ggplot(data = egp, mapping = aes(x = 'Content', y = count, fill = GeneType)) + geom_bar(stat = 'identity', position = 'stack',width = 1)
p + 
#  coord_polar(theta = 'y')+
  labs(x = '', y = '', title = '')+theme(axis.text = element_blank())+theme(axis.ticks = element_blank())+
#  geom_text(aes(y = df$nums/2 + c(0, cumsum(df$nums)[-length(df$nums)]), x = sum(df$nums)/150, label = label)
  geom_text(aes(y =  c(3, cumsum(egp$count)[-length(egp$count)]),x = sum(egp$count)/24399,
            label = label),size=4,position = position_dodge(0.9))
#############
#折线图
egene<-read.csv('~/Documents/图数据/egene折线图.csv',header = F)
egene<-read.csv('~/Documents/图数据/eqtl折线图.csv',header = F)

ggplot(data = egene, mapping = aes(x = V1, y = V2)) + geom_line()+
  scale_y_continuous(breaks=seq(0,0.45,0.05))+ylab("")+xlab('')+
  scale_x_continuous(breaks=seq(1,42,2),limits = c(1,42),expand=c(0,0))+
#xlim(1, 42)+
#  scale_x_discrete(expand=c(0,0))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

###############
#egene for article
eg<-read.csv('~/Documents/图数据/egenemap1ti.csv')
eg<-eg[-13,]
eg[6,1]<-'eGenes mapped in 1 tissue'
eg[10,1]<-'eGenes mapped in >1 tissues'
library(forcats)
library(dplyr)
factor(eg$x2)
eg$x1<-factor(eg$x1,levels = c('All tested genes','eGenes mapped in 1 tissue','eGenes mapped in >1 tissues'))
ggplot(data = eg, mapping = aes(x = x1, y = counts, fill = x2)) + geom_bar(stat = 'identity', position = 'fill')
library(stringr)
strcut<-function(x)  str_replace(x, "(.{13})", "\\1\n")
eg<-eg[-c(5,9),]
eg %>%
  mutate(x2 = fct_relevel(x2,"no hits","cis","trans","both")) %>%
ggplot(mapping = aes(x = x1, y = counts,fill= x2)) +
  theme_classic()+
  theme(legend.position = c(0.8,0.5))+
  geom_bar(stat = 'identity', position = 'dodge') + 
  geom_text(mapping = aes(label = counts), size = 3, colour = 'black', vjust = -1, hjust = .5, position = position_dodge(0.9))+
#  theme(legend.justification = c("right", "top"),legend.key.size = unit(10, "pt"),legend.title=element_blank(),
#        axis.text.x = element_text(size = 7, vjust = 0.5, hjust = 0.2, angle = 45),
#        axis.text.y = element_text(size = 5))+
  guides(fill=guide_legend(title=NULL))+
#  theme_minimal()+
  scale_fill_manual(breaks = c("no hits","cis","trans","both"), 
                    values = c('#C2C0A6','#A8B545','#6A8C69','#401E17'))+
  scale_x_discrete(breaks=c('All tested genes','eGenes mapped in 1 tissue','eGenes mapped in >1 tissues'), labels = c('All tested genes',strcut(c('eGenes mapped in 1 tissue','eGenes mapped in >1 tissues'))))+
  labs(x = "",y='Counts')
c('All tested genes','eGenes mapped in 1 tissue','eGenes mapped in >1 tissues')
###############
eg<-read.csv('~/Documents/图数据/eqtlmap.csv')
eg$x1<-factor(eg$x1,levels = c('eQTLs mapped in 1 tissue','eQTLs mapped in >1 tissues'))
eg[,1]<-'eQTLs mapped in 1 tissue'
eg[4:6,1]<-'eQTLs mapped in >1 tissues'
eg %>%
  mutate(x2 = fct_relevel(x2,"cis","trans","both")) %>%
  ggplot(mapping = aes(x = reorder(x1, -counts), y = counts,fill= x2)) +
  theme_classic()+
  theme(legend.position = c(0.87,0.5))+
  geom_bar(stat = 'identity', position = 'dodge') + 
  geom_text(mapping = aes(label = counts), size = 3, colour = 'black', vjust = -1, hjust = .5, position = position_dodge(0.9))+
  #  theme(legend.justification = c("right", "top"),legend.key.size = unit(10, "pt"),legend.title=element_blank(),
  #        axis.text.x = element_text(size = 7, vjust = 0.5, hjust = 0.2, angle = 45),
  #        axis.text.y = element_text(size = 5))+
  guides(fill=guide_legend(title=NULL))+
  #  theme_minimal()+
  scale_fill_manual(breaks = c("cis","trans","both"), 
                    values = c('#A8B545','#E3C75F','#CC8D1A'))+
  scale_x_discrete(breaks=c('eQTLs mapped in 1 tissue','eQTLs mapped in >1 tissues'), labels = c(strcut(c('eQTLs mapped in 1 tissue','eQTLs mapped in >1 tissues'))))+
  labs(x = "",y='Counts') 
#############egene freq
n_occur <- data.frame(table(l3basesnpsum3$V2))
n_occur <- data.frame(table(ldsnpsum3$V2))
n_occurld<-data.frame(table(n_occur[,2]))
egenefreq<-rbind(n_occurrt,n_occursh,n_occurlb,n_occurlt,n_occurlb,n_occurln)
write.csv(egenefreq,file = '~/Documents/图数据/egenefreq.csv')
egenefreq<-read.csv('~/Documents/图数据/egenefreq.csv')
factor(egenefreq$Tissue)
ggplot(data = egenefreq, mapping = aes(x = Var1, y = Freq, colour = Tissue))+ 
  geom_line(key_glyph = draw_key_rect) + geom_point()+#绘制线图和点图
  scale_x_continuous(breaks = c(1:6))+xlab('Numbers')+ ylab('Frequency')+
  theme_classic()+
  scale_color_manual(values = c('#401E17','#276573','#65BF8C','#F2B33D','#F2C1AE','#D95032'))+
  guides(colour=guide_legend(title=NULL))+
  theme(legend.position = 'left')
#  scale_y_continuous(breaks=1:3000, minor_breaks=NULL)+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
#        panel.border = element_blank(), #边框
#        axis.title = element_blank(),  #轴标题
#        axis.text = element_blank(), # 文本
        axis.ticks = element_blank())
  theme_bw()
  theme_minimal()
  theme_classic()
  theme_bw()
scale_linetype_manual(values = c(1,2)) #自定义线条类型+
scale_color_manual(values = c('steelblue','darkred')) #自定义颜色+ 
scale_shape_manual(values = c(21,23)) #自定义点形状+
scale_fill_manual(values = c('red','black')) #自定义点的填充色
##PVE
pve<-read.csv('~/Documents/图数据/pve.csv')
library(ggplot2)
p <- ggplot(pve,aes(x=reorder(Traits,X),y=PVE,fill=type))+geom_bar(position="dodge",stat="identity")
p+xlab("") + ylab("PVE") + labs(fill="")+theme_classic()+ scale_fill_manual(values = c('#7972A6','#F2B84B'))+
  theme(legend.position = 'bottom',axis.text.x=element_text(vjust=1,size=5,face = "bold"))

