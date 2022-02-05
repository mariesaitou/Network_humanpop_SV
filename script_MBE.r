library(ggplot2)
library("ggrepel")
library(rworldmap)
library(reshape2)
library(dplyr)
library(gplots)
library (tidyr)


## Figure 2B
ggplot(dat, aes((Allele.freq), y=log2(d_from_ave),alpha=d_from_ave,shapec=as.factor(kind4),color=as.factor(kind4))) +
  geom_point()+theme_bw() +geom_label_repel(data=MDSgene,aes(label = gene.name),  size =3.5, alpha=0.9,fontface = 'bold')

## Figure 3A
ggplot(hist2,aes(value,color=variable))+geom_histogram()+theme_bw() 

## Figure 3B
x<-MDS
x <- subset (x, kind2=="SV")

x<- x[order(x$Chromosome),]
nCHR <- length(unique(x$Chromosome))
x$BPcum <- NA
s <- 0
nbp <- c()

for (i in unique(x$Chromosome)){
  nbp[i] <- max(x[x$Chromosome == i,]$start)
  x[x$Chromosome == i,"BPcum"] <- x[x$Chromosome == i,"start"] + s
  s <- s + nbp[i]
}


axis.set <- x %>%
  group_by(Chromosome) %>%
  summarize(center = (max(BPcum) + min(BPcum)) / 2)
MDS<-subset(MDS, MDS$kind2=="SV")
quantile(MDS$d_from_ave, c(.10,.50,.90, .95, .99))
x1 <- subset (x, d_from_ave>4.23 )

mannplot <- ggplot(x, aes(x = BPcum, y = d_from_ave,
                          color = as.factor(Chromosome))) +
  geom_point(alpha = 0.75) +
  geom_hline(yintercept = 4.23, linetype = "dashed",color='coral', size=1) +
  scale_x_continuous(label = axis.set$Chromosome, breaks = axis.set$center)  +
  scale_color_manual(values = rep(c("gray", "black"), nCHR)) +
  scale_size_continuous(range = c(0.5,3))  +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
  )+geom_label_repel(data=x1,aes(label = gene.name), color = 'blue',
                     size = 3, alpha=0.8)
mannplot


## Figure 4D

input <- read.csv("sv.network.heatmap.csv", row.names=1,header = TRUE)
input1 <- as.matrix(input)
library(RColorBrewer)
my_palette <- rev(brewer.pal (8,"RdBu"))
heatmap.2 (input1, density.info="none", trace="none", col=my_palette, scale="column",
           Rowv = T, Colv=F,dendrogram = c("row"),cexRow=0.8, cexCol=1,  main="")

#############Alkan analysis

###prep Alkan data
a <- read.csv ("MDS_result_alkan2.csv")
b <- read.csv ("YRI.population_genes.csv")
b <- b [,c(1:5, 114,115)]
df <- merge (a, b, by.x="X", by.y="GENE")
#df <- subset (df, df$average<40)

df <- df[!grepl("^LINC", df$X),]
df <- df[!grepl("^MIR", df$X),]
df <- df[!grepl("^LOC", df$X),]
df <- df[!grepl("^SNORD", df$X),]
#install.packages("splitstackshape")
library (splitstackshape)
head(df)
df <- cSplit(df, "X", sep=";", stripWhite=TRUE)
df <- distinct(df, X_01, .keep_all = TRUE)
df <- subset (df, df$CHROM!="chrX")

df$CHROM<-gsub("chr","",as.character(df$CHROM))
df$chrom_order <- as.numeric(df$CHROM)

head(df)
### prep 1000 genomes data

a <- read.csv("STable.csv")
head(a)
a <- subset (a, a$d_from_ave>4.26)
m <- subset (a, a$gene.name!="" & a$kind!="SNP" & a$kind!="Known.SNP")
df$Presence_absence <- ifelse(df$X_01 %in% m$gene.name, 1, 0) # this does the trick

df$Presence_absence <- as.factor (df$Presence_absence)
df$CHROM = with(df, reorder(CHROM, -chrom_order))

aman <- subset(df, grepl(X_01, "AMY") )
aman <- df[grepl("^AMY", df$X_01),]
aman <- subset (df, df$MU>2)

png ("Canavar.png",  width = 7, height = 5, units = 'in', res = 300)
library (RColorBrewer)
library (ggridges)
ggplot(df, aes(x = MU, y = CHROM, label=X_01)) +
  geom_density_ridges_gradient( scale = 0.5, size = 0.5, rel_min_height = 0.01,
                                jittered_points=TRUE, position = "raincloud",
                                alpha = 0.5, point_size=0.8,point_color="#666666") +
  theme_bw()+
geom_text_repel(data= aman, max.overlaps = 100, size=3, 
                box.padding = unit(1, "lines"),
                point.padding = unit(0.2, "lines"))
dev.off()
?scale_colour_manual

###overlap ratios

a <- read.csv ("MDS_result_alkan2.csv")
b <- read.csv ("gene_with_protein_product.csv")
df <- subset (a, a$X %in% b$symbol)
#df <- a

c <- read.csv("STable.csv")
c <- subset (c, c$d_from_ave>4.26)
c <- subset (c, c$gene.name!="" & c$kind!="SNP" & 
               c$kind!="Known.SNP" & 
               c$kind!="ALU" & c$kind!="LINE1"&
             c$kind=="CNV")
c <- subset (c, c$gene.name %in% df$X)
c <- unique (c$gene.name)

quantile (a$MU, 0.50)

m <- subset (df, df$MU>0.163095)
n <- subset (c, c  %in% m$X)

######## Fst vs. Mu examples

install.packages ("rworldmap")
library (rworldmap)



df <- read.csv("Fst_mu_allele.csv")
df$LONG <- as.numeric (df$Population.longitude)
df$LAT <- as.numeric (df$Population.latitude)

df$Derived <- df$esv3643467
df$Ancestral <- df$pop_size - df$Derived

mapPies( df,nameX="LONG", nameY="LAT",
          nameZs=c('Ancestral', 'Derived'),mapRegion='world',
         zColours=c("light grey", "#D95F02"),
         landCol = "light grey")

######## distribution of allele counts

df <- read.csv ("allele.number.csv")
head(df)
hist(df$Alelle.number, col="#D95F02", main="Allele Number", xlab="Copy number of alleles")


######## comparision of Mu vs. allele frequency

df <- read.csv ("STable (1).csv")
head(df)
df <- subset (df, df$kind=="SNP")
ggplot(df, aes(x = Allele.freq, y = MU)) +
  geom_point(alpha=0.1) + geom_smooth() +
  theme_classic()
  

