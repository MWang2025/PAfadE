#this is the code for DE analysis and for plots


library(MSnbase)
f <- "P1021_100823_Proteins.txt"
getEcols(f, split = "\t")
e <- c(70:83,85)
x <- readMSnSet2(f, e, sep = "\t")
sampleNames(x) <- c("Glc_1","Glc_2","Glc_3","Glc_4",
                    "C8_1","C8_2","C8_3","C8_4",
                    "C16_1","C16_2","C16_3","C16_4",
                    "C18_1","C18_2","C18_4")
tmp <- data.frame(do.call(rbind, strsplit(sampleNames(x), "_")))
names(tmp) <- c("Sample", "Replicate")
rownames(tmp) <- sampleNames(x)
pData(x) <- tmp
featureNames(x) <- fData(x)$Accession
x<- x[!grepl( "^cRAP" , fData(x)$Accession ) ]
dev.new ()
library(pRoloc)
naplot(x)
x <- filterNA(x, 0)
x <- filterZero(x,0)

x<- log(x,2)
y<-t(x)
xx<- featureNames(y)
fData(y)$markers <- fData(y)$Sample
plot2D(y, cex = 1.5)
highlightOnPlot(y, foi = xx, labels = TRUE, pos = 4, cex = .8)
addLegend(y, where = "topleft",horiz = T)
boxplot(exprs(x), las=2, col = c(rep(c("blue", "yellow","green", "red"), times = c(4,4,4,3))))

x <- normalise(x, "diff.median")

exprs(x) <- exprs(x) + 2
boxplot(exprs(x), las=2, col = c(rep(c("blue", "yellow","green", "red"), times = c(4,4,4,3))))

y<-t(x)
xx <- featureNames(y)
fData(y)$markers <- fData(y)$Sample
plot2D(y, cex = 1.5)
addLegend(y, where = "topright",horiz = F)


library(ggplot2)
library(limma)
library(plotly)
mat <- model.matrix(~ 0 + pData(x)$Sample)
colnames(mat) <- c('C16','C18','C8','Glc')
pData(x)$cmat <- mat

#C8 vs Glc
dostats_Glc_C8 <- function(x) {
  cont.matrix <- makeContrasts(DE = C8 - Glc,
                               levels = pData(x)[, "cmat"])
  fit <- lmFit(exprs(x), pData(x)[, "cmat"])
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  tt1 <- topTable(fit2, adjust.method = "BH", coef = "DE", number = Inf)
  return(tt1)
}
res_Glc_C8<- dostats_Glc_C8(x)
fData(x)$DE_Glc_C8 <- res_Glc_C8[featureNames(x), ]

dfr_Glc_C8 <- data.frame(Desc = fData(x)$Desc,
                      round(exprs(x),2),
                      log2_FC = fData(x)$DE_Glc_C8[, "logFC"],
                      adj.P.Val = fData(x)$DE_Glc_C8[, "adj.P.Val"])


#Glc_C18
dostats_Glc_C18 <- function(x) {
  cont.matrix <- makeContrasts(DE = C18 - Glc,
                               levels = pData(x)[, "cmat"])
  fit <- lmFit(exprs(x), pData(x)[, "cmat"])
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  tt1 <- topTable(fit2, adjust.method = "BH", coef = "DE", number = Inf)
  return(tt1)
}
res_Glc_C18<- dostats_Glc_C18(x)
fData(x)$DE_Glc_C18 <- res_Glc_C18[featureNames(x), ]

dfr_Glc_C18 <- data.frame(Desc = fData(x)$Desc,
                         round(exprs(x),2),
                         log2_FC = fData(x)$DE_Glc_C18[, "logFC"],
                         adj.P.Val = fData(x)$DE_Glc_C18[, "adj.P.Val"])


#Glc_C16
dostats_Glc_C16 <- function(x) {
  cont.matrix <- makeContrasts(DE = C16 - Glc,
                               levels = pData(x)[, "cmat"])
  fit <- lmFit(exprs(x), pData(x)[, "cmat"])
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  tt1 <- topTable(fit2, adjust.method = "BH", coef = "DE", number = Inf)
  return(tt1)
}
res_Glc_C16<- dostats_Glc_C16(x)
fData(x)$DE_Glc_C16 <- res_Glc_C16[featureNames(x), ]

dfr_Glc_C16 <- data.frame(Desc = fData(x)$Desc,
                          round(exprs(x),2),
                          log2_FC = fData(x)$DE_Glc_C16[, "logFC"],
                          adj.P.Val = fData(x)$DE_Glc_C16[, "adj.P.Val"])


#differential expression data of all three fatty acids vs. Glc
df_1021<-data.frame(dfr_Glc_C8,dfr_Glc_C16[,c(17,18)],dfr_Glc_C18[,c(17,18)])
colnames(df_1021)[17:22]<- c("Glc_C8 Log2FC", "Glc_C8 adj.P.Val", "Glc_C16 Log2FC","Glc_C16 adj.P.Val","Glc_C18 Log2FC","Glc_C18 adj.P.Val")

#Add protein symbols
library(stringr)
Function<-sub("OS=.*", "", df_1021$Desc)
Locus_tag16<- str_extract(dfr_Glc_C16$Desc, "(?<=GN=).*(?= PE)")
str_sub(Locus_tag16, 1, 1) <- str_sub(Locus_tag16, 1, 1) %>% str_to_upper()

df_1021<-data.frame(Locus_tag16,Function,df_1021[,(2:22)])
colnames(df_1021)[1]<-c("Protein_Symbol")

#correct the protein symbols - add abbraviation to some specific proteins 
df_1021$Protein_Symbol<-replace(df_1021$Protein_Symbol,df_1021$Protein_Symbol=="PA2634","ICL")
df_1021$Protein_Symbol<-replace(df_1021$Protein_Symbol,df_1021$Protein_Symbol=="PA2957","PvrA")
df_1021$Protein_Symbol<-replace(df_1021$Protein_Symbol,df_1021$Protein_Symbol=="FadB","FadB5")
df_1021$Protein_Symbol<-replace(df_1021$Protein_Symbol,df_1021$Protein_Symbol=="FadA","FadA5")

for(i in 1:nrow(df_1021)){
  if(rownames(df_1021)[i]=='Q9HXM8'){df_1021[i,1]<-'PA3767'}
  else if(rownames(df_1021)[i]=='Q9HW98'){df_1021[i,1]<-'PA4302'}
  else if(rownames(df_1021)[i]=='Q9HTC7'){df_1021[i,1]<-'PA5439'}
}


df_C8<-data.frame(df_1021[1:19])
df_C16<-data.frame(df_1021[1:17],df_1021[20:21])
df_C18<-data.frame(df_1021[1:17],df_1021[22:23])


#venn diagram
library(ggvenn)
#1.1 filter out the up-regulated proteins
df8_UP <- df_C8 %>% filter(Glc_C8.Log2FC>=1, Glc_C8.adj.P.Val<=0.01)
df16_UP <- df_C16 %>% filter(Glc_C16.Log2FC>=1, Glc_C16.adj.P.Val<=0.01)
df18_UP <- df_C18 %>% filter(Glc_C18.Log2FC>=1, Glc_C18.adj.P.Val<=0.01)
#1.2 venn diagaram of the upregulated protein
x<-list(Octanoate=rownames(df8_UP),Palmitate=rownames(df16_UP),Oleate=rownames(df18_UP))
ggvenn(x,c("Octanoate","Palmitate","Oleate"))
col1<-c("#007ED3","#FF9D1E","#894FC6")
ggvenn(x,c("Octanoate","Palmitate","Oleate"),
       show_percentage=FALSE,
       fill_color=col1,
       set_name_size=10, 
       text_size=10,
       stroke_linetype = "solid",
       stroke_size=0.6) 
#2.1 filter out the down-regulated proteins
df8_Dn <- df_C8 %>% filter(Glc_C8.Log2FC<=-1, Glc_C8.adj.P.Val<=0.05)
df18_Dn <- df_C18 %>% filter(Glc_C18.Log2FC<=-1, Glc_C18.adj.P.Val<=0.05)
df16_Dn <- df_C16 %>% filter(Glc_C16.Log2FC<=-1, Glc_C16.adj.P.Val<=0.05)

#2.2 plot veen diagram of the down-regulated proteins
y<-list(Octanoate=rownames(df8_Dn),Palmitate=rownames(df16_Dn),Oleate=rownames(df18_Dn))
ggvenn(y,c("Octanoate","Palmitate","Oleate"))
col1<-c("#007ED3","#FF9D1E","#894FC6")
ggvenn(y,c("Octanoate","Palmitate","Oleate"),
       show_percentage=FALSE,
       fill_color=col1,
       set_name_size=10, 
       text_size=10,
       stroke_linetype = "solid",
       stroke_size=0.6) 


#volcanol plot
library(EnhancedVolcano)
dev.new (width=12,height=15)
keyvars<-c('PA0506','PA0508')

##1.1 filter those non-significant proteins for grey color
df8_NS<-df_C8 %>% filter(abs(Glc_C8.Log2FC)<1, Glc_C8.adj.P.Val<=0.01)
df16_NS<-df_C16 %>% filter(abs(Glc_C16.Log2FC)<1, Glc_C16.adj.P.Val<=0.01)
df18_NS<-df_C18 %>% filter(abs(Glc_C18.Log2FC)<1, Glc_C18.adj.P.Val<=0.01)

##1.2 filter those with p value >0.01
df8_NC<-df_C8 %>% filter(Glc_C8.adj.P.Val>0.01)
df18_NC<-df_C18 %>% filter(Glc_C18.adj.P.Val>0.01)
df16_NC<-df_C16 %>% filter(Glc_C16.adj.P.Val>0.01)


#figure2a-Glc vs. C8

keyvals.C8<-ifelse(
  df_C8$Protein_Symbol %in% df8_UP$Protein_Symbol, 'red3',(ifelse
                                                           (df_C8$Protein_Symbol %in% df8_Dn$Protein_Symbol,'royalblue3',ifelse
                                                             (df_C8$Protein_Symbol %in% df8_NS$Protein_Symbol,'grey','gold'
                                                             )
                                                           )
  )
)

names(keyvals.C8)[keyvals.C8=='red3']<-'df8_UP$Protein_Symbol'
names(keyvals.C8)[keyvals.C8=='royalblue3']<-'df8_Dn$Protein_Symbol'
names(keyvals.C8)[keyvals.C8=='grey']<-'df8_NS$Protein_Symbol'
names(keyvals.C8)[keyvals.C8=='gold']<-'df8_NC$Protein_Symbol'

EnhancedVolcano(df_C8,lab = df_C8$Protein_Symbol,
                    xlab = bquote(~log[2] ~ "fold change"),
                    ylab = bquote(~-log[10] ~"adj."~ italic(p)~"value"),
                    x='Glc_C8.Log2FC',y='Glc_C8.adj.P.Val',
                    #title = bquote(bold(atop("C8:0 vs Glucose",""))),
                    title="C8:0 vs Glucose\n    ",
                    legendPosition = "none",
                    caption = NULL,
                    #col = c("grey30", "grey30", "royalblue", "red2"),
                    #title=expression(paste("Glc vs C18:1",Delta^9)),                                             
                    pCutoff = 0.01,
                    xlim = c(-6.5, 6.5),
                    ylim = c(0, -log10(10e-18)),
                    #borderWidth = 1.0,
                    FCcutoff = 1.0,
                    pointSize = c(ifelse(df_C8$Protein_Symbol %in% keyvars, 6, 3)),
                    #pointSize =3.0,
                    labSize = 4.5,
                    subtitle=NULL,
                    legendLabSize = 12,legendIconSize = 4.0,
                    labhjust = 0.5,
                    labvjust = 1.5,
                    boxedLabels = T,
                    #selectLab = c(Glclab,FAlab),
                    #shapeCustom = keyvals,
                    selectLab = c('Pgl','Eda','Edd','Zwf','Glk','GntR','GltR','PA4435','FadA5','FadB5','GlcB','EtfA','EtfB','FadD1','PA0506','PA0508',"ICL","PA1631"),
                    drawConnectors = T,
                    widthConnectors = 0.05,
                    colCustom = keyvals.C8
)

#Figure 2b Glc vs. C16

keyvals.C16<-ifelse(
  df_C16$Protein_Symbol %in% df16_UP$Protein_Symbol, 'red3',(ifelse
                                                             (df_C16$Protein_Symbol %in% df16_Dn$Protein_Symbol,'royalblue3',ifelse
                                                               (df_C16$Protein_Symbol %in% df16_NS$Protein_Symbol,'grey','gold'
                                                               )
                                                             )
  )
)

names(keyvals.C16)[keyvals.C16=='red3']<-'df16_UP$Protein_Symbol'
names(keyvals.C16)[keyvals.C16=='royalblue3']<-'df16_Dn$Protein_Symbol'
names(keyvals.C16)[keyvals.C16=='grey']<-'df16_NS$Protein_Symbol'
names(keyvals.C16)[keyvals.C16=='gold']<-'df16_NC$Protein_Symbol'

EnhancedVolcano(df_C16,lab = df_C16$Protein_Symbol,
                xlab = bquote(~log[2] ~ "fold change"),
                ylab = bquote(~-log[10] ~"adj."~ italic(p)~"value"),
                x='Glc_C16.Log2FC',y='Glc_C16.adj.P.Val',
                #title=expression(bold(paste("C18:1",Delta^9, " vs Glucose"),"\n", sep="")), 
                title = bquote(bold(atop("C16:0 vs Glucose",""))),
                #title="C8:0 vs Glucose\n    ",
                legendPosition = "none",
                caption = NULL,
                #col = c("grey30", "grey30", "royalblue", "red2"),
                
                pCutoff = 0.01,
                xlim = c(-6.5, 6.5),
                ylim = c(0, -log10(10e-18)),
                #borderWidth = 1.0,
                FCcutoff = 1.0,
                pointSize = c(ifelse(df_C16$Protein_Symbol %in% keyvars, 6, 3)),
                #pointSize =3.0,
                labSize = 4.5,
                subtitle=NULL,
                legendLabSize = 12,legendIconSize = 4.0,
                labhjust = 0.5,
                labvjust = 1.5,
                boxedLabels = T,
                #selectLab = c(Glclab,FAlab),
                #shapeCustom = keyvals,
                selectLab = c('Pgl','Eda','Edd','Zwf','Glk','GntR','GltR','PA4435','FadA5','FadB5','GlcB','EtfA','EtfB','FadD1','PA0506','PA0508',"ICL"),
                drawConnectors = T,
                widthConnectors = 0.05,
                colCustom = keyvals.C16
)
#Figure2c- Glc vs. C18
keyvals.C18<-ifelse(
  df_C18$Protein_Symbol %in% df18_UP$Protein_Symbol, 'red3',(ifelse
                                                             (df_C18$Protein_Symbol %in% df18_Dn$Protein_Symbol,'royalblue3',ifelse
                                                               (df_C18$Protein_Symbol %in% df18_NS$Protein_Symbol,'grey','gold'
                                                               )
                                                             )
  )
)
names(keyvals.C18)[keyvals.C18=='red3']<-'df18_UP$Protein_Symbol'
names(keyvals.C18)[keyvals.C18=='royalblue3']<-'df18_Dn$Protein_Symbol'
names(keyvals.C18)[keyvals.C18=='grey']<-'df18_NS$Protein_Symbol'
names(keyvals.C18)[keyvals.C18=='gold']<-'df18_NC$Protein_Symbol'

EnhancedVolcano(df_C18,lab = df_C18$Protein_Symbol,
                xlab = bquote(~log[2] ~ "fold change"),
                ylab = bquote(~-log[10] ~"adj."~ italic(p)~"value"),
                x='Glc_C18.Log2FC',y='Glc_C18.adj.P.Val',
                title=expression(bold(paste("C18:1",Delta^9, " vs Glucose"),"\n", sep="")), 
                #title = bquote(bold(atop("C8:0 vs Glucose",""))),
                #title="C8:0 vs Glucose\n    ",
                legendPosition = "none",
                caption = NULL,
                #col = c("grey30", "grey30", "royalblue", "red2"),
                #title=expression(paste("Glc vs C18:1",Delta^9)),                                             
                pCutoff = 0.01,
                xlim = c(-6.5, 6.5),
                ylim = c(0, -log10(10e-18)),
                #borderWidth = 1.0,
                FCcutoff = 1.0,
                pointSize = c(ifelse(df_C18$Protein_Symbol %in% keyvars, 6, 3)),
                #pointSize =3.0,
                labSize = 4.5,
                subtitle=NULL,
                legendLabSize = 12,legendIconSize = 4.0,
                labhjust = 0.5,
                labvjust = 1.5,
                boxedLabels = T,
                #selectLab = c(Glclab,FAlab),
                #shapeCustom = keyvals,
                selectLab = c('Pgl','Eda','Edd','Zwf','Glk','GntR','GltR','PA4435','FadA5','FadB5','GlcB','EtfA','EtfB','FadD1','PA0506','PA0508',"ICL"),
                drawConnectors = T,
                widthConnectors = 0.05,
                colCustom = keyvals.C18
)
