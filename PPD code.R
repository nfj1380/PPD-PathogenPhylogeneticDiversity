#PPD code 

#Written by Nick Fountain-Jones (nfj@umn.edu)

library(ape)
library(picante)
library(pez)
library(treeio)
library(ggtree)
library(tidyverse)


#------------------------------------------------------------------------
##################Visualize MCC Tree##################
#------------------------------------------------------------------------

beast <- read.beast('FIV_DG_cauchy_MCC.tree')

get_taxa_name(tree_view = NULL, node = NULL)

#if there are problems in the tip names - can import new ones in a csv (old names first, new names second)
nn <- read.csv('newNames.csv', head=F)
beast <- taxa_rename(beast, nn)

#plot tree to see node position
p <- ggtree(beast) + geom_tiplab(size=2)+ geom_text2(aes(subset=!isTip, label=node), hjust=-.3) 
p

#add branch lengths etc
p1 <- ggtree(beast , mrsd="2015-01-01") + theme_tree2()+ 
  geom_tiplab(align=TRUE, linetype='dashed', linesize=.3, size=2) + 
  #geom_range("length_0.95_HPD", color='red', size=2, alpha=.5) + #gets length estimates 
  geom_text2(aes(label=round(as.numeric(posterior), 2), 
                 subset=as.numeric(posterior)> 0.9, 
                 x=branch), vjust=0)
p1
#can isolate particular sections of the tree
viewClade(p+geom_tiplab(size=2), node=218)

WStree <- tree_subset(beast, node='X1399_WS_2011-08-03', levels_back = 8)

as.phylo(WStree)
head(WStree)


write.beast(beast, file = "WSFR.tre", translate = TRUE,
            tree.name = "WSFR")

#have to reopen it to make it a phylo object - bit of bug at the moment
b <- read.nexus('WSFR.tre')
str(b)
plot(b)


#------------------------------------------------------------------------
##################Classify sequences into each contrast##################
#------------------------------------------------------------------------

comm <- as.data.frame(get_taxa_name(tree_view = NULL, node = NULL))
                      
names(comm) <- c('date') 
comm <- comm %>% separate(date,c("ID","Site","Year", "month", "day"))

#for FFV
noHunt <-comm %>% filter(Year < 2010)
Hunt <-comm %>% filter(Year > 2010)

noHunt <- add_column(noHunt, noHunt = 1, Hunt=0)
Hunt <- add_column(Hunt, noHunt = 0, Hunt=1)
Combined <- bind_rows(noHunt, Hunt,id = NULL)
commDate <- Combined %>% unite(date,Year, month, day, sep='-')
comm <- commDate %>% unite(seq,ID, Site, date )
rownames(comm) <- comm$seq
comm <- select(comm,-c(seq))
comm <-t(as.matrix(comm))

#forFIV
noHuntWS <-comm %>% filter(Year <= 2010, Site=='WS')
HuntWS <-comm %>% filter(Year > 2010, Site=='WS')

noHuntFR <-comm %>% filter(Year <= 2010, Site=='FR')
HuntFR <-comm %>% filter(Year > 2010, Site=='FR')

noHuntFR <- add_column(noHuntFR, noHuntFR = 1, HuntFR=0, noHuntWS=0,HuntWS=0 )
HuntFR <- add_column(HuntFR, noHuntFR = 0, HuntFR=1, noHuntWS=0,HuntWS=0 )
HuntWS <- add_column(HuntWS, noHuntFR = 0, HuntFR=0, noHuntWS=0,HuntWS=1 )
noHuntWS <- add_column(noHuntWS, noHuntFR = 0, HuntFR=0, noHuntWS=1,HuntWS=0 )

Combined <- bind_rows(noHuntWS, HuntWS,noHuntFR, HuntFR, id = NULL)
commDate <- Combined %>% unite(date,Year, month, day, sep='-')
comm <- commDate %>% unite(seq,ID, Site, date )
rownames(comm) <- comm$seq
comm <- select(comm,-c(seq))
comm <-t(as.matrix(comm))
#------------------------------------------------------------------------
##Extraction of the spatio-temporal information contained in 1,000 trees
#------------------------------------------------------------------------

library(diagram)
library(vioplot)
library(lubridate)
library(OutbreakTools)

allTrees = scan(file="FFVenvPol_relaxedC _GMRF_3.trees.txt", what="", sep="\n", quiet=TRUE)
burnIn = 10001; index1 = which(grepl("tree STATE_0 ",allTrees)) + burnIn
samplingFrequency = (length(allTrees)-index1)/1000
allTrees1 = allTrees[1:(which(grepl("tree STATE_0 ",allTrees))-1)]
allTrees2 = c(allTrees[seq(((index1-1)+samplingFrequency),length(allTrees),samplingFrequency)],"End;")
#can write this as an object which could be useful
write(c(allTrees1,allTrees2), "FFVenvPol_relaxedC _GMRF_3_1000.trees")

beast1000 <- read.nexus('FIV_DG_cauchy_1000.trees')
get_taxa_name(tree_view = NULL, node = NULL)
beast1000 <- taxa_rename(beast1000, nn)

#plot one tree to see if it works
plot(beast1000[[1]])

#####change names on a multiphylo object
nn <- read.csv('newNames.csv', head=F)
beast <- taxa_rename(beast, nn)

beast1000renamed<-lapply(beast1000,taxa_rename, name=nn)
plot(beast1000renamed[[1]])
#------------------------------------------------------------------------
##################MultiPhyloPD##################
#------------------------------------------------------------------------

MultiPhyloPD <- function(ph, comm,...){
  out <- list(NA, ncol = nrow(comm), nrow = length(ph)) 

  for(i in 1:length(ph)){
    prunedphy <- prune.sample(comm,ph[[i]])
    pd <- pd.query(prunedphy, comm, standardize = TRUE, 
                   null.model="uniform", reps=1000)
    
    out[[i]] <- pd
    
  }
  out <- t(as.data.frame(out))
  colnames(out) <- row.names(comm)
  row.names(out) <- c(1:length(ph))
  as.data.frame(out)
  

}

p <- MultiPhyloPD(beast1000renamed, comm)
str(p)

#make a tidy object
df <- gather(p, key, value) 
df$key <-as.factor(df$key)

#plot
plot <- ggplot(df, aes(y=value,x=key))
plot + geom_boxplot()+labs(title = 'FiV WSFR', y="SES.PD", x = 'Contast')+ theme_bw()


#------------------------------------------------------------------------
##################Metrics sepparately##################
#------------------------------------------------------------------------

faithsSES <- ses.pd(comm, b, null.model = c("taxa.labels"),
       runs = 9999, iterations = 10000)
faithsSES

#phylomeasure
library(PhyloMeasures)

pm <- pd.query(b, comm, standardize = TRUE, 
           null.model="uniform", reps=1000)
pm <-as.data.frame(pm)
row.names(pm) <-row.names(comm)
#pez results harder to interpret
pddata <- comparative.comm(b, comm)
faiths <- .pd(pddata, include.root = TRUE, abundance.weighted = FALSE)
faiths.ses <- generic.null(pddata, c(.pd), null.model=c("taxa.labels"), permute = 10000)
#row.names(faiths.ses) <-row.names(comm)
as.data.frame(faiths.ses)


