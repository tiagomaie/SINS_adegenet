#####################################################
#####################################################
##Testing parts of SINS_to_adegenet in the real simulation data
#####################################################
#####################################################

library("adegenet")
library("ade4")
library(hierfstat)
library(pegas)
library(ggplot2) #to plot
library(reshape2) #to melt/reshape tables
library(gridExtra) #to put several plots into 1 image, easily
library(cowplot) #publication ready plots?

##dir with original SINS_to_adegenet.r file
setwd("/home/tiago/Documents/Documents/R/R_scripts/SINS_adegenet/")
source(file = "SINS_to_adegenet.r")

##mount folder from elephant
##sshfs tmaie@elephant01a:/Users/tmaie/SINS_Sampler ~/server_folder

##dir with output from SINS_Sampler (mounted from elephant01a)
#setwd("/home/tiago/server_folder/dist/output/server_sim_Test_Clusters/server_sim/Adegenet_sim1/")
setwd("/home/tiago/server_folder/dist/output/server_sim_Test_Clusters_discrete_deme/server_sim_discrete_deme/Adegenet_sim1/")
getwd()

the_generation_list = Set_generations_list(10,6000,5)
the_layers_list = Set_layers_list("layer0")
the_markers_list = Set_markers_list("A1","A2","A3","A4","A5","A6","A7","A8","A9","A10")

the_file_names = Get_file_names(the_generation_list, the_layers_list)
raw_data = Get_data_from_files(the_file_names = the_file_names,
                               set_markers = the_markers_list,
                               generations_list = the_generation_list)



test.pop = raw_data[[50]][pop=c("pop6A","pop6B")]
Hs(raw_data[[50]][pop=c("pop6A","pop6B","pop6C","pop6D","pop6E")])
Hs(raw_data[[50]])
compute.Hs.df.colnames = c("pop1A","pop2A","pop3A","pop4A","pop4B","pop5A","pop5B","pop5C","pop6A","pop6B","pop6C","pop6D","pop6E")
compute.Hs.df = as.data.frame(matrix(data = NA,nrow = length(raw_data),ncol = length(compute.Hs.df.colnames)))
colnames(compute.Hs.df) = compute.Hs.df.colnames
head(compute.Hs.df["pop3A"])

##COMPUTE EXPECTED HETEROZYGOSITY PER DEME
compute.Hs.raw.data = NULL
compute.Hs.df = data.frame()


## iterate through the raw_data,
## if there are more than 1 allele in either of the locus
## compute Hs() and stack its values on top of each other
## bind the generation corresponding to each computed value
for(i in 1:length(raw_data)){
  if(length(raw_data[[i]]@all.names$loc1)>1 || length(raw_data[[i]]@all.names$loc2)>1){
    test.var = stack(Hs(raw_data[[i]]))
    test.var.cbind = cbind(test.var,gen = raw_data[[i]]@other$generation)
    compute.Hs.df = rbind(compute.Hs.df,test.var.cbind)
  }
}

Compute.Hs.make.df = function (X){
  compute.Hs.df = data.frame()
  for(i in 1:length(X)){
    if(length(X[[i]]@all.names$loc1)>1 || length(X[[i]]@all.names$loc2)>1){
      test.var = stack(Hs(X[[i]]))
      test.var.cbind = cbind(test.var,gen = X[[i]]@other$generation)
      compute.Hs.df = rbind(compute.Hs.df,test.var.cbind)
    }
  }
  compute.Hs.df.numeric.gen = compute.Hs.df
  ## transform gen col from factor to numeric so that we are able to use gen as a continuous scale for the plot
  compute.Hs.df.numeric.gen$gen = as.numeric(as.character(compute.Hs.df.numeric.gen$gen))
  return(compute.Hs.df.numeric.gen)
}

tail(compute.Hs.df)

compute.Hs.df.numeric.gen = compute.Hs.df
## transform gen col from factor to numeric so that we are able to use gen as a continuous scale for the plot
compute.Hs.df.numeric.gen$gen = as.numeric(as.character(compute.Hs.df.numeric.gen$gen))

ggplot(data = compute.Hs.df.numeric.gen)+
  geom_line(aes(x=gen, y=values,group=ind,color=ind), ## group by ind, color by ind
            subset(compute.Hs.df.numeric.gen,
                   ind=="pop1A" | ind=="pop2A" | ind=="pop3A" | 
                     ind=="pop4A" | ind=="pop4B" | ind=="pop5A" | 
                     ind=="pop5B" | ind=="pop5C" | ind=="pop6A" | 
                     ind=="pop6B" | ind=="pop6C" | ind=="pop6D" | 
                     ind=="pop6E"), alpha=1.0)+ ## subset the data, choose which ind(pops) we want to plot
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  ggtitle("Expected Heterozigosity for each deme")+
  xlab("Generations")+
  ylab("Expected Heterozigosity")+
  theme(legend.position="right")+
  scale_x_continuous(breaks=c(0,3000,6000))+
  #scale_x_discrete(labels=compute.Hs.df$gen,breaks=compute.Hs.df$gen, limits = compute.Hs.df$gen)+
  geom_vline(xintercept = 3000,linetype = "longdash")+
  coord_cartesian(ylim=c(0,1), xlim=c(0,6000))

## Plot Grouped demes (pop1,pop2,pop3,pop4,pop5,pop6)
Transform_into_simple_pop = function(genind_array_obj){
  # Transforms a genind array grouping populations into simpler(read broader) groups
  # pop1A will be pop1
  # its defined for the following population code pop[1-6][A-E]
  # TODO refactor so that it uses regex instead of having it explicitly defined
  
  Group.demes.function = function(X){
    ## Group demes according to their population code (eg "pop1A" will be "pop1", "pop6E" will be "pop6")
    varN = NULL
    for(i in popNames(X)){
      if(i=="pop1A"|| i == "pop1B" || i == "pop1C" || i == "pop1D" || i == "pop1E"){
        varN = c(varN, "pop1")
      } else if (i == "pop2A"|| i == "pop2B" || i == "pop2C" || i == "pop2D" || i == "pop2E"){
        varN = c(varN, "pop2")
      } else if (i == "pop3A"|| i == "pop3B" || i == "pop3C" || i == "pop3D" || i == "pop3E"){
        varN = c(varN, "pop3")
      } else if (i == "pop4A" || i == "pop4B" || i == "pop4C" || i == "pop4D" || i == "pop4E"){
        varN = c(varN, "pop4")
      } else if (i == "pop5A" || i == "pop5B" || i == "pop5C"|| i == "pop5D" || i == "pop5E"){
        varN = c(varN, "pop5")
      } else if (i == "pop6A" || i == "pop6B" || i == "pop6C" || i == "pop6D" || i == "pop6E"){
        varN = c(varN, "pop6")
      }
    }
    popNames(X) = varN
    return(X)
  }
  
  raw_data.grouped.demes = genind_array_obj
  
  raw_data.grouped.demes = lapply(X = raw_data.grouped.demes,FUN = Group.demes.function)
  
  return(raw_data.grouped.demes)
}

hs.df.grouped.demes = Compute.Hs.make.df(raw_data.grouped.demes)
tail(hs.df.grouped.demes)


ggplot(data = hs.df.grouped.demes)+
  geom_line(aes(x=gen, y=values,group=ind,color=ind), ## group by ind, color by ind
            subset(hs.df.grouped.demes,
                   ind=="pop1" | ind=="pop2" | ind=="pop3" | 
                     ind=="pop4" | ind=="pop5" | ind=="pop6"), alpha=1.0)+ ## subset the data, choose which ind(pops) we want to plot
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  ggtitle("Expected Heterozigosity for group of demes")+
  xlab("Generations")+
  ylab("Expected Heterozigosity")+
  theme(legend.position="right")+
  scale_x_continuous(breaks=c(0,3000,6000))+
  #scale_x_discrete(labels=compute.Hs.df$gen,breaks=compute.Hs.df$gen, limits = compute.Hs.df$gen)+
  geom_vline(xintercept = 3000,linetype = "longdash")+
  coord_cartesian(ylim=c(0,1), xlim=c(0,6000))



###################
###################
###################

Get_genind_from_gen = function(a_genind_array, generation){
  ## simple function to get a genind object from an array, from a given specified generation
  # Args:
  #     a_genind_array - an array of genind objects
  #     generation - a number with the generation that we want to search for
  #
  # Returns:
  #     The genind object with the given generation if it finds it or a message if it doesn't.
  
  for(i in 1:length(a_genind_array)){
    current.gen = as.numeric(a_genind_array[[i]]@other$generation)
    
    if(current.gen == generation){
      return(a_genind_array[[i]])
    } 
  }
  return(paste("Generation '",generation, "' does not exist in the given genind array",sep=""))
}



###################
####Fst#####
###################
Compute.Fst.plot = function(gen.ind.object,matrix.as.triangle = T, upper.triangle = F, has.in.graph.txt = F){
  compute.Fst = gen.ind.object
  compute.Fst.generation = gen.ind.object@other$generation
  #compute.Fst = compute.Fst[pop=c("pop1A","pop2A","pop3A","pop4A","pop4B","pop5A","pop5B","pop5C","pop6A","pop6B","pop6C","pop6D","pop6E")]
  popNames(compute.Fst)
  compute.Fst = pairwise.fst(compute.Fst,res.type = "matrix")
  
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  
  get.matrix.triangle=matrix.as.triangle
  get.upper.triangle=upper.triangle
  if(get.matrix.triangle){
    if(get.upper.triangle){
      # Skip this to get the whole matrix instead of just one of the triangles  
      compute.Fst.upper.tri=get_upper_tri(compute.Fst)
      melted.compute.Fst = melt(compute.Fst.upper.tri,id.vars = names(compute.Fst),na.rm = T,factorsAsStrings = T)
    } else {
      compute.Fst.lower.tri=get_lower_tri(compute.Fst)
      melted.compute.Fst = melt(compute.Fst.lower.tri,id.vars = names(compute.Fst),na.rm = T,factorsAsStrings = T)
    }
  } else {
    melted.compute.Fst = melt(compute.Fst,id.vars = names(compute.Fst),na.rm = T,factorsAsStrings = T)
  }
  
  ## Check min and max values to adjust color palette limits
  round(max(melted.compute.Fst$value)+0.05,digits = 1)
  min(melted.compute.Fst$value)
  ## round values in matrix
  melted.compute.Fst$value = round(melted.compute.Fst$value,3)
  
  fst.heatmap = ggplot(data = melted.compute.Fst, aes(Var1, Var2, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "yellow", high = "red", mid = "orange", 
                         midpoint = 0.15, limit = c(0.0,0.30001), space = "Lab", 
                         name="Pairwise\nFst values") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    geom_rect(inherit.aes = F,aes(xmin=1-0.4, xmax=1+0.4, ymin=1-0.4, ymax=1+0.4),color="#999999", linetype=3,fill=NA)+
    geom_rect(inherit.aes = F,aes(xmin=2-0.4, xmax=2+0.4, ymin=2-0.4, ymax=2+0.4),color="#999999", linetype=3,fill=NA)+
    geom_rect(inherit.aes = F,aes(xmin=3-0.4, xmax=3+0.4, ymin=3-0.4, ymax=3+0.4),color="#999999", linetype=3,fill=NA)+
    geom_rect(inherit.aes = F,aes(xmin=4-0.4, xmax=5+0.4, ymin=4-0.4, ymax=5+0.4),color="#999999", linetype=3,fill=NA)+
    geom_rect(inherit.aes = F,aes(xmin=6-0.4, xmax=8+0.4, ymin=6-0.4, ymax=8+0.4), color="#999999", linetype=3,fill=NA)+
    geom_rect(inherit.aes = F,aes(xmin=9-0.4, xmax=13+0.4, ymin=9-0.4, ymax=13+0.4),color="#999999", linetype=3,fill=NA)+
    coord_fixed()
  fst.heatmap
  
  fst.heatmap= fst.heatmap+
    #geom_text(aes( Var1,Var2, label = value), color = "black", size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      #panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.4, 0.7),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))+
    ggtitle(paste("Pairwise Fst - Generation ",compute.Fst.generation))
  
  if(has.in.graph.txt){
    fst.heatmap = fst.heatmap + geom_text(aes( Var1,Var2, label = value), color = "black", size = 4)
  }
  
  return(fst.heatmap)
}
#####
#####
#####

Compute.Fst.plot(raw_data[[598]],has.in.graph.txt = T)## just before environmental change
Compute.Fst.plot(raw_data[[700]],has.in.graph.txt = T)## just after environmental change

a = Compute.Fst.plot(raw_data[[1199]],has.in.graph.txt = T)
b = Compute.Fst.plot(raw_data[[1199]],has.in.graph.txt = F)
ggplotly(a)
ggplotly(b)

the.fst.plot = Compute.Fst.plot(raw_data[[50]])
raw_data[[99]]@other$generation
raw_data[[119]]@other$generation
raw_data[[11]]@pop
length(raw_data)

########################################################################
########################################################################
############  MAKE GRID OUT OF FST PLOTS
########################################################################
########################################################################

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

the_legend = get_legend(Compute.Fst.plot(Get_genind_from_gen(raw_data,100))+ theme(legend.position=c(0.5,0.5),legend.justification="center"))
fst.30 = Compute.Fst.plot(Get_genind_from_gen(raw_data,30))+ theme(legend.position="none",axis.text.x=element_text(vjust=1.5))
fst.50 = Compute.Fst.plot(Get_genind_from_gen(raw_data,50))+ theme(legend.position="none",axis.text.x=element_text(vjust=1.5))
fst.75 = Compute.Fst.plot(Get_genind_from_gen(raw_data,75))+ theme(legend.position="none",axis.text.x=element_text(vjust=1.5))
fst.100 = Compute.Fst.plot(Get_genind_from_gen(raw_data,100))+ theme(legend.position="none",axis.text.x=element_text(vjust=1.5))
fst.1000 = Compute.Fst.plot(Get_genind_from_gen(raw_data,1000))+ theme(legend.position="none",axis.text.x=element_text(vjust=1.5))
fst.3000 = Compute.Fst.plot(Get_genind_from_gen(raw_data,3000))+ theme(legend.position="none",axis.text.x=element_text(vjust=1.5))
fst.4000 = Compute.Fst.plot(Get_genind_from_gen(raw_data,4000))+ theme(legend.position="none",axis.text.x=element_text(vjust=1.5))
fst.5000 = Compute.Fst.plot(Get_genind_from_gen(raw_data,5000))+ theme(legend.position="none",axis.text.x=element_text(vjust=1.5))
fst.6000 = Compute.Fst.plot(Get_genind_from_gen(raw_data,6000))+ theme(legend.position="none",axis.text.x=element_text(vjust=1.5))
blank_plot = ggplot()+geom_blank(aes(1,1))+cowplot::theme_nothing()

fst.30$layers[[8]] = NULL
fst.50$layers[[8]] = NULL
fst.75$layers[[8]] = NULL
fst.100$layers[[8]] = NULL
fst.1000$layers[[8]] = NULL
fst.3000$layers[[8]] = NULL
fst.4000$layers[[8]] = NULL
fst.5000$layers[[8]] = NULL
fst.6000$layers[[8]] = NULL

a = grid.arrange(
  blank_plot,
  the_legend,
  blank_plot,
  fst.30,
  fst.50,
  fst.75,
  fst.100,
  fst.1000,
  fst.3000,
  fst.4000,
  fst.5000,
  fst.6000,
  ncol=3, nrow=4,
  widths = c(5, 5, 5), heights = c(0.7,3, 3, 3)
  
)
####################################################################
####################################################################
####################################################################


####################################################################
## Use code below to get all the images needed for a gif of Fst plots
## (potentially change the format to a more lossless format instead of jpeg)
## Save plots to list
plot_list = list()
initial.genind.obj = 99
for(i in seq(from = initial.genind.obj,to = length(raw_data),by = 20)){
  fst.plot.gen = raw_data[[i]]@other$generation
  path.to.folder = file.path("~","Documents","R_Plots","Fst_plots",paste(fst.plot.gen,"gen_fst_plot",".jpg",sep=""))
  
  jpeg(filename = path.to.folder)
  print(Compute.Fst.plot(raw_data[[i]]))
  #plot_list[[i-initial.genind.obj+1]] = Compute.Fst.plot(raw_data[[i]])
  dev.off()
  
}
####################################################################



#########################################################
#########################################################
###################ISOLATION BY DISTANCE###################
#########################################################
#########################################################

Make_pairFst_distance_plots = function(genind_array, generation){
  #genind_array=raw_data
  #generation=6000
  
  # we create two instances of data because the pairwise.fst function takes genind objects instead of genepop
  raw_data_ind = Get_genind_from_gen(a_genind_array = genind_array,generation = generation)
  raw_data_pop = genind2genpop(raw_data_ind,process.other = T)
  
  
  
  distance.geographic = dist(raw_data_pop@other$xy)
  
  ## make dist obj into dataframe (so that we can plot it)
  Dist.to.df <- function(inDist) {
    if (class(inDist) != "dist") stop("wrong input type")
    A <- attr(inDist, "Size")
    B <- if (is.null(attr(inDist, "Labels"))) sequence(A) else attr(inDist, "Labels")
    if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
    if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
    data.frame(
      row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
      col = rep(B[-length(B)], (length(B)-1):1),
      value = as.vector(inDist))
  }
  
  
  compute.Fst.ibd = pairwise.fst(raw_data_ind,res.type = "dist")
  #is.euclid(compute.Fst.ibd)
  #cailliez(compute.Fst.ibd)
  compute.Fst.ibd = Dist.to.df(compute.Fst.ibd)
  distance.geographic.ibd = Dist.to.df(distance.geographic)
  
  dist.df = data.frame(distance.geographic.ibd,compute.Fst.ibd)
  
  the_plot = ggplot(dist.df)+
    geom_point(aes(x = dist.df$value,y = dist.df$value.1))+
    geom_smooth(aes(x = value, y = value.1), method = 'lm', se = F, fullrange=T)+
    scale_x_continuous(expand=c(0.0,0.0), limits=c(0,12),breaks = seq(0,12,1)) +
    scale_y_continuous(expand=c(0.0,0.0), limits=c(0,0.4)) +
    coord_cartesian(xlim = c(0, 12), ylim = c(0, 0.4))+
    ylab("Paiwise Fst")+
    xlab("Geographical distance between demes")+
    ggtitle(paste("Generation ",raw_data_ind@other$generation))+
    cowplot::theme_cowplot()
  
  return(the_plot)
}

ggplotly()

dist.plot.gen3000 = Make_pairFst_distance_plots(genind_array = raw_data, generation = 3000)
dist.plot.gen3050 = Make_pairFst_distance_plots(genind_array = raw_data, generation = 3050)
dist.plot.gen3200 = Make_pairFst_distance_plots(genind_array = raw_data, generation = 3200)
dist.plot.gen3600 = Make_pairFst_distance_plots(genind_array = raw_data, generation = 3600)
dist.plot.gen4000 = Make_pairFst_distance_plots(genind_array = raw_data, generation = 4000)
dist.plot.gen6000 = Make_pairFst_distance_plots(genind_array = raw_data, generation = 6000)
grid.arrange(
  dist.plot.gen3000,
  dist.plot.gen3050,
  dist.plot.gen3200,
  dist.plot.gen3600,
  dist.plot.gen4000,
  dist.plot.gen6000,
  ncol=3,nrow=2)


#####TODO

raw_data_ind = Get_genind_from_gen(a_genind_array = raw_data,generation = 3000)
raw_data_pop = genind2genpop(Get_genind_from_gen(a_genind_array = raw_data,generation = 6000),process.other = T)

distance.genetic = dist.genpop(raw_data_pop, method = 3)
distance.geographic = dist(raw_data_pop@other$xy)

Dgeo = distance.geographic
Dgen = distance.genetic

isolation.by.distance = mantel.randtest(distance.genetic,distance.geographic,nrepet = 1000)

plot(isolation.by.distance)
plot(distance.geographic,distance.genetic)
abline(lm(distance.genetic~distance.geographic), col="red",lty=2)


library(MASS)
dens <- kde2d(Dgeo,Dgen, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo, Dgen, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(Dgen~Dgeo))
title("Isolation by distance plot")


connection.network = chooseCN(raw_data_pop@other$xy,ask=F,type = 2)
D = dist(raw_data_pop@tab)
mon1=monmonier(xy = raw_data_pop@other$xy,D,connection.network)

pco1 <- dudi.pco(D,scannf=FALSE,nf=1)
barplot(pco1$eig,main="Eigenvalues")
D=dist(pco1$li)

mon1=monmonier(xy = raw_data_pop@other$xy,D,connection.network)
coords.monmonier(mon1)
plot(mon1)


###################
###################
###################


#raw_data[[1199]]
a_summary = summary(raw_data[[1199]])$Hobs
# t() is the transpose function
data_to_plot = data.frame(t(a_summary),row.names = raw_data[[1199]]@other$generation) 

a=Do_summary_calcs(raw_data[[1199]])
a
length(levels(raw1199@pop))

###loop through gens

plot.raw.fst.loc1 = c()
plot.raw.fst.loc2 = c()
plot.raw.He = c()
options(error=recover)
options(error=NULL)

for(i in 1:length(raw_data)){
  
  test.genind2df.WCfst=gtrunchier
  head(test.genind2df.WCfst[,-2])
  test.genind2df.WCfst = genind2df(raw1199 ,sep = " ")
  allele.count(raw1199)
  allelic.richness(raw1199)
  pairwise.WCfst(dat=test.genind2df.WCfst, diploid = T)
  wc(ndat=test.genind2df.WCfst, diploid = T)
  length(dim(head(genind2df(raw1199,sep = ""))))
  
  plot.raw.fst.loc1 = c(plot.raw.fst.loc1, Fst(as.loci(raw_data[[i]]))[,"Fst"]["loc1"])
  plot.raw.fst.loc2 = c(plot.raw.fst.loc2, Fst(as.loci(raw_data[[i]]))[,"Fst"]["loc2"])
  if(length(levels(raw_data[[i]]@pop))>1){
    #plot.raw.He = c(plot.raw.He,Hs(raw_data[[i]]))
  }
  #fstat(raw_data[[5]])
  #Fst(as.loci(raw_data[[5]]))[,"Fst"]["loc2"]
  #Hs(raw_data[[5]])
}
par(mfrow=c(1,1))

raw.df = data.frame(gen = the_generation_list, fst=plot.raw.fst.loc1)

head(raw.df, n=100)

raw.ggplot = ggplot(data = raw.df)+
  geom_line(aes(x = gen, y = fst))+
  geom_vline(xintercept = 3000,linetype = "longdash")
raw.ggplot

###



head(raw_data[[1199]]@other$xy)
raw_data[[1199]][pop=6]@other$xy[,]
raw_data[[1199]][258:260,]@other$xy
raw_data[[1199]][258:260,]@tab
fstat(raw_data[[1199]][pop=1:6]) #Fst values less than 0.05 are usually considered negligible

the_pop_six=raw_data[[1199]][pop=6]
gen2ind_the_pop_six = genind2genpop(raw_data[[1199]],process.other = TRUE)
gen2ind_the_pop_six@tab
gen2ind_the_pop_six@other$xy
raw1199 = raw_data[[1199]]
setPop(raw1199)=as.factor("pop1")
raw1199 = pop(as.factor("pop1"))
ind1_raw1199 = raw1199[1]
pop(ind1_raw1199)=as.factor("pop123")
raw1199[1]@pop=as.factor("pop3")

###repooling pops
raw1199 -> raw1199.redefPops

popNames(raw1199)
popNames(raw1199.redefPops)
popNames(raw1199.redefPops) = c("pop1","pop2","pop3","pop4","pop4","pop5","pop5","pop5","pop6","pop6","pop6","pop6","pop6") 

raw1 = raw_data[[1]]
popNames(raw1)
popNames(raw1) =  c("pop1","pop2","pop3","pop4","pop4","pop5","pop5","pop5","pop6","pop6","pop6","pop6","pop6") 


popNames(raw1199.redefPops)=rep(as.factor("pop1"),times=2)
###

fstat_raw1199 = fstat(raw1199)
fst_raw1199 = fstat_raw1199[1]
fit_raw1199 = fstat_raw1199[3]
fis_raw1199 = fstat_raw1199[4]
hw.test(raw1199,B = 500)
?hw.test
?t.test
Hs(raw1199) #expected heterozygosity within populations
sumcalc_raw1199 = Do_summary_calcs(raw1199)
Fst(as.loci(raw1199))

#if p-value is less than 0.01 accept alternative hypothesis
#(default Ha is that our obs data is greater than the simulated)
Gtest_1199 <- gstat.randtest(raw1199,nsim=200)



raw1199


?gstat.randtest
Gtest_1199$obs

plot(Gtest_1199)
sumcalc_raw1199[[3]]$estimate

seppop(raw1199)$pop1A
raw1199.inbreeding= inbreeding(raw1199)
raw1199.Fbar = sapply(raw1199.inbreeding, mean)
hist(raw1199.Fbar)

raw1199.pairwise.fst = pairwise.fst(raw1199)
is.euclid(raw1199.pairwise.fst)
raw1199.lingoes = lingoes(raw1199.pairwise.fst, TRUE)
plot(raw1199.lingoes)

pop(raw1199) =  as.factor(c(rep(x = "pop2",times=nrow(raw1199@tab)-1),"pop1"))
indlist=NULL
for(i in 1:nrow(raw1199@other$xy)){
  if(identical(x = raw1199@other$xy[i,],y=as.integer(c(0,0)))){
    someind = raw1199[i,]
    pop(someind)=as.factor("00")
    
    #raw1199@other$xy[i,]
    
  }else if(identical(raw1199@other$xy[i,] , as.integer(c(0,2)))){
    someind = raw1199[i,]
    pop(someind)=as.factor("02")
    
    
  }else if(identical(raw1199@other$xy[i,] , as.integer(c(2,0)))){
    someind = raw1199[i,]
    pop(someind)=as.factor("20")
    
    
  }else if(identical(raw1199@other$xy[i,] , as.integer(c(4,4)))){
    someind = raw1199[i,]
    pop(someind)=as.factor("44")
    
    
  }else if(identical(raw1199@other$xy[i,] , as.integer(c(0,4)))){
    someind = raw1199[i,]
    pop(someind)=as.factor("04")
    
    
  }else if(identical(raw1199@other$xy[i,] , as.integer(c(4,0)))){
    someind = raw1199[i,]
    pop(someind)=as.factor("40")
    
    
  }else if(identical(raw1199@other$xy[i,] , as.integer(c(8,0)))){
    someind = raw1199[i,]
    pop(someind)=as.factor("80")
    
    
  }else if(identical(raw1199@other$xy[i,] , as.integer(c(0,8)))){
    someind = raw1199[i,]
    pop(someind)=as.factor("08")
    
    
  }else if(identical(raw1199@other$xy[i,] , as.integer(c(8,8)))){
    someind = raw1199[i,]
    pop(someind)=as.factor("88")
    
    
  }
  indlist = c(indlist,someind)
  
  
}
newgenind=repool(indlist)
newgenind@pop

tail(the_pop_six@other$xy[40:41,])
tail(the_pop_six@tab)
the_pop_six@tab[40:41,]
names(the_pop_six@tab[50,])

length(raw_data[[1199]]@other$xy)
raw_data[[1199]]@other$xy[,1]
raw_data[[1199]]@other[zxcasdouahsodhaoduhawoudaowuhdaowudhoawdu=2]
raw_data[[1199]][xy]

asdasd = raw_data$`6000`@other$xy[250,]



thepop1=NULL
popNames(raw_data[[1000]])
thepop1 = raw_data[[1000]][pop=6]

data_to_plot = data.frame()

pop1listHe = data.frame()

pop1sum=NULL
data_sum_list = list()
generation_list = NULL
for(i_data in seq(1:length(raw_data))){
  
  if(length(levels(raw_data[[i_data]]@pop))>1){ ##if there is more than ONE pop in that generation
    generation_list = c(generation_list,raw_data[[i_data]]@other$generation)
    data_sum = summary(raw_data[[i_data]])
    data_sum_list = rbind(data_sum_list,data_sum$Hexp)
  }
  
}
head(data_sum_list)
df_data_sum_list=data.frame(loc1 = as.numeric(data_sum_list[,1]),
                            loc2 = as.numeric(data_sum_list[,2]),
                            gen = as.numeric(generation_list))

#rownames(df_data_sum_list) =generation_list
#colnames(df_data_sum_list) = the_markers_list

head(df_data_sum_list)

str(as.numeric(data_sum_list[,1]))


rsh_df_data_sum_list = melt(data = df_data_sum_list, id.vars = c("gen")) #reshapes the dataframe 
head(rsh_df_data_sum_list)

the_plot = ggplot(data = df_data_sum_list) + 
  geom_line(aes(x = as.numeric(rownames(df_data_sum_list)), y=df_data_sum_list$A1, color="A1"))+
  geom_line(aes(x = as.numeric(rownames(df_data_sum_list)), y=df_data_sum_list$A2, color="A2"))+
  geom_vline(xintercept = 3000,linetype = "longdash")+
  xlab("Generations")+ylab("Hexp")+
  coord_cartesian(ylim=c(0,1), xlim=c(0,6000))
the_plot

ggplot(data = rsh_df_data_sum_list,
       aes(x= gen, y = value, colour=variable))+
  geom_line()+
  geom_vline(xintercept = 3000,linetype = "longdash")+
  xlab("Generations")+ylab("Hexp")+
  coord_cartesian(ylim=c(0,1), xlim=c(0,6000))
head(rsh_df_data_sum_list$gen)




plot_df = ggplot(data=values_All_df, aes(x=all_generations,y=He))+
  geom_jitter(aes(color=dataset))+
  geom_smooth(color="black")+
  #facet_grid(.~dataset)+
  facet_grid(dataset~.)+
  geom_vline(xintercept = 3000,linetype = "longdash")+
  #geom_text(aes(x=3000+700,y=0+0.1,label="fragmentation\nevent",family="sans"),size=5)+
  #theme_light()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  ggtitle("Expected Heterozigosity for each demographic cluster")+
  xlab("Generations")+
  ylab("Expected Heterozigosity")+
  theme(legend.position="none")+
  scale_x_continuous(breaks=c(0,3000,6000))+
  coord_cartesian(ylim=c(0,1), xlim=c(0,6000))


#######################################
#######################################
#######################################
#SPCA TUTORIAL
#######################################
#######################################
#######################################

grp <- find.clusters(raw1199, max.n.clust=40)
names(grp)
raw1199.spca.grp = table(pop(raw1199),grp$grp)
nrow(raw1199.spca.grp)
table.value(raw1199.spca.grp, col.lab=paste("inf", 1:ncol(raw1199.spca.grp)),
            row.lab=paste("ori", 1:nrow(raw1199.spca.grp)))

raw1199.dapc <- dapc(raw1199, grp$grp)

scatter(raw1199.dapc)

compoplot(raw1199.dapc, posi="bottomright",
          txt.leg=paste("Cluster", 1:6), lab="",
          ncol=1, xlab="individuals", col=funky(6))

adegenetServer(what = "DAPC")

#######################################
#######################################
#######################################
library(plotly)
