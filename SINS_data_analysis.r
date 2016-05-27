library(adegenet)
library(ade4)
library(hierfstat)
library(pegas)
library(ggplot2) #to plot
library(reshape2) #to melt/reshape tables
library(gridExtra) #to put several plots into 1 image, easily
library(cowplot) #publication ready plots?
library(plotly) # interactive plots

###########################################################################
###########################################################################
#++++SECTION 1 - FUNCTIONS TO GET DATA
###########################################################################
###########################################################################

################
#SET MARKERS
################
#+++++++++++++++
Set_markers_list = function(...){
  # Make list of the given markers (e.g. "A1","A2","A3", ...)
  the_markers_list=NULL
  the_markers_list = c(...)
  return(the_markers_list)
}

################
#SET GENERATIONS
################
#+++++++++++++++
Set_generations_list = function(initialTimeStep, finalTimeStep, timeStepIncrement){
  # Make sequence of all the generations that we want to study (e.g. from gen 1 to 100 by increments of 5)
  the_generations_list=NULL
  the_generations_list = seq(initialTimeStep,finalTimeStep,timeStepIncrement)
  return(the_generations_list)
}

################
#SET LAYERS
################
#+++++++++++++++
Set_layers_list = function(...){
  # Make list of the given layers (e.g. "layer0", "layer1")
  the_layers_list=NULL
  the_layers_list = c(...)
  return(the_layers_list)
}

################
#GET FILE NAMES
################
#+++++++++++++++
Get_file_names = function(generations_list, layers_list){
  # Build the names of the files that we want to analyse
  the_file_names = NULL
  
  for(j in generations_list){
    
    the_file_names = rbind(the_file_names, paste(j,"_",layers_list,"_",sep = ""))
    
  }
  return(the_file_names)
}

################
#READ DATA
################
#+++++++++++++++
Read_data_func = function(name_of_file, markers_list){
  
  current_generation = NULL
  current_generation = strsplit(x=name_of_file,split = "_")[[1]][1]
  
  data_file = NULL
  marker_data = NULL
  bind_markers = NULL
  a_file_name = NULL
  
  for(marker in markers_list){
    a_file_name =paste(name_of_file,marker,sep = "")
    
    #read first 5 rows of the table
    tab5rows <- read.table(file = a_file_name, header = FALSE, nrows = 5,row.names = 1)
    
    #get the classes of the columns of the table
    classes <- sapply(tab5rows, class)
    #read the rest (or a bigger part) of the table knowing the class of the columns speeds up the reading process
    
    
    data_file <- read.table(a_file_name, header = FALSE, colClasses = classes, row.names = 1)
    marker_data = as.matrix(paste(data_file$V2,data_file$V3,sep = "/"))
    bind_markers = cbind(bind_markers,marker_data)
  }
  
  # note ncode set to 3 - alleles are considered to be coded by three characters
  markers_df2genind = df2genind(X = bind_markers, ploidy = 2, ncode = 3,ind.names = rownames(data_file),sep = "/")
  
  # define pops
  markers_df2genind@pop = data_file$V6
  
  # define coords
  xy_coords = cbind(data_file$V4,data_file$V5)
  
  markers_df2genind@other$xy = xy_coords
  
  # define gen
  markers_df2genind@other$generation = current_generation
  
  # return genind object
  return(markers_df2genind)
}

################
#GET DATA FROM FILES
################
#+++++++++++++++
Get_data_from_files = function(the_file_names, set_markers,generations_list){
  # Gets data from the files with the parameters gotten from the previous functions
  the_data_list=NULL
  for(i in the_file_names){
    the_data_list = c(the_data_list,Read_data_func(i, set_markers))
  }
  #set names of elements in list to the generations that they belong to
  names(the_data_list) = generations_list
  return(the_data_list)
}




###########################################################################
###########################################################################
#++++SECTION 2 - FUNCTIONS TO ANALYSE DATA 
###########################################################################
###########################################################################


################
#DO SUMMARY CALCS ON GENIND OBJ
################
#+++++++++++++++
Do_summary_calcs = function(the_genind_obj){
  
  summary_genind = summary(the_genind_obj)
  names(summary_genind)
  
  par(mfrow=c(2,2),oma=c(0,0,1,0))#sets 2x2 spots for the graphics and sets space for outer title
  
  plot(summary_genind$n.by.pop, summary_genind$pop.n.all, xlab="Colonies sample size",
       ylab="Number of alleles",main="Alleles numbers and sample sizes",
       type="n",
       xlim=c(range(summary_genind$n.by.pop)[1]-5,range(summary_genind$n.by.pop)[2]+5),
       ylim=c(range(summary_genind$pop.n.all)[1]-5,range(summary_genind$pop.n.all)[2]+5))
  text(summary_genind$n.by.pop,summary_genind$pop.n.all,lab=names(summary_genind$n.by.pop))
  barplot(summary_genind$loc.n.all, ylab="Number of alleles",
          main="Number of alleles per locus")
  barplot(summary_genind$Hexp-summary_genind$Hobs, main="Heterozygosity: expected-observed",
          ylab="Hexp - Hobs")
  barplot(summary_genind$n.by.pop, main="Sample sizes per population",
          ylab="Number of genotypes",las=3)
  title(paste("Generation: ",the_genind_obj@other$generation), outer=TRUE)
  
  the_plot = recordPlot()
  
  bartlett_test = bartlett.test(list(summary_genind$Hexp,summary_genind$Hobs))
  
  the_ttest = t.test(summary_genind$Hexp,summary_genind$Hobs,pair=T,var.equal=TRUE,alter="greater")
  
  results = list(summary_genind,bartlett_test,the_ttest,the_plot)
  
  
  #return(summary_genind)
  return(results)
}

################
#COMPUTE HS ON GENIND ARRAY
################
#+++++++++++++++
Compute.Hs.make.df = function (genind_array){
  compute.Hs.df = data.frame()
  for(i in 1:length(genind_array)){
    # if there are more than 1 allele per locus compute Hs
    if(length(genind_array[[i]]@all.names$loc1)>1 || length(genind_array[[i]]@all.names$loc2)>1){
      # Stack helps shape the data so that we can plot it later
      test.var = stack(Hs(genind_array[[i]]))
      test.var.cbind = cbind(test.var,gen = genind_array[[i]]@other$generation)
      compute.Hs.df = rbind(compute.Hs.df,test.var.cbind)
    }
  }
  compute.Hs.df.numeric.gen = compute.Hs.df
  ## transform gen col from factor to numeric so that we are able to use gen as a continuous scale for the plot
  compute.Hs.df.numeric.gen$gen = as.numeric(as.character(compute.Hs.df.numeric.gen$gen))
  return(compute.Hs.df.numeric.gen)
}

################
#TRANSFORM INTO SIMPLE POP (SO THAT WE CAN THEN COMPUTE HS ON GENIND ARRAY)
################
#+++++++++++++++
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

################
#MAKE HS PLOTS
################
#+++++++++++++++
Make.HS.plots = function(compute.Hs.df.numeric.gen, isSimple = F){
  if(!isSimple){
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
  } else {
    ggplot(data = compute.Hs.df.numeric.gen)+
      geom_line(aes(x=gen, y=values,group=ind,color=ind), ## group by ind, color by ind
                subset(compute.Hs.df.numeric.gen,
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
    
  }
}

################
#GET GENIND OBJ FROM A GIVEN GENERATION
################
#+++++++++++++++
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

################
#COMPUTE FST PLOT
################
#+++++++++++++++
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

################
#COMPUTE PAIRWISE_FST/DISTANCE (IBD) PLOT
################
#+++++++++++++++
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
    scale_x_continuous(expand=c(0.0,0.0), limits=c(0,12),breaks = seq(0,12,1)) + #TODO change limits to dynamic values
    scale_y_continuous(expand=c(0.0,0.0), limits=c(0,0.4)) +
    coord_cartesian(xlim = c(0, 12), ylim = c(0, 0.4))+
    ylab("Paiwise Fst")+
    xlab("Geographical distance between demes")+
    ggtitle(paste("Generation ",raw_data_ind@other$generation))+
    cowplot::theme_cowplot()
  
  return(the_plot)
}


###########################################################################
###########################################################################
#++++SECTION 3 - DO DATA ANALYSIS 
###########################################################################
###########################################################################

# set working directory (where data is located)
setwd("/home/tiago/server_folder/dist/output/server_sim_Test_Clusters_discrete_deme/server_sim_discrete_deme/Adegenet_sim1/")

# set parameters
the_generation_list = Set_generations_list(10,6000,5)
the_layers_list = Set_layers_list("layer0")
the_markers_list = Set_markers_list("A1","A2")

# get names of the files and then get data from these files
the_file_names = Get_file_names(the_generation_list, the_layers_list)
raw_data = Get_data_from_files(the_file_names = the_file_names,
                               set_markers = the_markers_list,
                               generations_list = the_generation_list)


# perform analysis #######


sample.summary = Do_summary_calcs(the_genind_obj = Get_genind_from_gen(a_genind_array = raw_data,generation = 6000))

sample.HS = Compute.Hs.make.df(genind_array = raw_data)
sample.HS.simplePop = Transform_into_simple_pop(genind_array_obj = raw_data)
sample.HS.simplePop = Compute.Hs.make.df(genind_array = sample.HS.simplePop)
sample.HS = Make.HS.plots(compute.Hs.df.numeric.gen = sample.HS,isSimple = F)
sample.HS.simplePop = Make.HS.plots(compute.Hs.df.numeric.gen = sample.HS.simplePop,isSimple = T)

sample.Fst = Compute.Fst.plot(gen.ind.object = Get_genind_from_gen(a_genind_array = raw_data,generation = 6000),
                 matrix.as.triangle = T,
                 upper.triangle = F,
                 has.in.graph.txt = F)
sample.pFstDist = Make_pairFst_distance_plots(genind_array = raw_data,generation = 6000)



