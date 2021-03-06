if(!require("devtools")){
  install.packages("devtools", dependencies = TRUE)
  install_github("thibautjombart/adegenet")
  library("adegenet")
  }


  
if(!require("ade4")) install.packages("ade4", dependencies = TRUE)
if(!require("hierfstat")) install.packages("hierfstat", dependencies = TRUE)
if(!require("pegas")) install.packages("pegas", dependencies = TRUE)
if(!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)
if(!require("reshape2")) install.packages("reshape2", dependencies = TRUE)
if(!require("gridExtra")) install.packages("gridExtra", dependencies = TRUE)
#if(!require("cowplot")) install.packages("cowplot", dependencies = TRUE)
#if(!require("plotly")) install.packages("plotly", dependencies = TRUE)
if(!require("Hmisc")) install.packages("Hmisc", dependencies = TRUE)
if(!require("plyr")) install.packages("plyr", dependencies = TRUE)
if(!require("PopGenReport")) install.packages("PopGenReport", dependencies = TRUE)
if(!require("data.table")) install.packages("data.table", dependencies = TRUE)
  
  #library(adegenet)
  #library(ade4)
  #library(hierfstat)
  #library(pegas)
  #library(ggplot2) # to plot
  #library(reshape2) # to melt/reshape tables
  #library(gridExtra) # to put several plots into 1 image, easily
  #library(cowplot) # publication ready plots?
  #library(plotly) # interactive plots
  
  ###########################################################################
  ###########################################################################
  #++++SECTION 1 - FUNCTIONS TO GET DATA
  ###########################################################################
  ###########################################################################
  
  ################
  #SET NUMBER OF SIMULATIONS
  ################
  #+++++++++++++++
  Set_number_of_simulations = function(number_of_simulations){
    return(number_of_simulations)
  }
  
  
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
  Read_data_func = function(name_of_file, markers_list, is_marker_diploid = TRUE, simulation_ID = 0){
  
    #name_of_file = the_file_names[[9]]
    #markers_list = c("A1","A2","A3","A4","A5","A6","A7","A8","A9","A10")
    #is_marker_diploid=T
    #simulation_ID=1
    #setwd(dir = paste(the_path_to_data,"Adegenet_sim_1","/",sep = ""))
    
    if(FALSE){
      #***fread is much faster than an optimized read.table***
    #read first 5 rows of the table
    #tab5rows <- read.table(file = paste(name_of_file,markers_list[1],sep = ""), header = FALSE, nrows = 100,row.names = 1,na.strings = "NA")
    
    #get the classes of the columns of the table
    #tabClasses <- sapply(tab5rows, class)
    }
    data_file <- fread(input = paste(name_of_file,markers_list[1],sep = ""),encoding = "UTF-8", header = FALSE, na.strings = "NA",data.table = F)
    
    # get the number of rows that the table will have, use it later to speed up process
    # not sure if it is that much faster, consider erasing it
    sizeOfTable = nrow(data_file)
    
    markerCondition = markers_list == "MT" | markers_list == "Y"
    
    current_generation = NULL
    current_generation = strsplit(x=name_of_file,split = "_")[[1]][1]
    
    data_file = NULL
    marker_data.haplo = NULL
    marker_data.diplo = NULL

    a_file_name = NULL
    
    #bind_markers.diplo = NULL
    bind_markers.diplo = matrix(NA, nrow=sizeOfTable, ncol=length(markers_list))
    #bind_markers.haplo = NULL
    bind_markers.haplo = matrix(NA, nrow=sizeOfTable, ncol=length(markers_list))
    
    for(marker in 1:length(markers_list)){
      
      a_file_name = paste(name_of_file,markers_list[marker],sep = "")
      
      #read the rest (or a bigger part) of the table knowing the class of the columns speeds up the reading process
      #data_file <- read.table(a_file_name, header = FALSE, colClasses = tabClasses, nrows = sizeOfTable, row.names = 1,na.strings = "NA")
      #fread is much faster than an optimized read.table
      data_file <- fread(input = a_file_name, header = FALSE, encoding = "UTF-8", na.strings = "NA",data.table = F)
      
      #if marker is Y or MTdna then it is haploid data, save that markers col (instead of 2 cols for diploid)
      if(markerCondition[marker]){
        marker_data.haplo = as.matrix(data_file$V2)
        #bind_markers.haplo = cbind(bind_markers.haplo,marker_data.haplo)
        bind_markers.haplo[,marker] = marker_data.haplo
      }else{
        #marker_data.diplo = as.matrix(paste(data_file$V2,data_file$V3,sep = "/"))
        
        marker_data.diplo = as.vector(data_file$V2)
        
        
        #length(marker_data.diplo)
        #bind_markers.diplo = cbind(bind_markers.diplo,marker_data.diplo)
        bind_markers.diplo[,marker] = marker_data.diplo
      }
    }
    
    
    
    if(is_marker_diploid){
      # note ncode set to 3 - alleles are considered to be coded by three characters
      markers.diplo_df2genind = df2genind(X = bind_markers.diplo, ploidy = 2, ncode = 3,ind.names = data_file$V1, sep = "/") 
      # define pops
      markers.diplo_df2genind@pop = factor(data_file$V5)
      # define coords
      xy_coords = cbind(data_file$V3,data_file$V4)
      markers.diplo_df2genind@other$xy = xy_coords
      # define gen
      markers.diplo_df2genind@other$generation = current_generation
      # define simulation ID
      markers.diplo_df2genind@other$simulation.ID = simulation_ID
      # return genind object
      return(markers.diplo_df2genind)
    }else{
      markers.haplo_df2genind = df2genind(X = bind_markers.haplo, ploidy = 1, ncode = 3,ind.names = data_file$V1, sep = "/", NA.char = "NA")  
      # define pops
      markers.haplo_df2genind@pop = data_file$V5
      # define coords
      xy_coords = cbind(data_file$V3,data_file$V4)
      markers.haplo_df2genind@other$xy = xy_coords
      # define gen
      markers.haplo_df2genind@other$generation = current_generation
      # define simulation ID
      markers.diplo_df2genind@other$simulation.ID = simulation_ID
      # return genind object
      return(markers.haplo_df2genind)
    }
    
  }
  
  ################
  #GET DATA FROM FILES FROM A SINGLE SIMULATION
  ################
  #+++++++++++++++
  Get_data_from_files_single_sim = function(the_file_names,
                                            set_markers,generations_list,
                                            path_to_data, simulation_number,
                                            is_marker_diploid=TRUE){
    
    # Gets data from the files with the parameters gotten from the previous functions
    
      setwd(paste(path_to_data,"Adegenet_sim_",simulation_number,"/",sep=""))
    
      the_data_list=NULL
      for(i in the_file_names){
        
        the_data_list = c(the_data_list,Read_data_func(i, set_markers,is_marker_diploid = is_marker_diploid, simulation_ID = simulation_number))
      }
      #set names of elements in list to the generations that they belong to
      names(the_data_list) = generations_list
    
  
    return(the_data_list)
  }
  
  ################
  #GET DATA FROM FILES FROM MULTIPLE SIMULATIONS
  ################
  #+++++++++++++++
  Get_data_from_files_multi_sim = function(the_file_names,
                                           set_markers,generations_list, 
                                           path_to_data, 
                                           number_of_simulations,
                                           is_marker_diploid=TRUE){
    
    # Gets data from the files with the parameters gotten from the previous functions
    # Returns list of lists of genind objects
    
    ####test data
    #set_markers = c("A1")
    #generations_list = seq(0,100,10)
    #number_of_simulations=1
    #the_file_names= Get_file_names(generations_list,"layer0")
    #path_to_data = the_path_to_data
    #is_marker_diploid = TRUE
    ####
    
    all_sim=vector(mode = "list",length = number_of_simulations)
    
    
    
    for(j in 1:number_of_simulations){
      setwd(paste(path_to_data,"/","Adegenet_sim_",j,"/",sep=""))
      
      the_data_list = vector(mode = "list",length = length(the_file_names))
      for(i in 1:length(the_file_names)){
        
        #the_data_list = c(the_data_list,Read_data_func(i, set_markers,is_marker_diploid=is_marker_diploid, simulation_ID = j))
        the_data_list[[i]] <- Read_data_func(the_file_names[i], set_markers,is_marker_diploid=is_marker_diploid, simulation_ID = j)
        #the_data_list = lapply(the_file_names,FUN = Read_data_func, name_of_file = the_file_names, markers_list = set_markers, is_marker_diploid = T, simulation_ID = j)
        
      }
      #set names of elements in list to the generations that they belong to
      names(the_data_list) = generations_list
      
      #all_sim = cbind(all_sim,the_data_list)  
      #all_sim = c(all_sim,the_data_list)
      
      #all_sim[[length(all_sim)+1]] = the_data_list
      all_sim[[j]] = the_data_list
    }
    
    ###test
    #str(list(the_data_list,the_data_list),max.level = 1)
    #str(all_sim, max.level = 1)
    #dim(all_sim)
    #class(all_sim)
    #all_sim[[2]][[5]]@other$generation
    ###
    
    
    return(all_sim)
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
  Do_summary_calcs = function(simulation_name, the_genind_obj, the_path_to_plot_folder){
    
    
    summary_genind = summary(the_genind_obj)
    names(summary_genind)
    
    
    pdf(file = paste(the_path_to_plot_folder,"/",simulation_name,"_SummaryStats.pdf",sep=""),width = 20, height = 15)
    par(mfrow=c(2,2),oma=c(0,0,1,0))#sets 2x2 spots for the graphics and sets space for outer title
    
    plot(summary_genind$n.by.pop, summary_genind$pop.n.all, xlab="Colonies sample size",
         ylab="Number of alleles",main="Alleles numbers and sample sizes",
         type="n",
         xlim=c(range(summary_genind$n.by.pop)[1]-5,range(summary_genind$n.by.pop)[2]+5),
         ylim=c(range(summary_genind$pop.n.all)[1]-5,range(summary_genind$pop.n.all)[2]+5))
    text(summary_genind$n.by.pop,summary_genind$pop.n.all,lab=names(summary_genind$n.by.pop))
    barplot(summary_genind$loc.n.all, ylab="Number of alleles",
            main="Number of alleles per locus")
    barplot(summary_genind$Hexp-summary_genind$Hobs,
            ylab="Hexp - Hobs",col=ifelse(summary_genind$Hexp-summary_genind$Hobs>0,"red","blue"))
    title(expression(phantom("Heterozygosity: ") * "expected" * phantom(" - observed")),col.main="red")
    title(expression(phantom("Heterozygosity: expected - ") * "observed"),col.main="blue")
    title(expression("Heterozygosity:" * phantom(" expected ") * "-" * phantom(" observed"),col.main="black"))
    
    barplot(summary_genind$n.by.pop, main="Sample sizes per population",
            ylab="Number of genotypes",las=3)
    title(paste("Generation: ",the_genind_obj@other$generation), outer=TRUE)
  
    
    the_plot = recordPlot()
    dev.off()
    
    bartlett_test = bartlett.test(list(summary_genind$Hexp,summary_genind$Hobs))
    
    the_ttest = t.test(summary_genind$Hexp,summary_genind$Hobs,pair=T,var.equal=TRUE,alter="greater")
    
    results = list(summary_genind,bartlett_test,the_ttest,the_plot)
    
    print("Summary Calcs done.")
    #return(summary_genind)
    return(results)
  }
  
  Do_summary_calcs_updated = function(simulation_name, the_genind_array, the_path_to_plot_folder, groupBy="simulation", simGroup=NULL, genGroup=NULL){
    
    if(FALSE){
      
      simulation_name=name_of_simulation
      the_genind_array=raw_data_multi_sim
      the_path_to_plot_folder=the_path_to_plot_folder
      groupBy="simulation"
      simGroup=1
      genGroup=100
      
      
    }
    
    if(FALSE){
    mySummary= function(object, verbose = TRUE, ...){
      x <- object
      if(!is.genind(x)) stop("Provided object is not a valid genind.")
      
      
      if(is.null(pop(x))){
        pop(x) <- rep("P1", nInd(x))
      }
      
      ## BUILD THE OUTPUT ##
      ## type-independent stuff
      res <- list()
      
      res$n <- nrow(myTab(x))
      
      res$n.by.pop <- as.numeric(table(pop(x)))
      names(res$n.by.pop) <- popNames(x)
      
      ## PA case ##
      if(x@type=="PA"){
        ## % of missing data
        res$NA.perc <- 100*sum(is.na(myTab(x)))/prod(dim(myTab(x)))
        
        return(invisible(res))
      }
      
      
      ## codom case ##
      res$loc.n.all <- nAll(x)
      
      temp <- myTab(genind2genpop(x,quiet=TRUE))
      
      res$pop.n.all <- apply(temp,1,function(r) sum(r!=0,na.rm=TRUE))
      
      res$NA.perc <- 100*(1-mean(propTyped(x,by="both")))
      
      ## handle heterozygosity
      if(any(ploidy(x) > 1)){
        ## auxiliary function to compute observed heterozygosity
        temp <- lapply(seploc(x),tab, freq=TRUE)
        f1 <- function(tab){
          H <- apply(tab, 1, function(vec) any(vec > 0 & vec < 1))
          H <- mean(H,na.rm=TRUE)
          return(H)
        }
        
        res$Hobs <- unlist(lapply(temp,f1))
        
        ## auxiliary function to compute expected heterozygosity
        ## freq is a vector of frequencies
        f2 <- function(freq){
          H <- 1-sum(freq*freq,na.rm=TRUE)
          return(H)
        }
        
        temp <- genind2genpop(x,pop=rep(1,nInd(x)),quiet=TRUE)
        temp <- myTab(temp, freq=TRUE, quiet=TRUE)
        res$Hexp <-tapply(temp^2, locFac(x), function(e) 1-sum(e, na.rm=TRUE))
      } else { # no possible heterozygosity for haploid genotypes
        res$Hobs <- 0
        res$Xexp <- 0
      }
      
      ## add class and return
      class(res) <- "genindSummary"
      return(res)
    }  # end summary.genind
    }
    
    allData_summary = rapply(object = the_genind_array,f = function(X){return(summary(X))}, how = 'list')
    
    #allData_meanSummary = rapply(allData_summary, mean, how='list')
    
    allData_meanSummary = rapply(allData_summary, mean, how='replace')
    lapply(allData_meanSummary, function(x) split(x,f = '[['))
    
    #gen - stat
    
    
    lapply(X = allData_meanSummary,FUN = lapply,'[[',7)
    
    
    allData_meanSimSummary=NULL
    allData_meanSimSummary$Hexp=mean(unlist(lapply(lapply(allData_meanSummary, '[[',10), '[[',7)))
    
    asd = mean(unlist(lapply(lapply(allData_meanSummary, '[[',10), '[[',7)))
    str(asd)
    
    
    lapply(allData_meanSummary, function(x){lapply(x,FUN ='[[',2)},mean)
    
    allData_meanSummary_1 = rbindlist(allData_meanSummary,fill=T, use.names = T)
    
    
    unlist(allData_meanSummary[[1]][[10]], recursive = F, use.names = T)
    
    str(allData_meanSummary,max.level = 2)
    
    allData_meanSummary
    
    as.data.frame(allData_meanSummary[[1]][[1]]$pop.n.all)
    
    ggplot(allData_meanSummary)
    
    
    pdf(file = paste(the_path_to_plot_folder,"/",simulation_name,"_SummaryStats.pdf",sep=""),width = 20, height = 15)
    par(mfrow=c(2,2),oma=c(0,0,1,0))#sets 2x2 spots for the graphics and sets space for outer title
    
    plot(summary_genind$n.by.pop, summary_genind$pop.n.all, xlab="Colonies sample size",
         ylab="Number of alleles",main="Alleles numbers and sample sizes",
         type="n",
         xlim=c(range(summary_genind$n.by.pop)[1]-5,range(summary_genind$n.by.pop)[2]+5),
         ylim=c(range(summary_genind$pop.n.all)[1]-5,range(summary_genind$pop.n.all)[2]+5))
    text(summary_genind$n.by.pop,summary_genind$pop.n.all,lab=names(summary_genind$n.by.pop))
    barplot(summary_genind$loc.n.all, ylab="Number of alleles",
            main="Number of alleles per locus")
    barplot(summary_genind$Hexp-summary_genind$Hobs,
            ylab="Hexp - Hobs",col=ifelse(summary_genind$Hexp-summary_genind$Hobs>0,"red","blue"))
    title(expression(phantom("Heterozygosity: ") * "expected" * phantom(" - observed")),col.main="red")
    title(expression(phantom("Heterozygosity: expected - ") * "observed"),col.main="blue")
    title(expression("Heterozygosity:" * phantom(" expected ") * "-" * phantom(" observed"),col.main="black"))
    
    barplot(summary_genind$n.by.pop, main="Sample sizes per population",
            ylab="Number of genotypes",las=3)
    title(paste("Generation: ",the_genind_obj@other$generation), outer=TRUE)
    
    
    the_plot = recordPlot()
    dev.off()
    
    bartlett_test = bartlett.test(list(summary_genind$Hexp,summary_genind$Hobs))
    
    the_ttest = t.test(summary_genind$Hexp,summary_genind$Hobs,pair=T,var.equal=TRUE,alter="greater")
    
    results = list(summary_genind,bartlett_test,the_ttest,the_plot)
    
    print("Summary Calcs done.")
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
      print(i)
      # if there are more than 1 allele per locus compute Hs
      #if(length(genind_array[[i]]@all.names$loc1)>1 || length(genind_array[[i]]@all.names$loc2)>1){
        # Stack helps shape the data so that we can plot it later
        test.var = stack(Hs(genind_array[[i]]))
        test.var.cbind = cbind(test.var,gen = genind_array[[i]]@other$generation)
        compute.Hs.df = rbind(compute.Hs.df,test.var.cbind)
      #}
    }
    compute.Hs.df.numeric.gen = compute.Hs.df
    ## transform gen col from factor to numeric so that we are able to use gen as a continuous scale for the plot
    compute.Hs.df.numeric.gen$gen = as.numeric(as.character(compute.Hs.df.numeric.gen$gen))
    return(compute.Hs.df.numeric.gen)
  }
  
  ################
  #COMPUTE HS OVER SEVERAL SIMULATIONS
  ################
  #+++++++++++++++
  Compute.Hs.over.all.sims = function (sim_data){
  
    # Note that if for some reasons all the individuals are homozygous while analysing MORE than 1
    # marker, then this function will crash because for some reason, while analyzing 1 marker it can
    # compute that the heterozygosity can be 0, but while calculating the same for several marker
    # and if for one/several markers all individuals are homozygous, then the function will crash.
    all_sims=NULL
    
    for(j in 1:length(sim_data)){
      genind_array = sim_data[[j]]
      
      compute.Hs.df = data.frame()
      for(i in 1:length(genind_array)){
      # if there are more than 1 allele per locus compute Hs
        #if(length(genind_array[[i]]@all.names$loc1)>1 || length(genind_array[[i]]@all.names$loc2)>1){
          # Stack helps shape the data so that we can plot it later
        
        
        test.var = stack(Hs(genind_array[[i]]))
          
          test.var.cbind = cbind(test.var,gen = genind_array[[i]]@other$generation)
          compute.Hs.df = rbind(compute.Hs.df,test.var.cbind)
        #}
      }
      compute.Hs.df.numeric.gen = compute.Hs.df
      ## transform gen col from factor to numeric so that we are able to use gen as a continuous scale for the plot
      compute.Hs.df.numeric.gen$gen = as.numeric(as.character(compute.Hs.df.numeric.gen$gen))
    
      all_sims[[j]] = compute.Hs.df.numeric.gen$values
    }
    
    all_sims = do.call("cbind",all_sims)
  
    all_sims = data.frame(all_sims)
    all_sims$mean = apply(all_sims,MARGIN = 1,FUN = mean)
    compute.Hs.df.numeric.gen$gen
    all_sims$ind = compute.Hs.df.numeric.gen$ind
    all_sims$gen = compute.Hs.df.numeric.gen$gen
    
    return(all_sims)
  }
  
  ################
  #COMPUTE HS (PER LOCUS) OVER SEVERAL SIMULATIONS
  ################
  #+++++++++++++++
  Compute.Hs.over.all.sims.per.locus = function (sim_data){
    
    # Note that if for some reasons all the individuals are homozygous while analysing MORE than 1
    # marker, then this function will crash because for some reason, while analyzing 1 marker it can
    # compute that the heterozygosity can be 0, but while calculating the same for several marker
    # and if for one/several markers all individuals are homozygous, then the function will crash.
    
    #sim_data = raw_data_multi_sim
    all_sims=NULL
    
    for(j in 1:length(sim_data)){
      genind_array = sim_data[[j]]
      
      compute.Hs.df = data.frame()
      for(i in 1:length(genind_array)){
        # if there are more than 1 allele per locus compute Hs
        #if(length(genind_array[[i]]@all.names$loc1)>1 || length(genind_array[[i]]@all.names$loc2)>1){
        # Stack helps shape the data so that we can plot it later
        
        ###TESTING
        hs.sep.loc = lapply(X = seploc(genind_array[[i]]),FUN = Hs) #apply Hs to list
        hs.sep.loc.df = data.frame(hs.sep.loc)
        hs.sep.loc.df$values = apply(hs.sep.loc.df,MARGIN = 1,FUN = mean) #create new col "values" with the mean over all the loci
        new.test.var = subset(hs.sep.loc.df,select = "values") # new df w/ only "values" col
        new.test.var$ind = rownames(new.test.var) # transform row names into a collumn
        rownames(new.test.var) = NULL # transform row names into numbers
        ###/TESTING
        
        #test.var = stack(Hs(genind_array[[i]]))
        test.var = new.test.var
        
        test.var.cbind = cbind(test.var,gen = genind_array[[i]]@other$generation)
        compute.Hs.df = rbind(compute.Hs.df,test.var.cbind)
        #}
      }
      compute.Hs.df.numeric.gen = compute.Hs.df
      ## transform gen col from factor to numeric so that we are able to use gen as a continuous scale for the plot
      compute.Hs.df.numeric.gen$gen = as.numeric(as.character(compute.Hs.df.numeric.gen$gen))
      
      all_sims[[j]] = compute.Hs.df.numeric.gen$values
    }
    
    all_sims = do.call("cbind",all_sims)
    
    all_sims = data.frame(all_sims)
    all_sims$mean = apply(all_sims,MARGIN = 1,FUN = mean)
    #compute.Hs.df.numeric.gen$gen
    all_sims$ind = compute.Hs.df.numeric.gen$ind
    all_sims$gen = compute.Hs.df.numeric.gen$gen
    
    return(all_sims)
  }
  
  ################
  #TRANSFORM INTO SIMPLE POP (SO THAT WE CAN THEN COMPUTE HS ON GENIND ARRAY)
  ################
  #+++++++++++++++
  Transform_into_simple_pop = function(genind_array_obj){
    # Transforms a genind array grouping populations into simpler(read broader) groups
    # pop1A will be pop1
    # its defined for the following population code up to pop15
    Group.demes.function = function(X){
      ## Group demes according to their population code (eg "pop1A" will be "pop1", "pop6E" will be "pop6")
      X = raw_data_multi_sim[[1]][[1]]
      popNames(X)
      varN = NULL
      for(i in popNames(X)){
        
        if(strsplit(i,split = "pop1{1}")[[1]][1] == ""){
          varN = c(varN, "pop1")
        }
      }
      
      
      varN = NULL
      for(i in popNames(X)){
        
        if(strsplit(i,split = "pop1{1}")[[1]][1] == ""){
          varN = c(varN, "pop1")
        } else if (strsplit(i,split = "pop2{1}")[[1]][1] == ""){
          varN = c(varN, "pop2")
        } else if (strsplit(i,split = "pop3{1}")[[1]][1] == ""){
          varN = c(varN, "pop3")
        } else if (strsplit(i,split = "pop4{1}")[[1]][1] == ""){
          varN = c(varN, "pop4")
        } else if (strsplit(i,split = "pop5{1}")[[1]][1] == ""){
          varN = c(varN, "pop5")
        } else if (strsplit(i,split = "pop6{1}")[[1]][1] == ""){
          varN = c(varN, "pop6")
        } else if (strsplit(i,split = "pop7{1}")[[1]][1] == ""){
          varN = c(varN, "pop7")
        } else if (strsplit(i,split = "pop8{1}")[[1]][1] == ""){
          varN = c(varN, "pop8")
        } else if (strsplit(i,split = "pop9{1}")[[1]][1] == ""){
          varN = c(varN, "pop9")
        } else if (strsplit(i,split = "pop10{1}")[[1]][1] == ""){
          varN = c(varN, "pop10")
        } else if (strsplit(i,split = "pop11{1}")[[1]][1] == ""){
          varN = c(varN, "pop11")
        } else if (strsplit(i,split = "pop12{1}")[[1]][1] == ""){
          varN = c(varN, "pop12")
        } else if (strsplit(i,split = "pop13{1}")[[1]][1] == ""){
          varN = c(varN, "pop13")
        } else if (strsplit(i,split = "pop14{1}")[[1]][1] == ""){
          varN = c(varN, "pop14")
        } else if (strsplit(i,split = "pop15{1}")[[1]][1] == ""){
          varN = c(varN, "pop15")
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
  Make.HS.plots = function(compute.Hs.df.numeric.gen, isSimple = F, initial.gen=0, env.change=0, last.gen){
    if(!isSimple){
    hs.plot = ggplot(data = compute.Hs.df.numeric.gen)+
      geom_line(aes(x=gen, y=values,group=ind,color=ind), ## group by ind, color by ind
      alpha=1.0)+ 
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank())+
      ggtitle("Expected Heterozigosity for each deme")+
      xlab("Generations")+
      ylab("Expected Heterozigosity")+
      theme(legend.position="right")+
      scale_x_continuous(breaks=seq(initial.gen,last.gen,1000))+
      #scale_x_discrete(labels=compute.Hs.df$gen,breaks=compute.Hs.df$gen, limits = compute.Hs.df$gen)+
      geom_vline(xintercept = env.change,linetype = "longdash")+
      coord_cartesian(ylim=c(0,1), xlim=c(initial.gen,last.gen))
    } else {
      hs.plot = ggplot(data = compute.Hs.df.numeric.gen)+
        geom_line(aes(x=gen, y=values,group=ind,color=ind), ## group by ind, color by ind
                  subset(compute.Hs.df.numeric.gen,
                         ind=="pop1" | ind=="pop2" | ind=="pop3" | 
                           ind=="pop4" | ind=="pop5" | ind=="pop6" |
                         ind=="pop7" | ind=="pop8" | ind=="pop9"), alpha=1.0)+ ## subset the data, choose which ind(pops) we want to plot
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank())+
        ggtitle("Expected Heterozigosity for group of demes")+
        xlab("Generations")+
        ylab("Expected Heterozigosity")+
        theme(legend.position="right")+
        scale_x_continuous(breaks=seq(initial.gen,last.gen,1000))+
        #scale_x_discrete(labels=compute.Hs.df$gen,breaks=compute.Hs.df$gen, limits = compute.Hs.df$gen)+
        geom_vline(xintercept = env.change,linetype = "longdash")+
        coord_cartesian(ylim=c(0,1), xlim=c(initial.gen,last.gen))
      
    }
  
    print("HS plots done.")
    return(hs.plot)
  }
  
  ################
  #MAKE HS MEAN PLOT (MEAN OVER ALL SIMS)
  ################
  #+++++++++++++++
  Make.HS.mean.plots = function(compute.Hs.all.sims, initial.gen=0, env.change=0, last.gen, subtitle=NULL, showPerDeme=FALSE, number_of_simulations){
  
    summaryAllPop = summarySE(data = compute.Hs.all.sims, measurevar = "mean", groupvars = c("gen"))
    
      hs.plot = ggplot(summaryAllPop,aes(x = gen, y = mean))+geom_line()+geom_point()+
        geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
        ggtitle(label = paste("Mean Expected Heterozigosity over " ,number_of_simulations, " simulations", sep = ""),subtitle = subtitle)+
        xlab("Generations")+
        ylab("Expected Heterozigosity")+
        theme(legend.position="right")+
        scale_x_continuous(breaks=seq(initial.gen,last.gen,1000))+
        #scale_x_discrete(labels=compute.Hs.df$gen,breaks=compute.Hs.df$gen, limits = compute.Hs.df$gen)+
        geom_vline(xintercept = env.change,linetype = "longdash")+
        coord_cartesian(ylim=c(0,1), xlim=c(initial.gen,last.gen))
        #theme(plot.title = element_text(family = 'Helvetica',color = '#666666',face = 'bold',size = 18,hjust = 0.5))
    if(showPerDeme){
      hs.plot = hs.plot +geom_line(data = compute.Hs.all.sims,aes(x = gen, y = mean,
                                                                  #group=ind, 
                                                                  color=ind
                                                                  ),alpha=0.5)+
        #facet_grid(.~ind)+
        scale_colour_discrete(name = "Demes",guide=F)
        
    }
    print("HS mean plot done.")
    return(hs.plot)
  }
  

  ################
  #TOTAL NUMBER OF ALLELES PER LOCUS VS AVERAGE NUMBER OF ALLELES PER LOCUS PER POP (MEAN OVER ALL SIMS)
  ################
  #+++++++++++++++
  Tot.vs.Avg.N.Alleles = function(all.simulations.data, avg.over.all.loci=F){
    
    #all.simulations.data = raw_data_multi_sim
    
    save.generations = NULL
    mean.nbAllelePerSampledDeme = NULL
    tot.num.alleles = NULL
    
    for(j in 1:length(all.simulations.data[[1]])){
      for(i in 1:length(all.simulations.data)){  
        
        # t() transposes the matrix, that is, transforms columns in rows and viceversa  
        nbAllelePerSampledDeme = data.frame(t(nb.alleles(all.simulations.data[[i]][[j]])))
        mean.nbAllelePerSampledDeme = rbind(mean.nbAllelePerSampledDeme,colMeans(nbAllelePerSampledDeme))
        
        tot.num.alleles =rbind(tot.num.alleles, nAll(all.simulations.data[[i]][[j]]),deparse.level = 0)

        save.generations = rbind(save.generations, all.simulations.data[[i]][[j]]@other$generation)
        
      }
    }  
    
    
    numberOfLoci = length(levels(all.simulations.data[[1]][[1]]@loc.fac))
    
    mean.nbAllelePerSampledDeme = data.frame(mean.nbAllelePerSampledDeme)
    mean.nbAllelePerSampledDeme$generation = as.numeric(save.generations)
    
    mean.nbAllelePerSampledDeme.simAvg = aggregate(mean.nbAllelePerSampledDeme[,1:numberOfLoci], list(Generation = mean.nbAllelePerSampledDeme$generation), mean)
    
    
    tot.num.alleles = data.frame(tot.num.alleles)
    tot.num.alleles$generation = as.numeric(save.generations)
    
    tot.num.alleles.simAvg = aggregate(tot.num.alleles[,1:numberOfLoci], list(Generation = tot.num.alleles$generation), mean)
    
    #### TODO SAVE IMAGE TO FILE
    if(FALSE){
      ggplot2::ggsave(file=paste(the_path_to_plot_folder,"/",
                                 name_of_simulation,
                                 "_TotNumAllelesVSAvgNumAlleles",
                                 length(raw_data_multi_sim),
                                 "Sims.pdf",sep=""), width = 20, height = 15)
    }
    
    if(avg.over.all.loci){
      df.simAlleleAvg = data.frame(avg=rowMeans(mean.nbAllelePerSampledDeme.simAvg[,2:numberOfLoci+1]))
      df.simAlleleTotal = data.frame(tot=rowMeans(tot.num.alleles.simAvg[,2:numberOfLoci+1]))
      df.Avg.Vs.Total = data.frame(df.simAlleleAvg,df.simAlleleTotal)
      df.Avg.Vs.Total$gen = mean.nbAllelePerSampledDeme.simAvg$Generation
      return(ggplot(df.Avg.Vs.Total)+
               geom_line(aes(x = gen,y=avg,color="average"))+
               geom_line(aes(x = gen,y=tot,color="total"))+
               xlab("Generations")+
               ylab("Nb Alleles"))
    }
    
    #melting eases the plotting process with ggplot
    mean.nbAllelePerSampledDeme.simAvg.melt = melt(mean.nbAllelePerSampledDeme.simAvg,id.vars = "Generation")
    tot.num.alleles.simAvg.melt = melt(tot.num.alleles.simAvg,id.vars = "Generation")
    
    df.Avg.Vs.Total = data.frame(mean.nbAllelePerSampledDeme.simAvg.melt,totValue=tot.num.alleles.simAvg.melt$value)
    df.Avg.Vs.Total.melt = melt(data = df.Avg.Vs.Total,id.vars = c("Generation","variable"),measure.vars = c("value","totValue"),variable.name = "id.tag")
    
    return(
      ggplot(df.Avg.Vs.Total.melt,
             aes(x=Generation,
                 y=value,
                 color=variable,
                 linetype = id.tag))+
        geom_line()
           )
    
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
    compute.Fst = pairwise.fstb(compute.Fst)##TODO CHECK IS CHANGING FROM pairwise.fst TO pairwise.fstb BRINGS ANY ERROR
    
    
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
      #geom_rect(inherit.aes = F,aes(xmin=1-0.4, xmax=1+0.4, ymin=1-0.4, ymax=1+0.4),color="#999999", linetype=3,fill=NA)+
      #geom_rect(inherit.aes = F,aes(xmin=2-0.4, xmax=2+0.4, ymin=2-0.4, ymax=2+0.4),color="#999999", linetype=3,fill=NA)+
      #geom_rect(inherit.aes = F,aes(xmin=3-0.4, xmax=3+0.4, ymin=3-0.4, ymax=3+0.4),color="#999999", linetype=3,fill=NA)+
      #geom_rect(inherit.aes = F,aes(xmin=4-0.4, xmax=5+0.4, ymin=4-0.4, ymax=5+0.4),color="#999999", linetype=3,fill=NA)+
      #geom_rect(inherit.aes = F,aes(xmin=6-0.4, xmax=8+0.4, ymin=6-0.4, ymax=8+0.4), color="#999999", linetype=3,fill=NA)+
      #geom_rect(inherit.aes = F,aes(xmin=9-0.4, xmax=13+0.4, ymin=9-0.4, ymax=13+0.4),color="#999999", linetype=3,fill=NA)+
      coord_fixed()
    #fst.heatmap
    
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
  
    print("Pairwise Fst heatmap plot done.")
    return(fst.heatmap)
  }
  
  ################
  #COMPUTE PAIRWISE_FST/DISTANCE (IBD) PLOT
  ################
  #+++++++++++++++
  Make_pairFst_distance_plots = function(genind_array, generation,between_within_small_seperate_plot = F, vector_small_frags =NULL){
    #genind_array=raw_data
    #generation=6000
    #between_within_small_seperate_plot=T
    #vector_small_frags = c("pop1A", "pop2A")
    
    
    if(between_within_small_seperate_plot == TRUE && is.null(vector_small_frags)){
      stop("between_within_small_seperate_plot is TRUE but vector_small_frags is NULL")
    }
    
    # we create two instances of data because the pairwise.fst function takes genind objects instead of genepop
    raw_data_ind = Get_genind_from_gen(a_genind_array = genind_array,generation = generation)
    raw_data_pop = genind2genpop(raw_data_ind,process.other = T,quiet = T)
    
    
    
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
    
    
    compute.Fst.ibd = pairwise.fst(x = raw_data_ind, res.type = "matrix") # if we use the res.type="dist" we lose line/col names
    compute.Fst.ibd = as.dist(compute.Fst.ibd)
    #is.euclid(compute.Fst.ibd)
    #cailliez(compute.Fst.ibd)
    compute.Fst.ibd = Dist.to.df(compute.Fst.ibd)
    distance.geographic.ibd = Dist.to.df(distance.geographic)
    
    dist.df = data.frame(distance.geographic.ibd,compute.Fst.ibd)
    
    if(between_within_small_seperate_plot == FALSE){
    
    the_plot = ggplot(dist.df)+
      geom_point(aes(x = dist.df$value,y = dist.df$value.1))+
      geom_smooth(aes(x = value, y = value.1), method = 'lm', se = F, fullrange=T)+
      scale_x_continuous(expand=c(0.0,0.0), limits=c(0,max(dist.df$value)),breaks = seq(0,(max(dist.df$value)+1),1)) + #TODO change limits to dynamic values
      scale_y_continuous(expand=c(0.0,0.0), limits=c(0,max(dist.df$value.1)+0.11)) +
      coord_cartesian(xlim = c(0, max(dist.df$value)+1), ylim = c(0, 0.2))+
      ylab("Paiwise Fst")+
      xlab("Geographical distance between demes")+
      ggtitle(paste("Generation ",raw_data_ind@other$generation))#+
      #cowplot::theme_cowplot()
    } else {
      
      within_between.df = dist.df
      
  
        value.tag = NULL
        size.of.df = nrow(within_between.df)
        #within_between.df= within_between.df
        for(i in 1:size.of.df){
  
          if(any(within_between.df[i,2] == vector_small_frags)){
            value.tag = c(value.tag, "small_frag")
          } else if (as.numeric(strsplit(as.character(within_between.df[i,2]),split = "")[[1]][4])==as.numeric(strsplit(as.character(within_between.df[i,1]),split = "")[[1]][4])){
            value.tag = c(value.tag, "within_frag")
          } else {
            value.tag = c(value.tag, "between_frag")
          }
          
        }
  
      
      within_between.df$value.tag =value.tag
      
    
      the_plot = ggplot(within_between.df)+
        geom_point(aes(x = within_between.df$value,y = within_between.df$value.1,
                       color=within_between.df$value.tag,
                       shape=within_between.df$value.tag))+
        #geom_smooth(aes(x = value, y = value.1), method = 'lm', se = F, fullrange=T)+
        geom_smooth(data =subset(x = within_between.df,value.tag=="between_frag") ,aes(x = value, y = value.1, color=value.tag), method = 'lm', se = F, fullrange=T)+
        geom_smooth(data =subset(x = within_between.df,value.tag=="small_frag") ,aes(x = value, y = value.1, color=value.tag), method = 'lm', se = F, fullrange=T)+
        geom_smooth(data =subset(x = within_between.df,value.tag=="within_frag") ,aes(x = value, y = value.1, color=value.tag), method = 'lm', se = F, fullrange=T)+
        scale_x_continuous(expand=c(0.0,0.0), limits=c(0,max(dist.df$value)),breaks = seq(0,(max(dist.df$value)+1),1)) + #TODO change limits to dynamic values
        scale_y_continuous(expand=c(0.0,0.0), limits=c(0,max(dist.df$value.1)+0.11)) +
        coord_cartesian(xlim = c(0, max(dist.df$value)+1), ylim = c(0, 0.2))+
        ylab("Paiwise Fst")+
        xlab("Geographical distance between demes")+
        theme(legend.position=c(0.8, 0.8))+
        ggtitle(paste("Generation ",raw_data_ind@other$generation))
      the_plot
      
      
      if(FALSE){
      #########
      #########
      #########
      # TODO
      genind_array=raw_data
      generation=6000
      between_within_small_seperate_plot=T
      vector_small_frags = c("pop1A", "pop2A")
      
      
      head(within_between.df, n=10)
      
      coloring.pop = lapply(X=strsplit(x=as.character(within_between.df$col.1),split = ""),function(X){
        return(X[4])
      })
  
      within_between.df$the.coloring.pop = as.character(coloring.pop)
      within_between.df.subsetGeo = subset(x = within_between.df,value<14)
      
      the_plot = ggplot(within_between.df,aes(x = value,
                                              y = value.1,
                                              colour=value.tag
                                              ))+
        
        geom_point(aes(shape=the.coloring.pop),size = 3, alpha = 0.8)+
        
        scale_shape_manual(
          values = c(0:9)) +
        
        geom_smooth(data =subset(x = within_between.df,value.tag=="between_frag"),
                    inherit.aes = F,
                    aes(x = value, y = value.1, colour=value.tag),
                    method = 'lm', se = F, fullrange=T)+
        
        geom_smooth(data =subset(x = within_between.df,value.tag=="small_frag"),
                    inherit.aes = F,
                    aes(x = value, y = value.1, colour=value.tag),
                    method = 'lm', se = F, fullrange=T)+
        
        geom_smooth(data =subset(x = within_between.df,value.tag=="within_frag"),
                    inherit.aes = F,
                    aes(x = value, y = value.1, colour=value.tag),
                    method = 'lm', se = F, fullrange=T)+
        
        geom_smooth(data =subset(x = within_between.df.subsetGeo,value.tag=="between_frag"),
                    inherit.aes = F,
                    aes(x = value, y = value.1, colour="between_sub14"),
                    method = 'lm', se = F, fullrange=T)+
        
        scale_x_continuous(expand=c(0.0,0.0), limits=c(0,max(within_between.df$value)),
                           breaks = seq(0,(max(within_between.df$value)+1),1)) + 
        
        scale_y_continuous(expand=c(0.0,0.0), limits=c(0,max(within_between.df$value.1)+0.11)) +
        
        coord_cartesian(xlim = c(0, max(dist.df$value)+1), ylim = c(0, max(dist.df$value.1)+0.05))+
        ylab("Paiwise Fst")+
        xlab("Geographical distance between demes")+
        theme(legend.position=c(0.8, 0.8))+
        ggtitle(paste("Generation ",raw_data_ind@other$generation))
      the_plot
      #########
      #########
      #########
      }
    } 
    
    print("IBD plot done.")
    return(the_plot)
  }
  
  Make_pairFst_distance_plot_QUEMERE_Paper = function(){
    #doi: 10.1111/j.1365-294X.2010.04581.x
    
    df=data.frame("G","F",0.09,288.6,stringsAsFactors = F)
    colnames(x=df) = c("forest1","forest2","pFst","distPx")
    df
    
    df = rbind(df,c("G","H",0.04,206.0))
    df = rbind(df,c("G","B",0.10,224.0))
    df = rbind(df,c("G","A",0.15,322.2))
    df = rbind(df,c("G","C",0.11,259.7))
    df = rbind(df,c("G","D",0.10,210.5))
    df = rbind(df,c("G","E",0.09,119.6))
    #df = rbind(df,c("G","I",0.08,235.5))
    
    df = rbind(df,c("F","H",0.13,319.2))
    df = rbind(df,c("F","B",0.05,432.6))
    df = rbind(df,c("F","A",0.14,481.9))
    df = rbind(df,c("F","C",0.04,304.4))
    df = rbind(df,c("F","D",0.03,170.2))
    df = rbind(df,c("F","E",0.03,173.6))
    #df = rbind(df,c("F","I",0.19,430.8))
    
    df = rbind(df,c("H","B",0.15,427.7))
    df = rbind(df,c("H","A",0.15,528.9))
    df = rbind(df,c("H","C",0.16,444.5))
    df = rbind(df,c("H","D",0.15,346.6))
    df = rbind(df,c("H","E",0.15,232.0))
    #df = rbind(df,c("H","I",0.12,123.0))
    
    df = rbind(df,c("B","A",0.13,111.5))
    df = rbind(df,c("B","C",0.06,194.8))
    df = rbind(df,c("B","D",0.07,275.5))
    df = rbind(df,c("B","E",0.07,279.4))
    #df = rbind(df,c("B","I",0.21,426.6))
    
    df = rbind(df,c("A","C",0.13,190.4))
    df = rbind(df,c("A","D",0.14,313.1))
    df = rbind(df,c("A","E",0.15,352.0))
    #df = rbind(df,c("A","I",0.30,536.0))
    
    df = rbind(df,c("C","D",0.04,136.6))
    df = rbind(df,c("C","E",0.05,220.9))
    #df = rbind(df,c("C","I",0.23,496.5))
    
    df = rbind(df,c("D","E",0.01,116.0))
    #df = rbind(df,c("D","I",0.22,425.8))
    
    #df = rbind(df,c("E","I",0.22,311.3))
    
    df$forest1 = as.factor(df$forest1)
    df$forest2 = as.factor(df$forest2)
    df$pFst = as.numeric(df$pFst)
    df$distPx = as.numeric(df$distPx)
    
    #140px = 10km
    #
    #x=(px*10)/140
    
    df$distKm = as.numeric(lapply(X=df$distPx,function(x){return(newDist = (x*10)/140)}))
    #str(df)
    
    the_model = lm(data = df, formula = df$pFst ~ df$distKm)
    #smooth_values=data.frame(predict(the_model,se = T)) ## builds smooth values along with standard error
    
    library(ggplot2)
    library(plotly)
    p = ggplot(df) + geom_point(aes(x = distKm,y = pFst))+
      geom_label(aes(x = distKm,y = pFst,label=paste(forest1,forest2)))+
      geom_smooth(aes(x = distKm,y = pFst), method = 'lm', se = F, fullrange=T)
    p = p+  scale_x_continuous(expand=c(0.0,0.0), limits=c(0,max(df$distKm)+2),breaks = seq(0,(max(df$distKm)+2),1)) + #TODO change limits to dynamic values
      scale_y_continuous(expand=c(0.0,0.0), limits=c(0,max(df$pFst)+0.11)) +
      coord_cartesian(xlim = c(0, 27), ylim = c(0, 0.2))+
      #geom_abline(slope=the_model$coefficients[[2]],intercept = the_model$coefficients[[1]],color="red")+
      #xlim(0,26)+ylim(0,0.4)+
      ylab("Paiwise Fst")+
      xlab("Geographical distance forests")+
      ggtitle(label = paste("IBD plot - The Daraina Area"),subtitle = "Quéméré et al.")+
      annotate("text",x=7,y=c(0.19,0.18), label=c(paste("slope=",the_model$coefficients[[2]]),paste("intercept=",the_model$coefficients[[1]], sep = "")))
    p
    #ggplotly()
    return(p)
    
    
  }
  
  Make_pairFst_distance_plots_allSims = function(array_of_genind_array,
                                                 name_of_simulation, 
                                                 generation,
                                                 distanceMethod = "manhattan",
                                                 between_within_small_seperate_plot = F, 
                                                 vector_small_frags =NULL){
    
    ##test
    #array_of_genind_array=raw_data_multi_sim
    #str(raw_data_multi_sim,give.head = T,max.level = 2)
    #dim(raw_data_multi_sim)
    #dimnames(raw_data_multi_sim)
    #class(raw_data_multi_sim[1,2][[1]])
    #class(raw_data_single_sim[[1]]@other$generation)
    #raw_data_multi_sim[1,1]
    #raw_data_multi_sim[[1]][[2]]@other$generation
    
    #generation=6500
    #between_within_small_seperate_plot=F
    #vector_small_frags = c("pop1A", "pop2A")
    ##
    
    if(between_within_small_seperate_plot == TRUE && is.null(vector_small_frags)){
      stop("between_within_small_seperate_plot is TRUE but vector_small_frags is NULL")
    }
    
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
    
    
    raw_data_ind = Get_genind_from_gen(a_genind_array = array_of_genind_array[[1]],generation = generation)
    raw_data_pop = genind2genpop(raw_data_ind,process.other = T,quiet = T)
    
    
    distance.geographic = dist(raw_data_pop@other$xy,method = distanceMethod)
    
    Fst.ibd.list=data.frame()
    for(k in 1:length(array_of_genind_array))
    {
    # we create two instances of data because the pairwise.fst function takes genind objects instead of genepop
    raw_data_ind = Get_genind_from_gen(a_genind_array = array_of_genind_array[[k]],generation = generation)
    
  
    compute.Fst.ibd = pairwise.fstb(gsp = raw_data_ind) # if we use the res.type="dist" we lose line/col names
    compute.Fst.ibd = as.dist(compute.Fst.ibd)
    #is.euclid(compute.Fst.ibd)
    #cailliez(compute.Fst.ibd)
    compute.Fst.ibd = Dist.to.df(compute.Fst.ibd)
    #class(compute.Fst.ibd)
    if(empty(Fst.ibd.list))
      {Fst.ibd.list = cbind(as.numeric(compute.Fst.ibd$value))}
    else
      {Fst.ibd.list = cbind(Fst.ibd.list,as.numeric(compute.Fst.ibd$value))}
    
    }
    
    #head(data.frame(Fst.ibd.list))
    Fst.ibd.df = data.frame(Fst.ibd.list)
    Fst.ibd.df$mean = apply(Fst.ibd.df,1,mean,na.rm=T) #does the mean over all the col (simulations) values for each row
    #head(Fst.ibd.df)
    
    distance.geographic.ibd = Dist.to.df(distance.geographic)
    dist.df = data.frame(distance.geographic.ibd, Fst.ibd.df)
    #head(dist.df)
    
    if(between_within_small_seperate_plot == FALSE){
      
      the_plot = ggplot(dist.df)+
        geom_point(aes(x = dist.df$value,y = dist.df$mean))+
        geom_smooth(aes(x = value, y = mean), method = 'lm', se = F, fullrange=T)+
        scale_x_continuous(expand=c(0.0,0.0), limits=c(0,max(dist.df$value)),breaks = seq(0,(max(dist.df$value)+1),1)) + #TODO change limits to dynamic values
        scale_y_continuous(expand=c(0.0,0.0), limits=c(0,max(dist.df$mean)+0.1)) +
        coord_cartesian(xlim = c(0, max(dist.df$value)+1), ylim = NULL)+
        ylab("Paiwise Fst")+
        xlab("Geographical distance between demes")+
        ggtitle(paste("Generation ",raw_data_ind@other$generation),subtitle = name_of_simulation)#+
      #cowplot::theme_cowplot()
    } else {
      
      within_between.df = dist.df
      
      
      value.tag = NULL
      size.of.df = nrow(within_between.df)
      #within_between.df= within_between.df
      for(i in 1:size.of.df){
        
        if(any(within_between.df[i,2] == vector_small_frags)){
          value.tag = c(value.tag, "small_frag")
        } else if (as.numeric(strsplit(as.character(within_between.df[i,2]),split = "")[[1]][4])==as.numeric(strsplit(as.character(within_between.df[i,1]),split = "")[[1]][4])){
          value.tag = c(value.tag, "within_frag")
        } else {
          value.tag = c(value.tag, "between_frag")
        }
        
      }
      
      
      within_between.df$value.tag =value.tag
      #head(within_between.df, n = 100)
      
      the_plot = ggplot(within_between.df)+
        geom_point(aes(x = within_between.df$value,y = within_between.df$mean,
                       color=within_between.df$value.tag,
                       shape=within_between.df$value.tag))+
        #geom_smooth(aes(x = value, y = value.1), method = 'lm', se = F, fullrange=T)+
        geom_smooth(data =subset(x = within_between.df,value.tag=="between_frag") ,aes(x = value, y = mean, color=value.tag), method = 'lm', se = F, fullrange=T)+
        geom_smooth(data =subset(x = within_between.df,value.tag=="small_frag") ,aes(x = value, y = mean, color=value.tag), method = 'lm', se = F, fullrange=T)+
        geom_smooth(data =subset(x = within_between.df,value.tag=="within_frag") ,aes(x = value, y = mean, color=value.tag), method = 'lm', se = F, fullrange=T)+
        scale_x_continuous(expand=c(0.0,0.0), limits=c(0,max(dist.df$value)),breaks = seq(0,(max(dist.df$value)+1),1)) + #TODO change limits to dynamic values
        scale_y_continuous(expand=c(0.0,0.0), limits=c(0,max(dist.df$mean)+0.11)) +
        coord_cartesian(xlim = c(0, max(dist.df$value)+1), ylim = c(0, 0.2))+
        ylab("Paiwise Fst")+
        xlab("Geographical distance between demes")+
        theme(legend.position=c(0.8, 0.8))+
        ggtitle(paste("Generation ",raw_data_ind@other$generation),subtitle = name_of_simulation)
      the_plot
      
    } 
    
    print("IBD plot done.")
    return(the_plot)
  }
  
  pairFst_dist_plots_allSims_fourFragClasses = function(array_of_genind_array,
                                                        generation,
                                                        name_of_simulation,
                                                        distanceMethod = "manhattan",
                                                        vector_small_frags = c(2,4,5,7,14,15,16,17,23,26,27,28,33,34,35,41,46,47,48,49),
                                                        vector_mediumSmall_frags = c(1,3,6,8,12,13,18,19,21,22,29,30,31,32,37,39,42,43,44,45),
                                                        vector_mediumLarge_frags= c(9,10,11,24,25,36,38,40),
                                                        vector_large_frags = c(20)){
    
    ##test
    #array_of_genind_array=raw_data_multi_sim[1:2]
    #str(array_of_genind_array,give.head = T,max.level = 2)
    #dim(raw_data_multi_sim)
    #dimnames(raw_data_multi_sim)
    #class(raw_data_multi_sim[1,2][[1]])
    #class(raw_data_single_sim[[1]]@other$generation)
    #raw_data_multi_sim[1,1]
    #raw_data_multi_sim[[1]][[2]]@other$generation
    
    #generation=6000
    between_within_small_seperate_plot=T
    #vector_small_frags = c(2,4,5,7,14,15,16,17,23,26,27,28,33,34,35,41,46,47,48,49)
    #vector_mediumSmall_frags = c(1,3,6,8,12,13,18,19,21,22,29,30,31,32,37,39,42,43,44,45)
    #vector_mediumLarge_frags= c(9,10,11,24,25,36,38,40)
    #vector_large_frags = c(20)
    ##
    
    if(between_within_small_seperate_plot == TRUE && is.null(vector_small_frags)){
      stop("between_within_small_seperate_plot is TRUE but vector_small_frags is NULL")
    }
    

    
    raw_data_ind = Get_genind_from_gen(a_genind_array = array_of_genind_array[[1]],generation = generation)
    raw_data_pop = genind2genpop(raw_data_ind,process.other = T,quiet = T)

    #raw_data_ind@other$xy
    distance.geographic = dist(raw_data_pop@other$xy,method = distanceMethod)
    
    Fst.ibd.list=data.frame()
    
    Dgen.list = NULL
    
    raw_data_pop = NULL
    
    pairwiseFstMatrixList = vector(mode = "list", length = length(array_of_genind_array))
    
    for(k in 1:length(array_of_genind_array))
    {
      # we create two instances of data because the pairwise.fst function takes genind objects instead of genepop
      raw_data_ind = Get_genind_from_gen(a_genind_array = array_of_genind_array[[k]],generation = generation)
      
      if(FALSE){#used to calc mantelrandtest over all sims
        raw_data_pop = genind2genpop(raw_data_ind,process.other = T,quiet = T)
        compute.gen.dist = dist(raw_data_pop)
        Dgen.list[[k]] = (compute.gen.dist)
      }
      
      compute.Fst.ibd_newMethod = pairwise.fstb(gsp = raw_data_ind)
      compute.Fst.ibd_newMethod = as.dist(compute.Fst.ibd_newMethod,diag = F,upper = F)
      pairwiseFstMatrixList[[k]] = (compute.Fst.ibd_newMethod)
      
      ####compute.Fst.ibd = pairwise.fst(x = raw_data_ind, res.type = "matrix") # if we use the res.type="dist" we lose line/col names
      ####compute.Fst.ibd = as.dist(compute.Fst.ibd)
      #is.euclid(compute.Fst.ibd)
      #cailliez(compute.Fst.ibd)
      ####compute.Fst.ibd = Dist.to.df(compute.Fst.ibd)
      #class(compute.Fst.ibd)
      
      if(FALSE){
      if(empty(Fst.ibd.list)){
        Fst.ibd.list = cbind(as.numeric(compute.Fst.ibd$value))
      }else{
        Fst.ibd.list = cbind(Fst.ibd.list,as.numeric(compute.Fst.ibd$value))
      }
      }
  }

    if(FALSE){
    # Reduces all the lists to 1 by adding their individual elements
    # we then divide that value with the length of the original list to get the mean
    # could possible use this method to some of the other objects where im also calculating over lists of matrices
    Dgen.reduce = Reduce('+', Dgen.list)
    Dgen.mean = Dgen.reduce/length(Dgen.list)
    
    
    ## Calc mantel randtest
    ibd_mantel_randtest = mantel.randtest(Dgen.mean, distance.geographic)
    }
    
    
    # head(data.frame(Fst.ibd.list))
    # Mean over each of the pairwise values for all the simulations
    pairwiseFstReduce = Reduce('+',pairwiseFstMatrixList)
    pairwiseFstMean = pairwiseFstReduce/length(pairwiseFstMatrixList)
    # head(pairwiseFstMean)
    pairwiseFstMean_df = Dist.to.df(pairwiseFstMean)
    
    #Fst.ibd.df = data.frame(Fst.ibd.list)
    #Fst.ibd.df$mean = apply(Fst.ibd.df,1,mean,na.rm=T) #does the mean over all the col (simulations) values for each row
    #head(Fst.ibd.df)
    
    distance.geographic.ibd = Dist.to.df(distance.geographic)
    #dist.df = data.frame(distance.geographic.ibd,Fst.ibd.df)
    
    #geographic and genetic data into one dataframe. Mean values of pairwiseFst renamed to "mean"
    dist.df = data.frame(distance.geographic.ibd,mean = pairwiseFstMean_df$value)
    
   
    
      
      within_between.df = dist.df
      #within_between.df[,2]
      #head(within_between.df)
      
      # the following line separates the subpops into pops, for example pop13A will be "13" and pop1B will be "1"
      #strsplit(as.character(within_between.df[1432,2]),split = "[aA-zZ]+")[[1]][2]
      #strsplit(as.character(within_between.df[1432,1]),split = "[aA-zZ]+")[[1]][2]
      
      
      value.tag.between = NULL
      value.tag.within = NULL
      value.tag.withinBetween = NULL
      
      size.of.df = nrow(within_between.df)
      #within_between.df= within_between.df
      
      
      
      for(i in 1:size.of.df){
        
        # between the fragment contained in the vector and everyone else
        if (as.numeric(strsplit(as.character(within_between.df[i,2]),split = "[aA-zZ]+")[[1]][2])==as.numeric(strsplit(as.character(within_between.df[i,1]),split = "[aA-zZ]+")[[1]][2])){
          value.tag.between = c(value.tag.between, "within_frag")
        }else if(any(as.numeric(strsplit(as.character(within_between.df[i,2]),split = "[aA-zZ]+")[[1]][2]) == as.numeric(vector_mediumSmall_frags))){
          value.tag.between = c(value.tag.between, "mediumSmall_frag")
        }else if(any(as.numeric(strsplit(as.character(within_between.df[i,2]),split = "[aA-zZ]+")[[1]][2]) == as.numeric(vector_mediumLarge_frags))){
          value.tag.between = c(value.tag.between, "mediumLarge_frag")
        }else if(any(as.numeric(strsplit(as.character(within_between.df[i,2]),split = "[aA-zZ]+")[[1]][2]) == as.numeric(vector_large_frags))){
          value.tag.between = c(value.tag.between, "large_frag")
        }else if(any(as.numeric(strsplit(as.character(within_between.df[i,2]),split = "[aA-zZ]+")[[1]][2]) == as.numeric(vector_small_frags))){
          value.tag.between = c(value.tag.between, "small_frag")
        }
        
        # between the sampled demes belonging to the same fragment (within fragment)
        if (as.numeric(strsplit(as.character(within_between.df[i,2]),split = "[aA-zZ]+")[[1]][2])==as.numeric(strsplit(as.character(within_between.df[i,1]),split = "[aA-zZ]+")[[1]][2])){
          value.tag.within = c(value.tag.within, "within_frag")
        } else{
          value.tag.within = c(value.tag.within, "between_frag")
        }
        
        # between the small and everyone else
        if(any(as.numeric(strsplit(as.character(within_between.df[i,2]),split = "[aA-zZ]+")[[1]][2]) == as.numeric(vector_small_frags))){
          value.tag.withinBetween = c(value.tag.withinBetween, "small_frag")
        }
        else if (as.numeric(strsplit(as.character(within_between.df[i,2]),split = "[aA-zZ]+")[[1]][2])==as.numeric(strsplit(as.character(within_between.df[i,1]),split = "[aA-zZ]+")[[1]][2])){
          value.tag.withinBetween = c(value.tag.withinBetween, "within_frag")
        } 
        else {
          value.tag.withinBetween = c(value.tag.withinBetween, "between_frag")
        }
        
      }
      
      
      within_between.df$value.tag.between = value.tag.between
      within_between.df$value.tag.within = value.tag.within
      within_between.df$value.tag.withinBetween = value.tag.withinBetween
      #tail(within_between.df, n = 999)
      
      the_plot_all = ggplot(within_between.df)+
        geom_point(aes(x = within_between.df$value,y = within_between.df$mean,
                       color=within_between.df$value.tag.between,
                       shape=within_between.df$value.tag.between))+
        geom_smooth(color="black",aes(x = value, y = mean, color="Overall"), method = 'lm', se = F, fullrange=T)
      
      
      the_plot_between = ggplot() + 
        stat_smooth(data =subset(x = within_between.df,value.tag.between=="within_frag") ,aes(x = value, y = mean, color=value.tag.between), method = 'lm', se = F, fullrange=T)+
        geom_smooth(data =subset(x = within_between.df,value.tag.between=="mediumSmall_frag") ,aes(x = value, y = mean, color=value.tag.between), method = 'lm', se = F, fullrange=T)+
        geom_smooth(data =subset(x = within_between.df,value.tag.between=="mediumLarge_frag") ,aes(x = value, y = mean, color=value.tag.between), method = 'lm', se = F, fullrange=T)+
        geom_smooth(data =subset(x = within_between.df,value.tag.between=="large_frag") ,aes(x = value, y = mean, color=value.tag.between), method = 'lm', se = F, fullrange=T)+
        geom_smooth(data =subset(x = within_between.df,value.tag.between=="small_frag") ,aes(x = value, y = mean, color=value.tag.between), method = 'lm', se = F, fullrange=T)
      
      
      if(FALSE){
      new_data_test = subset(x = within_between.df,value.tag.between=="within_frag")## subset data
      new_model = lm(data = new_data_test,
                     formula = new_data_test$mean ~ new_data_test$value) ## create linear regression model gen ~ geo
      
      
      
        ## method to build geom_smooth of the data without actually using geom_smooth
        
        smooth_values=data.frame(predict(new_model,se = T)) ## builds smooth values along with standard error
      
      
        
        df <- data.frame(cbind(value = new_data_test$value,
                               mean = new_data_test$mean,
                               fit = smooth_values$fit,
                               upperBound = smooth_values$fit + 2 * smooth_values$se.fit,
                               lowerBound = smooth_values$fit - 2 * smooth_values$se.fit
        ))
        
        g <- ggplot(df, aes(value, mean))
        #g <- g + geom_point()
        g <- g + geom_linerange(aes(ymin = lowerBound, ymax = upperBound))
        g <- g + geom_point(aes(value, fit))
        g <- g + geom_smooth(method = "lm",color="green", size=2)
        g <- g + geom_abline(slope = new_model$coefficients[[2]],intercept = new_model$coefficients[[1]], color="red")
        g
      }
      
      the_plot_within = ggplot() + 
        geom_smooth(data =subset(x = within_between.df,value.tag.within=="within_frag") ,aes(x = value, y = mean, color=value.tag.within), method = 'lm', se = F, fullrange=T)+
        geom_smooth(data =subset(x = within_between.df,value.tag.within=="between_frag") ,aes(x = value, y = mean, color=value.tag.within), method = 'lm', se = F, fullrange=T)
      
      the_plot_withinBetween = ggplot() + 
        geom_smooth(data =subset(x = within_between.df,value.tag.withinBetween=="within_frag") ,aes(x = value, y = mean, color=value.tag.withinBetween), method = 'lm', se = F, fullrange=T)+
        geom_smooth(data =subset(x = within_between.df,value.tag.withinBetween=="between_frag") ,aes(x = value, y = mean, color=value.tag.withinBetween), method = 'lm', se = F, fullrange=T)+
        geom_smooth(data =subset(x = within_between.df,value.tag.withinBetween=="small_frag") ,aes(x = value, y = mean, color=value.tag.withinBetween), method = 'lm', se = F, fullrange=T)
      
      the_plots_list = list(the_plot_all,the_plot_between,the_plot_within,the_plot_withinBetween)
      
      for(i in 1:4){
        the_plots_list[[i]]=
          the_plots_list[[i]]+
          scale_x_continuous(expand=c(0.0,0.0), limits=c(0,max(dist.df$value)+1),breaks = seq(0,(max(dist.df$value)+1),1)) + 
          scale_y_continuous(expand=c(0.0,0.0), limits=c(0,max(dist.df$mean)+0.11)) +
          coord_cartesian(xlim = c(0, max(dist.df$value)+1), ylim = c(0, 1.0))+
          ylab("Paiwise Fst")+
          xlab("Geographical distance between demes")+
          theme(legend.position=c(0.8, 0.8))+
          ggtitle(paste("Generation ",raw_data_ind@other$generation),subtitle = name_of_simulation)
          
      }
      
     
    
    print("IBD plot done.")
    
    return(the_plots_list)
  }
  
  
  ################
  #COMPUTE PAIRWISE_FST/DISTANCE (IBD) PLOT SLOPE OVER TIME
  ################
  #+++++++++++++++ 
  Compute.PairwiseFst.IBD.Slope.OverTime = function(array_of_genind_array,name_of_simulation, distanceMethod = "manhattan",isPairwiseFstdist = TRUE){
    
    

    #array_of_genind_array=raw_data_multi_sim
    #distanceMethod = "manhattan"
    #isPairwiseFstdist = TRUE
    #*******
    
    generation_list=NULL
    Dgen.mean = NULL
    #length(levels(array_of_genind_array[[1]][[50]]@pop))
    
    for(j in 1:length(array_of_genind_array[[1]])){ # 1:length(array_of_genind_array[[1]])
      Dgen.list = NULL
    for(i in 1:length(array_of_genind_array)){
      
      #drop POPs that have less than 2 individuals (related to problems with pairwise.fstb)
      data_with_dropPop = selPopSize(array_of_genind_array[[i]][[j]],nMin=2)
      
      # choose if we want dist vs pairwiseFst method
      if(isPairwiseFstdist){
        compute.gen.dist = pairwise.fstb(gsp = data_with_dropPop)
        compute.gen.dist = as.dist(compute.gen.dist ,diag = F,upper = F)
      }else{
        raw_data_pop = genind2genpop(data_with_dropPop,process.other = T,quiet = T)
        compute.gen.dist = dist(raw_data_pop)
      }
      Dgen.list[[i]] = (compute.gen.dist)
      
      #str(Dgen.list)
      #apply(simplify2array(Dgen.list), c(1,2),mean)
      
    }
      
      Dgen.reduce = Reduce('+', Dgen.list)
      Dgen.mean[[j]] = Dgen.reduce/length(Dgen.list)
      generation_list = c(generation_list,as.numeric(array_of_genind_array[[1]][[j]]@other$generation))
      
    }
    
    raw_data_pop = genind2genpop(array_of_genind_array[[1]][[1]],process.other = T,quiet = T)
    distance.geographic = dist(raw_data_pop@other$xy,method = distanceMethod)
    distance.geographic = Dist.to.df(distance.geographic)
    
    Dgen.to.df.allSims=NULL
    
    for(i in 1:length(Dgen.mean)){
      Dgen.mean[[i]] = as.dist(Dgen.mean[[i]])
      if(empty(Dgen.to.df.allSims)){
        Dgen.mean.df = Dist.to.df(Dgen.mean[[i]])
      Dgen.to.df.allSims = data.frame(row=Dgen.mean.df$row,col=Dgen.mean.df$col,generationMean=Dgen.mean.df$value)
    }else{
      Dgen.mean.df = Dist.to.df(Dgen.mean[[i]])
      Dgen.to.df.allSims = data.frame(Dgen.to.df.allSims, generationMean=Dgen.mean.df$value)
      }
    }
    
    
    new_data_test = data.frame(distance.geographic)
    new_data_test = data.frame(new_data_test,Dgen.to.df.allSims[,3:length(Dgen.to.df.allSims)])
    
    
    #plot data pFst vs geodist as in normal IBD plot
    if(FALSE){
    #just testing stuff
      ggplot(new_data_test)+
      geom_smooth(aes(x=new_data_test$value,y=new_data_test$generationMean.19, col="19"),method = "lm",se = T)+
      geom_point(aes(x=new_data_test$value,y=new_data_test$generationMean.19, col="19"))+
      geom_smooth(aes(x=new_data_test$value,y=new_data_test$generationMean, col="1",size=2),method = "lm",se = T)+
      geom_point(aes(x=new_data_test$value,y=new_data_test$generationMean, col="1"))+
      geom_abline(slope = 0.0006182,intercept = 0.0351354, color = "black")
    }
    
    
    slope.df = NULL
    intercept.df = NULL
    all_model_information = NULL
    # excluding first 4 cols because they correspond to "row", "col", "value" and the first generation
    # first gen full of NaNs
    for(i in 5:(length(new_data_test))){
    new_model = lm(data = new_data_test,
                   formula = new_data_test[,i] ~ new_data_test$value) ## create linear regression model gen ~ geo
    all_model_information[[i-3]] = new_model
    slope.df[i-4] = new_model$coefficients[[2]]
    intercept.df[i-4] = new_model$coefficients[[1]]
    }
    
    
    #head(slope.df[1:20,])
    slope.df = data.frame(mean = slope.df,generation_list=generation_list[2:length(generation_list)])
    intercept.df= data.frame(mean = intercept.df,generation_list=generation_list[2:length(generation_list)])
    #head(slope.df)
    #length(all_model_information)
    
    p.slope = ggplot(slope.df)+
      #geom_line(aes(x=slope.df$generation_list, y=slope.df$mean))+
      geom_point(aes(x=slope.df$generation_list, y=slope.df$mean))+
      #geom_abline(data = all_model_information[[2]],slope = all_model_information[[2]]$coefficients[[2]],intercept = all_model_information[[2]]$coefficients[[1]], color="red")+
      #geom_abline(data = all_model_information[[1201]],slope = all_model_information[[1201]]$coefficients[[2]],intercept = all_model_information[[1201]]$coefficients[[1]], color="blue")+
      geom_smooth(aes(x=slope.df$generation_list, y=slope.df$mean))+
      coord_cartesian(xlim = NULL, ylim = NULL, expand = T)+ggtitle(label = "IBD slope over time",subtitle = name_of_simulation)
    
    p.intercept = ggplot(intercept.df)+
      #geom_line(aes(x=intercept.df$generation_list, y=intercept.df$mean))+
      geom_point(aes(x=intercept.df$generation_list, y=intercept.df$mean))+
      #geom_abline(data = all_model_information[[2]],slope = all_model_information[[2]]$coefficients[[2]],intercept = all_model_information[[2]]$coefficients[[1]], color="red")+
      #geom_abline(data = all_model_information[[1201]],slope = all_model_information[[1201]]$coefficients[[2]],intercept = all_model_information[[1201]]$coefficients[[1]], color="blue")+
      geom_smooth(aes(x=intercept.df$generation_list, y=intercept.df$mean))+
      coord_cartesian(xlim = NULL, ylim = NULL, expand = T)+ggtitle(label = "IBD intercept over time",subtitle = name_of_simulation)
    
    
    
    if(FALSE){
    ggplot()+geom_abline(data = all_model_information[[2]],slope = all_model_information[[2]]$coefficients[[2]],intercept = all_model_information[[2]]$coefficients[[1]], color="red")+
      geom_abline(data = all_model_information[[1197]],slope = all_model_information[[1197]]$coefficients[[2]],intercept = all_model_information[[1197]]$coefficients[[1]], color="blue")+
      coord_cartesian(xlim = c(-1,1), ylim = c(-1,1), expand = F)
    }
      
    if(FALSE){
    ## method to build geom_smooth of the data without actually using geom_smooth
    smooth_values=data.frame(predict(new_model,se = T)) ## builds smooth values along with standard error
    
    
    
    df <- data.frame(cbind(value = new_data_test$value,
                           mean = new_data_test[[53]],
                           fit = smooth_values$fit,
                           upperBound = smooth_values$fit + 2 * smooth_values$se.fit,
                           lowerBound = smooth_values$fit - 2 * smooth_values$se.fit
    ))
    
    g <- ggplot(df, aes(value, mean))
    #g <- g + geom_point()
    g <- g + geom_linerange(aes(ymin = lowerBound, ymax = upperBound))
    g <- g + geom_point(aes(value, fit))
    g <- g + geom_smooth(method = "lm",color="green", size=2)
    g <- g + geom_abline(slope = new_model$coefficients[[2]],intercept = new_model$coefficients[[1]], color="red")
    g
    }
    
    return(list(p.slope, p.intercept))
    
  }
  
  
  
  
  
  ################
  #TEST AND PLOT IBD HYPOTHESIS
  ################
  #+++++++++++++++
  Test_IBD_hypothesis = function(genind_object){
  
    # The original value of the correlation between the distance matrices is represented by the dot,
    # while histograms represent permuted values (i.e., under the absence of spatial structure).
    # Significant spatial structure would therefore result in the original value being out of the
    # reference distribution.
    Dgen = dist(genind_object)
    Dgeo = dist(genind_object$other$xy)
    ibd = mantel.randtest(Dgen,Dgeo)

  
    Plot_IBD_hypothesis = function(mantelRandTestResult){
      plot(mantelRandTestResult)
    }
    Plot_IBD_hypothesis(ibd)
    return(ibd)
  }
  
  Test_IBD_hypothesis_ggplot = function(genind_object,distanceMethod = "manhattan"){
    
    # The original value of the correlation between the distance matrices is represented by the dot,
    # while histograms represent permuted values (i.e., under the absence of spatial structure).
    # Significant spatial structure would therefore result in the original value being out of the
    # reference distribution.
    
    genind_object = raw_data_multi_sim[[1]][[10]]
    genind_to_pop = genind2genpop(genind_object,process.other = T,quiet = T)
    genind_object@other$xy
    Dgen = dist(genind_to_pop)
    Dgeo = dist(genind_to_pop$other$xy,method = distanceMethod)
    df.Dgen = Dist.to.df(Dgen)
    df.Dgeo = Dist.to.df(Dgeo)
    ibd = mantel.randtest(Dgen,Dgeo)
    as.matrix(as.data.frame(lapply(df.Dgeo, as.numeric)))
    mantel.randtest(df.Dgen,df.Dgeo)
    ibd
    
    df.ibd = data.frame(ibd$sim)
    head(df.ibd)
    ggplot(data=df.ibd) + geom_histogram(aes(x = ibd.sim),
                                         boundary=0,
                                         binwidth = 0.005,
                                         colour="black", fill="white")+
      geom_density(aes(x = ibd.sim),alpha=.2, fill="#FF6666")+
      geom_vline(xintercept = ibd$obs, color="red")
    
    Plot_IBD_hypothesis = function(mantelRandTestResult){
      plot(mantelRandTestResult)
    }
    Plot_IBD_hypothesis(ibd)
    return(ibd)
  }
  
  ###########################################################################
  ###########################################################################
  #++++SECTION 3 - HELPER FUNCTIONS 
  ###########################################################################
  ###########################################################################
  
  summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                        conf.interval=.95, .drop=TRUE) {
    ## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
    ##   data: a data frame.
    ##   measurevar: the name of a column that contains the variable to be summariezed
    ##   groupvars: a vector containing names of columns that contain grouping variables
    ##   na.rm: a boolean that indicates whether to ignore NA's
    ##   conf.interval: the percent range of the confidence interval (default is 95%)
    library(plyr)
    
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
      if (na.rm) sum(!is.na(x))
      else       length(x)
    }
    
    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun = function(xx, col) {
                     c(N    = length2(xx[[col]], na.rm=na.rm),
                       mean = mean   (xx[[col]], na.rm=na.rm),
                       sd   = sd     (xx[[col]], na.rm=na.rm)
                     )
                   },
                   measurevar
    )
    
    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))
    
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    
    return(datac)
  }
  
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
  
  
  ## save plots with set size
  save_the_plot = function(plotObject, pathToPlot = "~/Documents/",nameOfPlotFile){
    ggplot2::ggsave(plot = plotObject, file = paste(pathToPlot,"/", nameOfPlotFile,".pdf",sep=""), units = "cm", width = 18, height = 13)
    
  }
  
  
  
  ###########################################################################
  ###########################################################################
  #++++SECTION 4 - DO DATA ANALYSIS 
  ###########################################################################
  ###########################################################################
  
  
  # MOVED TO "Run_SINS_data_analysis.r"
  
  
  
