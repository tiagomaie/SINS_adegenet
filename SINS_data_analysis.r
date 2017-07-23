
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Functions in this file are used to get data from SINS into R so that
# we can run analysis on this data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++






# Use code below only if development version of adegenet is necessary
if(!require("devtools")){
  install.packages("devtools", dependencies = TRUE)
  install_github("thibautjombart/adegenet")
  library("adegenet")
  }


if (suppressPackageStartupMessages(!require("adegenet")))
  {install.packages("adegenet", dependencies = TRUE)}
if (suppressPackageStartupMessages(!require("ade4")))
  {install.packages("ade4", dependencies = TRUE)}
if (suppressPackageStartupMessages(!require("hierfstat")))
  {install.packages("hierfstat", dependencies = TRUE)}
if (suppressPackageStartupMessages(!require("pegas")))
  {install.packages("pegas", dependencies = TRUE)}
if (suppressPackageStartupMessages(!require("ggplot2")))
  {install.packages("ggplot2", dependencies = TRUE)}
if (suppressPackageStartupMessages(!require("reshape2")))
  {install.packages("reshape2", dependencies = TRUE)}
if (suppressPackageStartupMessages(!require("gridExtra")))
  {install.packages("gridExtra", dependencies = TRUE)}
#if(!require("cowplot")) install.packages("cowplot", dependencies = TRUE)
#if(!require("plotly")) install.packages("plotly", dependencies = TRUE)
if (suppressPackageStartupMessages(!require("Hmisc")))
  {install.packages("Hmisc", dependencies = TRUE)}
if (suppressPackageStartupMessages(!require("plyr")))
  {install.packages("plyr", dependencies = TRUE)}
if (suppressPackageStartupMessages(!require("PopGenReport")))
  {install.packages("PopGenReport", dependencies = TRUE)}
if (suppressPackageStartupMessages(!require("data.table")))
{install.packages("data.table", dependencies = TRUE)}
if (suppressPackageStartupMessages(!require("StAMPP")))
{install.packages("StAMPP", dependencies = TRUE)}
if (suppressPackageStartupMessages(!require("gridExtra")))
{install.packages("gridExtra", dependencies = TRUE)}
if (suppressPackageStartupMessages(!require("doParallel")))
{install.packages("doParallel", dependencies = TRUE)}
if (suppressPackageStartupMessages(!require("foreach")))
{install.packages("foreach", dependencies = TRUE)}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##### SECTION 1 - FUNCTIONS TO GET DATA ###################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++
### SET NUMBER OF SIMULATIONS ################
#+++++++++++++++
Set_number_of_simulations = function(number_of_simulations) {
  return(number_of_simulations)
}


#+++++++++++++++
# SET MARKERS ####
#+++++++++++++++
#+++++++++++++++
Set_markers_list = function(...) {
  # Make list of the given markers (e.g. "A1","A2","A3", ...)
  the_markers_list = NULL
  the_markers_list = c(...)
  return(the_markers_list)
}

################
#SET GENERATIONS
################
#+++++++++++++++
Set_generations_list = function(initialTimeStep,
                                finalTimeStep,
                                timeStepIncrement) {
  # Make sequence of all the generations that we want to study (e.g. from gen 1 to 100 by increments of 5)
  the_generations_list = NULL
  the_generations_list = seq(initialTimeStep, finalTimeStep, timeStepIncrement)
  return(the_generations_list)
}

################
#SET LAYERS
################
#+++++++++++++++
Set_layers_list = function(...) {
  # Make list of the given layers (e.g. "layer0", "layer1")
  the_layers_list = NULL
  the_layers_list = c(...)
  return(the_layers_list)
}

################
#GET FILE NAMES
################
#+++++++++++++++
Get_file_names = function(generations_list, layers_list) {
  # Build the names of the files that we want to analyse
  the_file_names = NULL

  for (j in generations_list) {
    the_file_names = rbind(the_file_names, paste(j, "_", layers_list, "_", sep = ""))

  }
  return(the_file_names)
}

################
#READ DATA
################
#+++++++++++++++
Read_data_func = function(name_of_file,
                          markers_list,
                          is_marker_diploid = TRUE,
                          simulation_ID = 0) {
  #name_of_file = the_file_names[[9]]
  #markers_list = c("A1","A2","A3","A4","A5","A6","A7","A8","A9","A10")
  #markers_list = c("Y")
  #is_marker_diploid=F
  #simulation_ID=1
  #setwd(dir = paste(the_path_to_data,"Adegenet_sim_1","/",sep = ""))

  if (FALSE) {
    #***fread is much faster than an optimized read.table***
    #read first 5 rows of the table
    #tab5rows <- read.table(file = paste(name_of_file,markers_list[1],sep = ""), header = FALSE, nrows = 100,row.names = 1,na.strings = "NA")

    #get the classes of the columns of the table
    #tabClasses <- sapply(tab5rows, class)
  }

  data_file <-
    fread(
      input = paste(name_of_file, markers_list[1], sep = ""),
      encoding = "UTF-8",
      header = FALSE,
      na.strings = "NA",
      data.table = F
    )

  # get the number of rows that the table will have, use it later to speed up process
  # not sure if it is that much faster, consider erasing it
  sizeOfTable = nrow(data_file)

  markerCondition = markers_list == "MT" | markers_list == "Y"

  current_generation = NULL
  current_generation = strsplit(x = name_of_file, split = "_")[[1]][1]

  data_file = NULL
  marker_data.haplo = NULL
  marker_data.diplo = NULL

  a_file_name = NULL

  #bind_markers.diplo = NULL
  bind_markers.diplo = matrix(NA, nrow = sizeOfTable, ncol = length(markers_list))
  #bind_markers.haplo = NULL
  bind_markers.haplo = matrix(NA, nrow = sizeOfTable, ncol = length(markers_list))

  for (marker in 1:length(markers_list)) {
    a_file_name = paste(name_of_file, markers_list[marker], sep = "")

    #read the rest (or a bigger part) of the table knowing the class of the columns speeds up the reading process
    #data_file <- read.table(a_file_name, header = FALSE, colClasses = tabClasses, nrows = sizeOfTable, row.names = 1,na.strings = "NA")
    #fread is much faster than an optimized read.table
    data_file <-
      fread(
        input = a_file_name,
        header = FALSE,
        encoding = "UTF-8",
        na.strings = "NA",
        data.table = F
      )

    #if marker is Y or MTdna then it is haploid data, save that markers col (instead of 2 cols for diploid)
    if (markerCondition[marker]) {
      marker_data.haplo = as.matrix(data_file$V2)
      #bind_markers.haplo = cbind(bind_markers.haplo,marker_data.haplo)
      bind_markers.haplo[, marker] = marker_data.haplo
    } else{
      #marker_data.diplo = as.matrix(paste(data_file$V2,data_file$V3,sep = "/"))

      marker_data.diplo = as.vector(data_file$V2)


      #length(marker_data.diplo)
      #bind_markers.diplo = cbind(bind_markers.diplo,marker_data.diplo)
      bind_markers.diplo[, marker] = marker_data.diplo
    }
  }



  if (is_marker_diploid) {
    # note ncode set to 3 - alleles are considered to be coded by three characters
    markers.diplo_df2genind = df2genind(
      X = bind_markers.diplo,
      ploidy = 2,
      ncode = 3,
      ind.names = data_file$V1,
      sep = "/"
    )
    # define pops
    markers.diplo_df2genind@pop = factor(data_file$V5)
    # define coords
    xy_coords = cbind(data_file$V3, data_file$V4)
    markers.diplo_df2genind@other$xy = xy_coords
    # define gen
    markers.diplo_df2genind@other$generation = current_generation
    # define simulation ID
    markers.diplo_df2genind@other$simulation.ID = simulation_ID
    # return genind object
    return(markers.diplo_df2genind)
  } else{
    markers.haplo_df2genind = df2genind(
      X = bind_markers.haplo,
      ploidy = 1,
      ncode = 3,
      ind.names = data_file$V1,
      sep = "/",
      NA.char = "NA",
      pop=factor(data_file$V5)
    )
    # define pops
    markers.haplo_df2genind@pop = factor(data_file$V5)
    # define coords
    xy_coords = cbind(data_file$V3, data_file$V4)
    markers.haplo_df2genind@other$xy = xy_coords
    # define gen
    markers.haplo_df2genind@other$generation = current_generation
    # define simulation ID
    markers.haplo_df2genind@other$simulation.ID = simulation_ID
    # return genind object
    return(markers.haplo_df2genind)
  }

}


### READ SINS DATA ####
#much improved version of the read_data_func+multi_sim functions
#to get data from sins to adegenet
#TODO: choose which generations to get
#TODO: enable option to load haploid data
readSinsData = function(dataDirectory, simulationName, nSim=c(1:2), markers = c(...),layerName,dataFormat,n_cores=2){

  if (suppressPackageStartupMessages(!require("doParallel")))
  {install.packages("doParallel", dependencies = TRUE)}
  if (suppressPackageStartupMessages(!require("foreach")))
  {install.packages("foreach", dependencies = TRUE)}

  #dataDirectory = "~/ServerFolders/Haystack_ServerFolder/SINS_Output/"
  #simulationName = "13x13_K50_m01_hlf_30STR"
  #simulationName = "25x25_25kGen"
  #nSim = 2
  #markers = c("A1","A3","A5")
  #generations = 25000
  #indPerGen = 4
  #layerName = "layer0"
  #n_cores = 2
  #dataFormat = "adegenet"

  fileName = paste(dataDirectory,
                   simulationName,"/simulation_",nSim,"/",
                   sep="")


  #register parallel parameters
  #to use the SNOW parallel system instead of the MULTICORE parallel system,
  #uncomment the line below and feed that variable to registerDoParallel()
  #also uncomment stopCluster() after the foreach loops. SNOW is slower than MULTICORE
  #mkCluster = makeCluster(n_cores)
  registerDoParallel(cores = n_cores)

  myData =
    foreach(i=fileName) %:%
    foreach(m=markers,.packages = c('data.table','adegenet')) %dopar% {

      dataFile = fread(input = paste(i,layerName,"_",m,".txt",sep=""),
                       sep = "auto", header = FALSE, na.strings = "NA",
                       encoding = "UTF-8")

      if(dataFormat == "sins"){
        colnames(dataFile) <- c("gen","X","Y","id","mom","dad",m)
        dataFile$pop = paste0("pop_",dataFile$X,"-",dataFile$Y)
      }else{
        colnames(dataFile) = c("id", m, "X","Y","pop","gen")
      }
      #create a unique ID key for all lines
      dataFile$key = seq_len(nrow(dataFile))
      dataFile = data.table(dataFile, key ="key")

      return(dataFile)

    }

  #stopCluster(mkCluster)

  #merge data so that all markers are in one data structure
  if(dataFormat == "sins"){
    mergeAllMarkers = foreach(i=1:length(myData)) %do%{
      Reduce(function(...) merge(...,by=c("key","id", "X","Y","pop","gen","mom","dad")),myData[[i]])
    }
  }else{
    mergeAllMarkers = foreach(i=1:length(myData)) %do%{
      Reduce(function(...) merge(...,by=c("key","id", "X","Y","pop","gen")),myData[[i]])
    }
  }

  rm(myData)

  numMarkers = length(markers)


  genIndObj<-
    foreach(j=1:length(mergeAllMarkers)) %do%{
      #split data by generation
      splitPerGen = split(mergeAllMarkers[[j]],f = mergeAllMarkers[[j]][[6]])
      #apply genind function on split data so that we have a list of
      #gen ind obj per generation
      mkCluster = makeCluster(n_cores,type = "FORK")
      clusterExport(cl = mkCluster, varlist = c("j","numMarkers","splitPerGen","dataFormat","df2genind"), envir = environment())
      genIndArrPerGen = parLapply(cl = mkCluster, X = splitPerGen, fun = function(df){

        firstMarkerIdx = ncol(df) - numMarkers + 1
        lastMarkerIdx = ncol(df)

        #df2genind takes single matrix of markers only
        genObj = df2genind(X = df[,firstMarkerIdx:lastMarkerIdx],
                           sep = "/", ploidy = 2, ncode = 3,
                           ind.names = df[[2]],
                           pop = df[[5]], type = "codom"
        )
        #add misc information
        if(dataFormat == "sins"){
          genObj@other$mom = df[[7]]
          genObj@other$dad = df[[8]]
        }
        xy_coords = cbind(df[[3]],df[[4]])
        genObj@other$generation = unique(df[[6]])
        genObj@other$xy = xy_coords
        genObj@other$simID = j
        return(genObj)
      })

      stopCluster(mkCluster)
      return(genIndArrPerGen)
    }

  return(genIndObj)
}



#***************
#GET DATA FROM FILES FROM A SINGLE SIMULATION
#***************
#+++++++++++++++
Get_data_from_files_single_sim = function(the_file_names,
                                          set_markers,
                                          generations_list,
                                          path_to_data,
                                          simulation_number,
                                          is_marker_diploid = TRUE) {
  # Gets data from the files with the parameters gotten from the previous functions

  setwd(paste(path_to_data, "Adegenet_sim_", simulation_number, "/", sep =
                ""))

  the_data_list = NULL
  for (i in the_file_names) {
    the_data_list = c(
      the_data_list,
      Read_data_func(
        i,
        set_markers,
        is_marker_diploid = is_marker_diploid,
        simulation_ID = simulation_number
      )
    )
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
                                         set_markers,
                                         generations_list,
                                         path_to_data,
                                         number_of_simulations,
                                         is_marker_diploid = TRUE) {
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

  all_sim = vector(mode = "list", length = number_of_simulations)



  for (j in 1:number_of_simulations) {
    setwd(paste(path_to_data, "/", "Adegenet_sim_", j, "/", sep = ""))

    the_data_list = vector(mode = "list", length = length(the_file_names))
    for (i in 1:length(the_file_names)) {
      #the_data_list = c(the_data_list,Read_data_func(i, set_markers,is_marker_diploid=is_marker_diploid, simulation_ID = j))
      the_data_list[[i]] <-
        Read_data_func(
          the_file_names[i],
          set_markers,
          is_marker_diploid = is_marker_diploid,
          simulation_ID = j
        )
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
Do_summary_calcs = function(simulation_name,
                            the_genind_obj,
                            the_path_to_plot_folder) {
  summary_genind = summary(the_genind_obj)
  names(summary_genind)


  pdf(
    file = paste(
      the_path_to_plot_folder,
      "/",
      simulation_name,
      "_SummaryStats.pdf",
      sep = ""
    ),
    width = 20,
    height = 15
  )
  par(mfrow = c(2, 2), oma = c(0, 0, 1, 0))#sets 2x2 spots for the graphics and sets space for outer title

  plot(
    summary_genind$n.by.pop,
    summary_genind$pop.n.all,
    xlab = "Colonies sample size",
    ylab = "Number of alleles",
    main = "Alleles numbers and sample sizes",
    type = "n",
    xlim = c(
      range(summary_genind$n.by.pop)[1] - 5,
      range(summary_genind$n.by.pop)[2] + 5
    ),
    ylim = c(
      range(summary_genind$pop.n.all)[1] - 5,
      range(summary_genind$pop.n.all)[2] + 5
    )
  )
  text(
    summary_genind$n.by.pop,
    summary_genind$pop.n.all,
    lab = names(summary_genind$n.by.pop)
  )
  barplot(summary_genind$loc.n.all, ylab = "Number of alleles",
          main = "Number of alleles per locus")
  barplot(
    summary_genind$Hexp - summary_genind$Hobs,
    ylab = "Hexp - Hobs",
    col = ifelse(summary_genind$Hexp - summary_genind$Hobs > 0, "red", "blue")
  )
  title(expression(
    phantom("Heterozygosity: ") * "expected" * phantom(" - observed")
  ), col.main = "red")
  title(expression(phantom("Heterozygosity: expected - ") * "observed"), col.main =
          "blue")
  title(expression(
    "Heterozygosity:" * phantom(" expected ") * "-" * phantom(" observed"),
    col.main = "black"
  ))

  barplot(
    summary_genind$n.by.pop,
    main = "Sample sizes per population",
    ylab = "Number of genotypes",
    las = 3
  )
  title(paste("Generation: ", the_genind_obj@other$generation),
        outer = TRUE)


  the_plot = recordPlot()
  dev.off()

  bartlett_test = bartlett.test(list(summary_genind$Hexp, summary_genind$Hobs))

  the_ttest = t.test(
    summary_genind$Hexp,
    summary_genind$Hobs,
    pair = T,
    var.equal = TRUE,
    alter = "greater"
  )

  results = list(summary_genind, bartlett_test, the_ttest, the_plot)

  print("Summary Calcs done.")
  #return(summary_genind)
  return(results)
}

Do_summary_calcs_updated = function(simulation_name,
                                    the_genind_array,
                                    the_path_to_plot_folder,
                                    groupBy = "simulation",
                                    simGroup = NULL,
                                    genGroup = NULL) {
  if (FALSE) {
    simulation_name = name_of_simulation
    the_genind_array = raw_data_multi_sim
    the_path_to_plot_folder = the_path_to_plot_folder
    groupBy = "simulation"
    simGroup = 1
    genGroup = 100


  }

  if (FALSE) {
    mySummary = function(object, verbose = TRUE, ...) {
      x <- object
      if (!is.genind(x))
        stop("Provided object is not a valid genind.")


      if (is.null(pop(x))) {
        pop(x) <- rep("P1", nInd(x))
      }

      ## BUILD THE OUTPUT ##
      ## type-independent stuff
      res <- list()

      res$n <- nrow(myTab(x))

      res$n.by.pop <- as.numeric(table(pop(x)))
      names(res$n.by.pop) <- popNames(x)

      ## PA case ##
      if (x@type == "PA") {
        ## % of missing data
        res$NA.perc <- 100 * sum(is.na(myTab(x))) / prod(dim(myTab(x)))

        return(invisible(res))
      }


      ## codom case ##
      res$loc.n.all <- nAll(x)

      temp <- myTab(genind2genpop(x, quiet = TRUE))

      res$pop.n.all <-
        apply(temp, 1, function(r)
          sum(r != 0, na.rm = TRUE))

      res$NA.perc <- 100 * (1 - mean(propTyped(x, by = "both")))

      ## handle heterozygosity
      if (any(ploidy(x) > 1)) {
        ## auxiliary function to compute observed heterozygosity
        temp <- lapply(seploc(x), tab, freq = TRUE)
        f1 <- function(tab) {
          H <- apply(tab, 1, function(vec)
            any(vec > 0 & vec < 1))
          H <- mean(H, na.rm = TRUE)
          return(H)
        }

        res$Hobs <- unlist(lapply(temp, f1))

        ## auxiliary function to compute expected heterozygosity
        ## freq is a vector of frequencies
        f2 <- function(freq) {
          H <- 1 - sum(freq * freq, na.rm = TRUE)
          return(H)
        }

        temp <- genind2genpop(x, pop = rep(1, nInd(x)), quiet = TRUE)
        temp <- myTab(temp, freq = TRUE, quiet = TRUE)
        res$Hexp <-
          tapply(temp ^ 2, locFac(x), function(e)
            1 - sum(e, na.rm = TRUE))
      } else {
        # no possible heterozygosity for haploid genotypes
        res$Hobs <- 0
        res$Xexp <- 0
      }

      ## add class and return
      class(res) <- "genindSummary"
      return(res)
    }  # end summary.genind
  }

  allData_summary = rapply(
    object = the_genind_array,
    f = function(X) {
      return(summary(X))
    },
    how = 'list'
  )

  #allData_meanSummary = rapply(allData_summary, mean, how='list')

  allData_meanSummary = rapply(allData_summary, mean, how = 'replace')
  lapply(allData_meanSummary, function(x)
    split(x, f = '[['))

  #gen - stat


  lapply(X = allData_meanSummary, FUN = lapply, '[[', 7)


  allData_meanSimSummary = NULL
  allData_meanSimSummary$Hexp = mean(unlist(lapply(
    lapply(allData_meanSummary, '[[', 10), '[[', 7
  )))

  asd = mean(unlist(lapply(
    lapply(allData_meanSummary, '[[', 10), '[[', 7
  )))
  str(asd)


  lapply(allData_meanSummary, function(x) {
    lapply(x, FUN = '[[', 2)
  }, mean)

  allData_meanSummary_1 = rbindlist(allData_meanSummary, fill = T, use.names = T)


  unlist(allData_meanSummary[[1]][[10]],
         recursive = F,
         use.names = T)

  str(allData_meanSummary, max.level = 2)

  allData_meanSummary

  as.data.frame(allData_meanSummary[[1]][[1]]$pop.n.all)

  ggplot(allData_meanSummary)


  pdf(
    file = paste(
      the_path_to_plot_folder,
      "/",
      simulation_name,
      "_SummaryStats.pdf",
      sep = ""
    ),
    width = 20,
    height = 15
  )
  par(mfrow = c(2, 2), oma = c(0, 0, 1, 0))#sets 2x2 spots for the graphics and sets space for outer title

  plot(
    summary_genind$n.by.pop,
    summary_genind$pop.n.all,
    xlab = "Colonies sample size",
    ylab = "Number of alleles",
    main = "Alleles numbers and sample sizes",
    type = "n",
    xlim = c(
      range(summary_genind$n.by.pop)[1] - 5,
      range(summary_genind$n.by.pop)[2] + 5
    ),
    ylim = c(
      range(summary_genind$pop.n.all)[1] - 5,
      range(summary_genind$pop.n.all)[2] + 5
    )
  )
  text(
    summary_genind$n.by.pop,
    summary_genind$pop.n.all,
    lab = names(summary_genind$n.by.pop)
  )
  barplot(summary_genind$loc.n.all, ylab = "Number of alleles",
          main = "Number of alleles per locus")
  barplot(
    summary_genind$Hexp - summary_genind$Hobs,
    ylab = "Hexp - Hobs",
    col = ifelse(summary_genind$Hexp - summary_genind$Hobs > 0, "red", "blue")
  )
  title(expression(
    phantom("Heterozygosity: ") * "expected" * phantom(" - observed")
  ), col.main = "red")
  title(expression(phantom("Heterozygosity: expected - ") * "observed"), col.main =
          "blue")
  title(expression(
    "Heterozygosity:" * phantom(" expected ") * "-" * phantom(" observed"),
    col.main = "black"
  ))

  barplot(
    summary_genind$n.by.pop,
    main = "Sample sizes per population",
    ylab = "Number of genotypes",
    las = 3
  )
  title(paste("Generation: ", the_genind_obj@other$generation),
        outer = TRUE)


  the_plot = recordPlot()
  dev.off()

  bartlett_test = bartlett.test(list(summary_genind$Hexp, summary_genind$Hobs))

  the_ttest = t.test(
    summary_genind$Hexp,
    summary_genind$Hobs,
    pair = T,
    var.equal = TRUE,
    alter = "greater"
  )

  results = list(summary_genind, bartlett_test, the_ttest, the_plot)

  print("Summary Calcs done.")
  #return(summary_genind)
  return(results)
}

################
#COMPUTE HS ON GENIND ARRAY
################
#+++++++++++++++
Compute.Hs.make.df = function (genind_array) {
  compute.Hs.df = data.frame()
  for (i in 1:length(genind_array)) {
    print(i)
    # if there are more than 1 allele per locus compute Hs
    #if(length(genind_array[[i]]@all.names$loc1)>1 || length(genind_array[[i]]@all.names$loc2)>1){
    # Stack helps shape the data so that we can plot it later
    test.var = stack(Hs(genind_array[[i]]))
    test.var.cbind = cbind(test.var, gen = genind_array[[i]]@other$generation)
    compute.Hs.df = rbind(compute.Hs.df, test.var.cbind)
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
Compute.Hs.over.all.sims = function (sim_data) {
  # Note that if for some reasons all the individuals are homozygous while analysing MORE than 1
  # marker, then this function will crash because for some reason, while analyzing 1 marker it can
  # compute that the heterozygosity can be 0, but while calculating the same for several marker
  # and if for one/several markers all individuals are homozygous, then the function will crash.
  all_sims = NULL

  for (j in 1:length(sim_data)) {
    genind_array = sim_data[[j]]

    compute.Hs.df = data.frame()
    for (i in 1:length(genind_array)) {
      # if there are more than 1 allele per locus compute Hs
      #if(length(genind_array[[i]]@all.names$loc1)>1 || length(genind_array[[i]]@all.names$loc2)>1){
      # Stack helps shape the data so that we can plot it later


      test.var = stack(Hs(genind_array[[i]]))

      test.var.cbind = cbind(test.var, gen = genind_array[[i]]@other$generation)
      compute.Hs.df = rbind(compute.Hs.df, test.var.cbind)
      #}
    }
    compute.Hs.df.numeric.gen = compute.Hs.df
    ## transform gen col from factor to numeric so that we are able to use gen as a continuous scale for the plot
    compute.Hs.df.numeric.gen$gen = as.numeric(as.character(compute.Hs.df.numeric.gen$gen))

    all_sims[[j]] = compute.Hs.df.numeric.gen$values
  }

  all_sims = do.call("cbind", all_sims)

  all_sims = data.frame(all_sims)
  all_sims$mean = apply(all_sims, MARGIN = 1, FUN = mean)
  compute.Hs.df.numeric.gen$gen
  all_sims$ind = compute.Hs.df.numeric.gen$ind
  all_sims$gen = compute.Hs.df.numeric.gen$gen

  return(all_sims)
}

################
#COMPUTE HS (PER LOCUS) OVER SEVERAL SIMULATIONS
################
#+++++++++++++++
Compute.Hs.over.all.sims.per.locus = function (sim_data) {
  # Note that if for some reasons all the individuals are homozygous while analysing MORE than 1
  # marker, then this function will crash because for some reason, while analyzing 1 marker it can
  # compute that the heterozygosity can be 0, but while calculating the same for several marker
  # and if for one/several markers all individuals are homozygous, then the function will crash.

  #sim_data = raw_data_multi_sim
  all_sims = NULL

  for (j in 1:length(sim_data)) {

    genind_array = sim_data[[j]]

    compute.Hs.df = data.frame()
    for (i in 1:length(genind_array)) {
      # if there are more than 1 allele per locus compute Hs
      #if(length(genind_array[[i]]@all.names$loc1)>1 || length(genind_array[[i]]@all.names$loc2)>1){
      # Stack helps shape the data so that we can plot it later

      ###TESTING
      hs.sep.loc = lapply(X = seploc(genind_array[[i]]), FUN = Hs) #apply Hs to list
      hs.sep.loc.df = data.frame(hs.sep.loc)
      hs.sep.loc.df$values = apply(hs.sep.loc.df, MARGIN = 1, FUN = mean) #create new col "values" with the mean over all the loci
      new.test.var = subset(hs.sep.loc.df, select = "values") # new df w/ only "values" col
      new.test.var$ind = rownames(new.test.var) # transform row names into a collumn
      rownames(new.test.var) = NULL # transform row names into numbers
      ###/TESTING

      #test.var = stack(Hs(genind_array[[i]]))
      test.var = new.test.var

      test.var.cbind = cbind(test.var, gen = genind_array[[i]]@other$generation)
      compute.Hs.df = rbind(compute.Hs.df, test.var.cbind)
      #}
    }
    compute.Hs.df.numeric.gen = compute.Hs.df
    ## transform gen col from factor to numeric so that we are able to use gen as a continuous scale for the plot
    compute.Hs.df.numeric.gen$gen = as.numeric(as.character(compute.Hs.df.numeric.gen$gen))

    all_sims[[j]] = compute.Hs.df.numeric.gen$values


  }



  all_sims = do.call("cbind", all_sims)

  all_sims = data.frame(all_sims)
  all_sims$mean = apply(all_sims, MARGIN = 1, FUN = mean)
  #compute.Hs.df.numeric.gen$gen
  all_sims$ind = compute.Hs.df.numeric.gen$ind
  all_sims$gen = compute.Hs.df.numeric.gen$gen

  return(all_sims)
}

################
#TRANSFORM INTO SIMPLE POP (SO THAT WE CAN THEN COMPUTE HS ON GENIND ARRAY)
################
#+++++++++++++++
Transform_into_simple_pop = function(genind_array_obj) {
  # Transforms a genind array grouping populations into simpler(read broader) groups
  # pop1A will be pop1
  # its defined for the following population code up to pop15
  Group.demes.function = function(X) {
    ## Group demes according to their population code (eg "pop1A" will be "pop1", "pop6E" will be "pop6")
    X = raw_data_multi_sim[[1]][[1]]
    popNames(X)
    varN = NULL
    for (i in popNames(X)) {
      if (strsplit(i, split = "pop1{1}")[[1]][1] == "") {
        varN = c(varN, "pop1")
      }
    }


    varN = NULL
    for (i in popNames(X)) {
      if (strsplit(i, split = "pop1{1}")[[1]][1] == "") {
        varN = c(varN, "pop1")
      } else if (strsplit(i, split = "pop2{1}")[[1]][1] == "") {
        varN = c(varN, "pop2")
      } else if (strsplit(i, split = "pop3{1}")[[1]][1] == "") {
        varN = c(varN, "pop3")
      } else if (strsplit(i, split = "pop4{1}")[[1]][1] == "") {
        varN = c(varN, "pop4")
      } else if (strsplit(i, split = "pop5{1}")[[1]][1] == "") {
        varN = c(varN, "pop5")
      } else if (strsplit(i, split = "pop6{1}")[[1]][1] == "") {
        varN = c(varN, "pop6")
      } else if (strsplit(i, split = "pop7{1}")[[1]][1] == "") {
        varN = c(varN, "pop7")
      } else if (strsplit(i, split = "pop8{1}")[[1]][1] == "") {
        varN = c(varN, "pop8")
      } else if (strsplit(i, split = "pop9{1}")[[1]][1] == "") {
        varN = c(varN, "pop9")
      } else if (strsplit(i, split = "pop10{1}")[[1]][1] == "") {
        varN = c(varN, "pop10")
      } else if (strsplit(i, split = "pop11{1}")[[1]][1] == "") {
        varN = c(varN, "pop11")
      } else if (strsplit(i, split = "pop12{1}")[[1]][1] == "") {
        varN = c(varN, "pop12")
      } else if (strsplit(i, split = "pop13{1}")[[1]][1] == "") {
        varN = c(varN, "pop13")
      } else if (strsplit(i, split = "pop14{1}")[[1]][1] == "") {
        varN = c(varN, "pop14")
      } else if (strsplit(i, split = "pop15{1}")[[1]][1] == "") {
        varN = c(varN, "pop15")
      }
    }
    popNames(X) = varN
    return(X)
  }

  raw_data.grouped.demes = genind_array_obj

  raw_data.grouped.demes = lapply(X = raw_data.grouped.demes, FUN = Group.demes.function)

  return(raw_data.grouped.demes)
}

################
#MAKE HS PLOTS
################
#+++++++++++++++
Make.HS.plots = function(compute.Hs.df.numeric.gen,
                         isSimple = F,
                         initial.gen = 0,
                         env.change = 0,
                         last.gen) {
  if (!isSimple) {
    hs.plot = ggplot(data = compute.Hs.df.numeric.gen) +
      geom_line(aes(
        x = gen,
        y = values,
        group = ind,
        color = ind
      ), ## group by ind, color by ind
      alpha = 1.0) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      ggtitle("Expected Heterozigosity for each deme") +
      xlab("Generations") +
      ylab("Expected Heterozigosity") +
      theme(legend.position = "right") +
      scale_x_continuous(breaks = seq(initial.gen, last.gen, 1000)) +
      #scale_x_discrete(labels=compute.Hs.df$gen,breaks=compute.Hs.df$gen, limits = compute.Hs.df$gen)+
      geom_vline(xintercept = env.change, linetype = "longdash") +
      coord_cartesian(ylim = c(0, 1),
                      xlim = c(initial.gen, last.gen))
  } else {
    hs.plot = ggplot(data = compute.Hs.df.numeric.gen) +
      geom_line(
        aes(
          x = gen,
          y = values,
          group = ind,
          color = ind
        ),
        ## group by ind, color by ind
        subset(
          compute.Hs.df.numeric.gen,
          ind == "pop1" |
            ind == "pop2" | ind == "pop3" |
            ind == "pop4" |
            ind == "pop5" | ind == "pop6" |
            ind == "pop7" |
            ind == "pop8" |
            ind == "pop9"
        ),
        alpha = 1.0
      ) + ## subset the data, choose which ind(pops) we want to plot
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      ggtitle("Expected Heterozigosity for group of demes") +
      xlab("Generations") +
      ylab("Expected Heterozigosity") +
      theme(legend.position = "right") +
      scale_x_continuous(breaks = seq(initial.gen, last.gen, 1000)) +
      #scale_x_discrete(labels=compute.Hs.df$gen,breaks=compute.Hs.df$gen, limits = compute.Hs.df$gen)+
      geom_vline(xintercept = env.change, linetype = "longdash") +
      coord_cartesian(ylim = c(0, 1),
                      xlim = c(initial.gen, last.gen))

  }

  print("HS plots done.")
  return(hs.plot)
}

################
#MAKE HS MEAN PLOT (MEAN OVER ALL SIMS)
################
#+++++++++++++++
Make.HS.mean.plots = function(compute.Hs.all.sims,
                              initial.gen = 0,
                              env.change = 0,
                              last.gen,
                              subtitle = NULL,
                              showPerDeme = FALSE,
                              number_of_simulations,
                              vector_small_frags,
                              vector_mediumSmall_frags,
                              vector_mediumLarge_frags,
                              vector_large_frags) {

  #compute.Hs.all.sims=hs.allsim.per.locus


  #summarizes data for a given measure by a given group
  summaryAllPop = summarySE(data = compute.Hs.all.sims,
                            measurevar = "mean",
                            groupvars = c("gen"))

  hs.plot = ggplot(summaryAllPop, aes(x = gen, y = mean))

  if (showPerDeme) {
    #remove last character from the pop eg: pop1A turns into pop1
    compute.Hs.all.sims$ind=substr(compute.Hs.all.sims$ind, 1, nchar(compute.Hs.all.sims$ind)-1)

    #atribute to each pop name a given forest frag type
    compute.Hs.all.sims$fSizeType = PopNameDf_to_PopSizeTypeDf(
      vectorPopNames = compute.Hs.all.sims$ind,
      vector_small_frags = vector_small_frags,
      vector_mediumSmall_frags = vector_mediumSmall_frags,
      vector_mediumLarge_frags = vector_mediumLarge_frags,
      vector_large_frags = vector_large_frags)

    #mean over generation and forest fragment size type
    summaryPerPop = summarySE(data = compute.Hs.all.sims,
                              measurevar = "mean",
                              groupvars = c("gen","fSizeType"))


    #define colours for forest frag groups
    cbPalette = c("Large"="#00CC00", "MedLarge"="#0066CC","MedSmall" = "#FFCC00","Small" = "#FF0000")


    hs.plot = hs.plot + geom_line(data = summaryPerPop,
                                  aes(x = gen, y = mean,
                                      #group=ind,
                                      color = fSizeType),
                                  alpha = 0.5) +
      #facet_grid(.~ind)+
      #scale_colour_discrete(name = "Forest \nFragments")+
      scale_colour_manual(name = "Forest \nFragments", values=cbPalette)
  }

  hs.plot = hs.plot + geom_line(size = .5) +
    #geom_point() +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .1, alpha = 0.1) +
    ggtitle(
      label = paste(
        "Mean Expected Heterozigosity over " ,
        number_of_simulations,
        " simulations",
        sep = ""
      ),
      subtitle = subtitle
    ) +
    xlab("Generations") +
    ylab("Expected Heterozigosity") +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          #legend.background = element_rect(size=0.5, linetype="solid",colour ="grey"),
          legend.title = element_text(colour="black", size=8,
                                      face="plain"),
          legend.text = element_text(colour="black", size=7,
                                     face="plain"),
          legend.margin = margin(unit = "cm",t = 0,r = 0,b = 0,l = 0)
          ) +
    scale_x_continuous(breaks = seq(initial.gen, last.gen, 1000)) +
    #scale_x_discrete(labels=compute.Hs.df$gen,breaks=compute.Hs.df$gen, limits = compute.Hs.df$gen)+
    geom_vline(xintercept = env.change, linetype = "longdash") +
    coord_cartesian(ylim = c(0, 1),
                    xlim = c(initial.gen, last.gen))
  #theme(plot.title = element_text(family = 'Helvetica',color = '#666666',face = 'bold',size = 18,hjust = 0.5))

  print("HS mean plot done.")
  return(hs.plot)
}

####### Compute.Hs.perSim.allLocus.PerPop ######
# Compute heterozigosity over time per simulation for single populations
Compute.Hs.perSim.allLocus.PerPop = function(
  name_of_simulation,
  sim_data,
  vector_small_frags = c(2, 4, 5, 7, 14, 15, 16, 17, 23, 26, 27, 28, 33, 34, 35, 41, 46, 47, 48, 49),
  vector_mediumSmall_frags = c(1, 3, 6, 8, 12, 13, 18, 19, 21, 22, 29, 30, 31, 32, 37, 39, 42, 43, 44, 45),
  vector_mediumLarge_frags = c(9, 10, 11, 24, 25, 36, 38, 40),
  vector_large_frags = c(20)){

  # sim_data = raw_data_multi_sim
  # popsToSample = c("pop17A","pop18B","pop24B","pop20K")
  #vector_small_frags = c(2, 4, 5, 7, 14, 15, 16, 17, 23, 26, 27, 28, 33, 34, 35, 41, 46, 47, 48, 49)
  #vector_mediumSmall_frags = c(1, 3, 6, 8, 12, 13, 18, 19, 21, 22, 29, 30, 31, 32, 37, 39, 42, 43, 44, 45)
  #vector_mediumLarge_frags = c(9, 10, 11, 24, 25, 36, 38, 40)
  #vector_large_frags = c(20)



  #S_Frag = paste("pop",vector_small_frags,"A",sep = "")
  #MS_Frag = paste("pop",vector_mediumSmall_frags,"A",sep = "")
  #ML_Frag = paste("pop",vector_mediumLarge_frags,"C",sep = "")
  #L_Frag = paste("pop",vector_large_frags,"M",sep = "")


  S_Frag = paste("pop",vector_small_frags,"A",sep = "")
  MS_Frag = paste("pop",vector_mediumSmall_frags,"B",sep = "")
  ML_Frag = paste("pop",vector_mediumLarge_frags,"C",sep = "")
  L_Frag = paste("pop",vector_large_frags,"H",sep = "")

  popsToSample = c(S_Frag,MS_Frag,ML_Frag,L_Frag)

  # color palette to be used in the plots
  cbPalette = c(L_Frag = "#00CC00",
                ML_Frag = "#0066CC",
                MS_Frag = "#FFCC00",
                S_Frag = "#FF0000")

  all_sims = NULL

  for(i in 1:length(sim_data)){

    genind_array = sim_data[[i]]

    compute.Hs.df = data.frame()
    for(j in 1:length(genind_array)){

      genind_array[[j]] = genind_array[[j]][pop=popsToSample]

      hs.sep.loc = lapply(X = seploc(genind_array[[j]]), FUN = Hs) #apply Hs to list
      hs.sep.loc.df = data.frame(hs.sep.loc)
      hs.sep.loc.df$values = apply(hs.sep.loc.df, MARGIN = 1, FUN = mean) #create new col "values" with the mean over all the loci
      new.test.var = subset(hs.sep.loc.df, select = "values") # new df w/ only "values" col
      new.test.var$ind = rownames(new.test.var) # transform row names into a collumn
      rownames(new.test.var) = NULL # transform row names into numbers

      new.test.var = cbind(new.test.var, gen = as.numeric(as.character(genind_array[[j]]@other$generation)))

      compute.Hs.df = rbind(compute.Hs.df, new.test.var)
    }
    all_sims[[i]] = compute.Hs.df
  }


  names(all_sims) = c(1:length(all_sims))

  all_sims = data.table::rbindlist(l = all_sims,use.names = TRUE,idcol = "SIMID")

  # Assign forest tags and sizes to populations
  forestFrag = NULL
  forestSize = NULL
  carryingCapacity = 100
  for(i in 1:nrow(all_sims)){

    if(any(all_sims[i][[3]] == S_Frag)){
      forestFrag = c(forestFrag, "S_Frag")
      forestSize = c(forestSize, 2*carryingCapacity)
    }else if(any(all_sims[i][[3]] == L_Frag)){
      forestFrag = c(forestFrag, "L_Frag")
      forestSize = c(forestSize, 99*carryingCapacity)
    }else if(any(all_sims[i][[3]] == MS_Frag)){
      forestFrag = c(forestFrag, "MS_Frag")
      forestSize = c(forestSize, 12*carryingCapacity)
    }else if(any(all_sims[i][[3]] == ML_Frag)){
      forestFrag = c(forestFrag, "ML_Frag")
      forestSize = c(forestSize, 36*carryingCapacity)
    }
  }

  all_sims=cbind(all_sims,forestFrag,forestSize)


  # CORRELATION LEVEL PEARSON CORR (MANTEL TEST)
  uniqGen = unique(all_sims$gen)
  uniqSim = unique(all_sims$SIMID)

  perGenCorr = lapply(X= uniqGen, function(Y){
    perSimGenCorr = lapply(X = uniqSim, function(i){
      #print(Y)
      #print(i)
      # Correlation between heterozigosity and forest size
      a = cor.test(all_sims[(all_sims$gen==Y & all_sims$SIMID==i)]$values, all_sims[(all_sims$gen==Y & all_sims$SIMID==i)]$forestSize, exact = TRUE, alternative = "greater", method = "pearson")
      a$gen = Y
      a$simID = i

      corrRes = data.frame(a$p.value,(a$p.value)<0.05,a$gen,a$estimate, a$simID)
      colnames(corrRes) = c("pvalue","p<0.05","generation","correlation", "simID")

      return(corrRes)
    })
    return(perSimGenCorr)
  })


  #merge data into single dataframe
  merged.data.frame = lapply(X = perGenCorr,function(X){ Reduce(function(...) merge(..., all=T), X)})
  merged.data.frame = Reduce(function(...) merge(..., all=T), merged.data.frame)

  #sort df by generation and simID columns
  merged.data.frame = merged.data.frame[order(merged.data.frame[,3], merged.data.frame[,5]),]


  #View(all_sims)

  myPlotCorr = ggplot(merged.data.frame, aes(x=as.factor(generation),y=correlation, colour=simID, fill=`p<0.05`))+
    geom_col(position="dodge",size=1)+
    coord_cartesian(ylim = c(0, 1.01)) +
    ylab("Pearson correlation") +
    xlab("Generations") +
    ggtitle(paste("Correlation between forest size and heterozigosity (",length(unique(all_sims$ind))," populations)",sep = ""),
            subtitle = name_of_simulation)+
    scale_fill_manual(name = "p<0.05",
                        values=c("TRUE" = "green", "FALSE"="red")
                        #values=rainbow(n = length(popsToSample))
    )
  myPlotCorr





  myplotBP = ggplot(all_sims, aes(x=as.factor(gen),y=values,color=forestFrag,linetype=SIMID))+
    geom_boxplot()+
    coord_cartesian(ylim = c(0, 1.01)) +
    ylab("Expected Heterozigosity") +
    xlab("Geographical distance between demes") +
    ggtitle(paste("Heterozigosity over time per simulation for single populations"),
            subtitle = name_of_simulation)+
    scale_colour_manual(name = "Forest \nFragments",
                        values=cbPalette
                         #values=rainbow(n = length(popsToSample))
                        )



  myplotP = ggplot(all_sims, aes(x=as.factor(gen),y=values,color=forestFrag))+
    geom_jitter()+
    coord_cartesian(ylim = c(0, 1.01)) +
    ylab("Expected Heterozigosity") +
    xlab("Geographical distance between demes") +
    ggtitle(paste("Heterozigosity over time per simulation for single populations"),
            subtitle = name_of_simulation)+
    scale_colour_manual(name = "Forest \nFragments",
                        values=cbPalette
                        #values=rainbow(n = length(popsToSample))
    )
  plotList = list(myplotBP,myplotP,myPlotCorr)


  return(plotList)

}



#### TOTAL NUMBER OF ALLELES PER LOCUS VS AVERAGE NUMBER OF ALLELES PER LOCUS PER POP (MEAN OVER ALL SIMS) #######
#+++++++++++++++
Tot.vs.Avg.N.Alleles = function(all.simulations.data,
                                avg.over.all.loci = F) {
  #all.simulations.data = raw_data_multi_sim

  save.generations = NULL
  mean.nbAllelePerSampledDeme = NULL
  tot.num.alleles = NULL

  for (j in 1:length(all.simulations.data[[1]])) {
    for (i in 1:length(all.simulations.data)) {
      # t() transposes the matrix, that is, transforms columns in rows and viceversa
      nbAllelePerSampledDeme = data.frame(t(nb.alleles(all.simulations.data[[i]][[j]])))
      mean.nbAllelePerSampledDeme = rbind(mean.nbAllelePerSampledDeme,
                                          colMeans(nbAllelePerSampledDeme))

      tot.num.alleles = rbind(tot.num.alleles,
                              nAll(all.simulations.data[[i]][[j]]),
                              deparse.level = 0)

      save.generations = rbind(save.generations,
                               all.simulations.data[[i]][[j]]@other$generation)

    }
  }


  numberOfLoci = length(levels(all.simulations.data[[1]][[1]]@loc.fac))

  mean.nbAllelePerSampledDeme = data.frame(mean.nbAllelePerSampledDeme)
  mean.nbAllelePerSampledDeme$generation = as.numeric(save.generations)

  mean.nbAllelePerSampledDeme.simAvg = aggregate(
    mean.nbAllelePerSampledDeme[, 1:numberOfLoci],
    list(Generation = mean.nbAllelePerSampledDeme$generation),
    mean
  )


  tot.num.alleles = data.frame(tot.num.alleles)
  tot.num.alleles$generation = as.numeric(save.generations)

  tot.num.alleles.simAvg = aggregate(tot.num.alleles[, 1:numberOfLoci],
                                     list(Generation = tot.num.alleles$generation),
                                     mean)

  #### TODO SAVE IMAGE TO FILE
  if (FALSE) {
    ggplot2::ggsave(
      file = paste(
        the_path_to_plot_folder,
        "/",
        name_of_simulation,
        "_TotNumAllelesVSAvgNumAlleles",
        length(raw_data_multi_sim),
        "Sims.pdf",
        sep = ""
      ),
      width = 20,
      height = 15
    )
  }

  if (avg.over.all.loci) {
    df.simAlleleAvg = data.frame(avg = rowMeans(mean.nbAllelePerSampledDeme.simAvg[, 2:numberOfLoci +
                                                                                     1]))
    df.simAlleleTotal = data.frame(tot = rowMeans(tot.num.alleles.simAvg[, 2:numberOfLoci +
                                                                           1]))
    df.Avg.Vs.Total = data.frame(df.simAlleleAvg, df.simAlleleTotal)
    df.Avg.Vs.Total$gen = mean.nbAllelePerSampledDeme.simAvg$Generation
    return(
      ggplot(df.Avg.Vs.Total) +
        geom_line(aes(
          x = gen, y = avg, color = "average"
        )) +
        geom_line(aes(
          x = gen, y = tot, color = "total"
        )) +
        xlab("Generations") +
        ylab("Nb Alleles")
    )
  }

  #melting eases the plotting process with ggplot
  mean.nbAllelePerSampledDeme.simAvg.melt = melt(mean.nbAllelePerSampledDeme.simAvg, id.vars = "Generation")
  tot.num.alleles.simAvg.melt = melt(tot.num.alleles.simAvg, id.vars = "Generation")

  df.Avg.Vs.Total = data.frame(mean.nbAllelePerSampledDeme.simAvg.melt, totValue =
                                 tot.num.alleles.simAvg.melt$value)
  df.Avg.Vs.Total.melt = melt(
    data = df.Avg.Vs.Total,
    id.vars = c("Generation", "variable"),
    measure.vars = c("value", "totValue"),
    variable.name = "id.tag"
  )

  return(ggplot(
    df.Avg.Vs.Total.melt,
    aes(
      x = Generation,
      y = value,
      color = variable,
      linetype = id.tag
    )
  ) +
    geom_line())

}

################
#GET GENIND OBJ FROM A GIVEN GENERATION
################
#+++++++++++++++
Get_genind_from_gen = function(a_genind_array, generation) {
  ## simple function to get a genind object from an array, from a given specified generation
  # Args:
  #     a_genind_array - an array of genind objects
  #     generation - a number with the generation that we want to search for
  #
  # Returns:
  #     The genind object with the given generation if it finds it or a message if it doesn't.

  for (i in 1:length(a_genind_array)) {
    current.gen = as.numeric(a_genind_array[[i]]@other$generation)

    if (current.gen == generation) {
      return(a_genind_array[[i]])
    }
  }
  return(paste(
    "Generation '",
    generation,
    "' does not exist in the given genind array",
    sep = ""
  ))
}

################
#COMPUTE FST PLOT
################
#+++++++++++++++
Compute.Fst.plot = function(gen.ind.object,
                            matrix.as.triangle = T,
                            upper.triangle = F,
                            has.in.graph.txt = F) {
  compute.Fst = gen.ind.object
  compute.Fst.generation = gen.ind.object@other$generation
  #compute.Fst = compute.Fst[pop=c("pop1A","pop2A","pop3A","pop4A","pop4B","pop5A","pop5B","pop5C","pop6A","pop6B","pop6C","pop6D","pop6E")]
  popNames(compute.Fst)
  compute.Fst = pairwise.fstb(compute.Fst)##TODO CHECK IS CHANGING FROM pairwise.fst TO pairwise.fstb BRINGS ANY ERROR


  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat) {
    cormat[lower.tri(cormat)] <- NA
    return(cormat)
  }

  # Get lower triangle of the correlation matrix
  get_lower_tri <- function(cormat) {
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }

  get.matrix.triangle = matrix.as.triangle
  get.upper.triangle = upper.triangle
  if (get.matrix.triangle) {
    if (get.upper.triangle) {
      # Skip this to get the whole matrix instead of just one of the triangles
      compute.Fst.upper.tri = get_upper_tri(compute.Fst)
      melted.compute.Fst = melt(
        compute.Fst.upper.tri,
        id.vars = names(compute.Fst),
        na.rm = T,
        factorsAsStrings = T
      )
    } else {
      compute.Fst.lower.tri = get_lower_tri(compute.Fst)
      melted.compute.Fst = melt(
        compute.Fst.lower.tri,
        id.vars = names(compute.Fst),
        na.rm = T,
        factorsAsStrings = T
      )
    }
  } else {
    melted.compute.Fst = melt(
      compute.Fst,
      id.vars = names(compute.Fst),
      na.rm = T,
      factorsAsStrings = T
    )
  }

  ## Check min and max values to adjust color palette limits
  round(max(melted.compute.Fst$value) + 0.05, digits = 1)
  min(melted.compute.Fst$value)
  ## round values in matrix
  melted.compute.Fst$value = round(melted.compute.Fst$value, 3)

  fst.heatmap = ggplot(data = melted.compute.Fst, aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "yellow",
      high = "red",
      mid = "orange",
      midpoint = 0.15,
      limit = c(0.0, 0.30001),
      space = "Lab",
      name = "Pairwise\nFst values"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      size = 12,
      hjust = 1
    )) +
    #geom_rect(inherit.aes = F,aes(xmin=1-0.4, xmax=1+0.4, ymin=1-0.4, ymax=1+0.4),color="#999999", linetype=3,fill=NA)+
    #geom_rect(inherit.aes = F,aes(xmin=2-0.4, xmax=2+0.4, ymin=2-0.4, ymax=2+0.4),color="#999999", linetype=3,fill=NA)+
    #geom_rect(inherit.aes = F,aes(xmin=3-0.4, xmax=3+0.4, ymin=3-0.4, ymax=3+0.4),color="#999999", linetype=3,fill=NA)+
    #geom_rect(inherit.aes = F,aes(xmin=4-0.4, xmax=5+0.4, ymin=4-0.4, ymax=5+0.4),color="#999999", linetype=3,fill=NA)+
    #geom_rect(inherit.aes = F,aes(xmin=6-0.4, xmax=8+0.4, ymin=6-0.4, ymax=8+0.4), color="#999999", linetype=3,fill=NA)+
    #geom_rect(inherit.aes = F,aes(xmin=9-0.4, xmax=13+0.4, ymin=9-0.4, ymax=13+0.4),color="#999999", linetype=3,fill=NA)+
    coord_fixed()
  #fst.heatmap

  fst.heatmap = fst.heatmap +
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
      legend.direction = "horizontal"
    ) +
    guides(fill = guide_colorbar(
      barwidth = 7,
      barheight = 1,
      title.position = "top",
      title.hjust = 0.5
    )) +
    ggtitle(paste("Pairwise Fst - Generation ", compute.Fst.generation))

  if (has.in.graph.txt) {
    fst.heatmap = fst.heatmap + geom_text(aes(Var1, Var2, label = value),
                                          color = "black",
                                          size = 4)
  }

  print("Pairwise Fst heatmap plot done.")
  return(fst.heatmap)
}

################
#COMPUTE PAIRWISE_FST/DISTANCE (IBD) PLOT
################
#+++++++++++++++
Make_pairFst_distance_plots = function(genind_array,
                                       generation,
                                       between_within_small_seperate_plot = F,
                                       vector_small_frags = NULL) {
  #genind_array=raw_data
  #generation=6000
  #between_within_small_seperate_plot=T
  #vector_small_frags = c("pop1A", "pop2A")


  if (between_within_small_seperate_plot == TRUE &&
      is.null(vector_small_frags)) {
    stop("between_within_small_seperate_plot is TRUE but vector_small_frags is NULL")
  }

  # we create two instances of data because the pairwise.fst function takes genind objects instead of genepop
  raw_data_ind = Get_genind_from_gen(a_genind_array = genind_array, generation = generation)
  raw_data_pop = genind2genpop(raw_data_ind, process.other = T, quiet = T)



  distance.geographic = dist(raw_data_pop@other$xy)

  ## make dist obj into dataframe (so that we can plot it)
  Dist.to.df <- function(inDist) {
    if (class(inDist) != "dist")
      stop("wrong input type")
    A <- attr(inDist, "Size")
    B <-
      if (is.null(attr(inDist, "Labels")))
        sequence(A)
    else
      attr(inDist, "Labels")
    if (isTRUE(attr(inDist, "Diag")))
      attr(inDist, "Diag") <- FALSE
    if (isTRUE(attr(inDist, "Upper")))
      attr(inDist, "Upper") <- FALSE
    data.frame(row = B[unlist(lapply(sequence(A)[-1], function(x)
      x:A))],
      col = rep(B[-length(B)], (length(B) - 1):1),
      value = as.vector(inDist))
  }


  compute.Fst.ibd = pairwise.fst(x = raw_data_ind, res.type = "matrix") # if we use the res.type="dist" we lose line/col names
  compute.Fst.ibd = as.dist(compute.Fst.ibd)
  #is.euclid(compute.Fst.ibd)
  #cailliez(compute.Fst.ibd)
  compute.Fst.ibd = Dist.to.df(compute.Fst.ibd)
  distance.geographic.ibd = Dist.to.df(distance.geographic)

  dist.df = data.frame(distance.geographic.ibd, compute.Fst.ibd)

  if (between_within_small_seperate_plot == FALSE) {
    the_plot = ggplot(dist.df) +
      geom_point(aes(x = dist.df$value, y = dist.df$value.1)) +
      geom_smooth(
        aes(x = value, y = value.1),
        method = 'lm',
        se = F,
        fullrange = T
      ) +
      scale_x_continuous(
        expand = c(0.0, 0.0),
        limits = c(0, max(dist.df$value)),
        breaks = seq(0, (max(dist.df$value) + 1), 1)
      ) + #TODO change limits to dynamic values
      scale_y_continuous(expand = c(0.0, 0.0),
                         limits = c(0, max(dist.df$value.1) + 0.11)) +
      coord_cartesian(xlim = c(0, max(dist.df$value) + 1), ylim = c(0, 0.2)) +
      ylab("Paiwise Fst") +
      xlab("Geographical distance between demes") +
      ggtitle(paste("Generation ", raw_data_ind@other$generation))#+
    #cowplot::theme_cowplot()
  } else {
    within_between.df = dist.df


    value.tag = NULL
    size.of.df = nrow(within_between.df)
    #within_between.df= within_between.df
    for (i in 1:size.of.df) {
      if (any(within_between.df[i, 2] == vector_small_frags)) {
        value.tag = c(value.tag, "small_frag")
      } else if (as.numeric(strsplit(as.character(within_between.df[i, 2]), split = "")[[1]][4]) ==
                 as.numeric(strsplit(as.character(within_between.df[i, 1]), split = "")[[1]][4])) {
        value.tag = c(value.tag, "within_frag")
      } else {
        value.tag = c(value.tag, "between_frag")
      }

    }


    within_between.df$value.tag = value.tag


    the_plot = ggplot(within_between.df) +
      geom_point(
        aes(
          x = within_between.df$value,
          y = within_between.df$value.1,
          color = within_between.df$value.tag,
          shape = within_between.df$value.tag
        )
      ) +
      #geom_smooth(aes(x = value, y = value.1), method = 'lm', se = F, fullrange=T)+
      geom_smooth(
        data = subset(x = within_between.df, value.tag == "between_frag") ,
        aes(x = value, y = value.1, color = value.tag),
        method = 'lm',
        se = F,
        fullrange = T
      ) +
      geom_smooth(
        data = subset(x = within_between.df, value.tag == "small_frag") ,
        aes(x = value, y = value.1, color = value.tag),
        method = 'lm',
        se = F,
        fullrange = T
      ) +
      geom_smooth(
        data = subset(x = within_between.df, value.tag == "within_frag") ,
        aes(x = value, y = value.1, color = value.tag),
        method = 'lm',
        se = F,
        fullrange = T
      ) +
      scale_x_continuous(
        expand = c(0.0, 0.0),
        limits = c(0, max(dist.df$value)),
        breaks = seq(0, (max(dist.df$value) + 1), 1)
      ) + #TODO change limits to dynamic values
      scale_y_continuous(expand = c(0.0, 0.0),
                         limits = c(0, max(dist.df$value.1) + 0.11)) +
      coord_cartesian(xlim = c(0, max(dist.df$value) + 1), ylim = c(0, 0.2)) +
      ylab("Paiwise Fst") +
      xlab("Geographical distance between demes") +
      theme(legend.position = c(0.8, 0.8)) +
      ggtitle(paste("Generation ", raw_data_ind@other$generation))
    the_plot


    if (FALSE) {
      #########
      #########
      #########
      # TODO
      genind_array = raw_data
      generation = 6000
      between_within_small_seperate_plot = T
      vector_small_frags = c("pop1A", "pop2A")


      head(within_between.df, n = 10)

      coloring.pop = lapply(X = strsplit(
        x = as.character(within_between.df$col.1),
        split = ""
      ), function(X) {
        return(X[4])
      })

      within_between.df$the.coloring.pop = as.character(coloring.pop)
      within_between.df.subsetGeo = subset(x = within_between.df, value <
                                             14)

      the_plot = ggplot(within_between.df,
                        aes(
                          x = value,
                          y = value.1,
                          colour = value.tag
                        )) +

        geom_point(aes(shape = the.coloring.pop),
                   size = 3,
                   alpha = 0.8) +

        scale_shape_manual(values = c(0:9)) +

        geom_smooth(
          data = subset(x = within_between.df, value.tag == "between_frag"),
          inherit.aes = F,
          aes(
            x = value,
            y = value.1,
            colour = value.tag
          ),
          method = 'lm',
          se = F,
          fullrange = T
        ) +

        geom_smooth(
          data = subset(x = within_between.df, value.tag == "small_frag"),
          inherit.aes = F,
          aes(
            x = value,
            y = value.1,
            colour = value.tag
          ),
          method = 'lm',
          se = F,
          fullrange = T
        ) +

        geom_smooth(
          data = subset(x = within_between.df, value.tag == "within_frag"),
          inherit.aes = F,
          aes(
            x = value,
            y = value.1,
            colour = value.tag
          ),
          method = 'lm',
          se = F,
          fullrange = T
        ) +

        geom_smooth(
          data = subset(x = within_between.df.subsetGeo, value.tag == "between_frag"),
          inherit.aes = F,
          aes(
            x = value,
            y = value.1,
            colour = "between_sub14"
          ),
          method = 'lm',
          se = F,
          fullrange = T
        ) +

        scale_x_continuous(
          expand = c(0.0, 0.0),
          limits = c(0, max(within_between.df$value)),
          breaks = seq(0, (max(
            within_between.df$value
          ) + 1), 1)
        ) +

        scale_y_continuous(expand = c(0.0, 0.0),
                           limits = c(0, max(within_between.df$value.1) + 0.11)) +

        coord_cartesian(xlim = c(0, max(dist.df$value) + 1),
                        ylim = c(0, max(dist.df$value.1) + 0.05)) +
        ylab("Paiwise Fst") +
        xlab("Geographical distance between demes") +
        theme(legend.position = c(0.8, 0.8)) +
        ggtitle(paste("Generation ", raw_data_ind@other$generation))
      the_plot
      #########
      #########
      #########
    }
  }

  print("IBD plot done.")
  return(the_plot)
}

Make_pairFst_distance_plot_QUEMERE_Paper = function() {
  #doi: 10.1111/j.1365-294X.2010.04581.x

  df = data.frame("G", "F", 0.09, 288.6, stringsAsFactors = F)
  colnames(x = df) = c("forest1", "forest2", "pFst", "distPx")
  df

  df = rbind(df, c("G", "H", 0.04, 206.0))
  df = rbind(df, c("G", "B", 0.10, 224.0))
  df = rbind(df, c("G", "A", 0.15, 322.2))
  df = rbind(df, c("G", "C", 0.11, 259.7))
  df = rbind(df, c("G", "D", 0.10, 210.5))
  df = rbind(df, c("G", "E", 0.09, 119.6))
  #df = rbind(df,c("G","I",0.08,235.5))

  df = rbind(df, c("F", "H", 0.13, 319.2))
  df = rbind(df, c("F", "B", 0.05, 432.6))
  df = rbind(df, c("F", "A", 0.14, 481.9))
  df = rbind(df, c("F", "C", 0.04, 304.4))
  df = rbind(df, c("F", "D", 0.03, 170.2))
  df = rbind(df, c("F", "E", 0.03, 173.6))
  #df = rbind(df,c("F","I",0.19,430.8))

  df = rbind(df, c("H", "B", 0.15, 427.7))
  df = rbind(df, c("H", "A", 0.15, 528.9))
  df = rbind(df, c("H", "C", 0.16, 444.5))
  df = rbind(df, c("H", "D", 0.15, 346.6))
  df = rbind(df, c("H", "E", 0.15, 232.0))
  #df = rbind(df,c("H","I",0.12,123.0))

  df = rbind(df, c("B", "A", 0.13, 111.5))
  df = rbind(df, c("B", "C", 0.06, 194.8))
  df = rbind(df, c("B", "D", 0.07, 275.5))
  df = rbind(df, c("B", "E", 0.07, 279.4))
  #df = rbind(df,c("B","I",0.21,426.6))

  df = rbind(df, c("A", "C", 0.13, 190.4))
  df = rbind(df, c("A", "D", 0.14, 313.1))
  df = rbind(df, c("A", "E", 0.15, 352.0))
  #df = rbind(df,c("A","I",0.30,536.0))

  df = rbind(df, c("C", "D", 0.04, 136.6))
  df = rbind(df, c("C", "E", 0.05, 220.9))
  #df = rbind(df,c("C","I",0.23,496.5))

  df = rbind(df, c("D", "E", 0.01, 116.0))
  #df = rbind(df,c("D","I",0.22,425.8))

  #df = rbind(df,c("E","I",0.22,311.3))

  df$forest1 = as.factor(df$forest1)
  df$forest2 = as.factor(df$forest2)
  df$pFst = as.numeric(df$pFst)
  df$distPx = as.numeric(df$distPx)

  #140px = 10km
  #
  #x=(px*10)/140

  df$distKm = as.numeric(lapply(X = df$distPx, function(x) {
    return(newDist = (x * 10) / 140)
  }))
  #str(df)

  the_model = lm(data = df, formula = df$pFst ~ df$distKm)
  #smooth_values=data.frame(predict(the_model,se = T)) ## builds smooth values along with standard error

  library(ggplot2)
  library(plotly)
  p = ggplot(df) + geom_point(aes(x = distKm, y = pFst)) +
    geom_label(aes(
      x = distKm,
      y = pFst,
      label = paste(forest1, forest2)
    )) +
    geom_smooth(
      aes(x = distKm, y = pFst),
      method = 'lm',
      se = F,
      fullrange = T
    )
  p = p +  scale_x_continuous(
    expand = c(0.0, 0.0),
    limits = c(0, max(df$distKm) + 2),
    breaks = seq(0, (max(df$distKm) + 2), 1)
  ) + #TODO change limits to dynamic values
    scale_y_continuous(expand = c(0.0, 0.0),
                       limits = c(0, max(df$pFst) + 0.11)) +
    coord_cartesian(xlim = c(0, 27), ylim = c(0, 0.2)) +
    #geom_abline(slope=the_model$coefficients[[2]],intercept = the_model$coefficients[[1]],color="red")+
    #xlim(0,26)+ylim(0,0.4)+
    ylab("Paiwise Fst") +
    xlab("Geographical distance forests") +
    ggtitle(label = paste("IBD plot - The Daraina Area"),
            subtitle = "Quéméré et al.") +
    annotate(
      "text",
      x = 7,
      y = c(0.19, 0.18),
      label = c(
        paste("slope=", the_model$coefficients[[2]]),
        paste("intercept=", the_model$coefficients[[1]], sep = "")
      )
    )
  p
  #ggplotly()
  return(p)


}

Make_pairFst_distance_plots_allSims = function(array_of_genind_array,
                                               name_of_simulation,
                                               generation,
                                               distanceMethod = "manhattan",
                                               between_within_small_seperate_plot = F,
                                               vector_small_frags = NULL) {
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

  if (between_within_small_seperate_plot == TRUE &&
      is.null(vector_small_frags)) {
    stop("between_within_small_seperate_plot is TRUE but vector_small_frags is NULL")
  }

  ## make dist obj into dataframe (so that we can plot it)
  Dist.to.df <- function(inDist) {
    if (class(inDist) != "dist")
      stop("wrong input type")
    A <- attr(inDist, "Size")
    B <-
      if (is.null(attr(inDist, "Labels")))
        sequence(A)
    else
      attr(inDist, "Labels")
    if (isTRUE(attr(inDist, "Diag")))
      attr(inDist, "Diag") <- FALSE
    if (isTRUE(attr(inDist, "Upper")))
      attr(inDist, "Upper") <- FALSE
    data.frame(row = B[unlist(lapply(sequence(A)[-1], function(x)
      x:A))],
      col = rep(B[-length(B)], (length(B) - 1):1),
      value = as.vector(inDist))
  }


  raw_data_ind = Get_genind_from_gen(a_genind_array = array_of_genind_array[[1]], generation = generation)
  raw_data_pop = genind2genpop(raw_data_ind, process.other = T, quiet = T)


  distance.geographic = dist(raw_data_pop@other$xy, method = distanceMethod)

  Fst.ibd.list = data.frame()
  for (k in 1:length(array_of_genind_array))
  {
    # we create two instances of data because the pairwise.fst function takes genind objects instead of genepop
    raw_data_ind = Get_genind_from_gen(a_genind_array = array_of_genind_array[[k]], generation = generation)


    compute.Fst.ibd = pairwise.fstb(gsp = raw_data_ind) # if we use the res.type="dist" we lose line/col names
    compute.Fst.ibd = as.dist(compute.Fst.ibd)
    #is.euclid(compute.Fst.ibd)
    #cailliez(compute.Fst.ibd)
    compute.Fst.ibd = Dist.to.df(compute.Fst.ibd)
    #class(compute.Fst.ibd)
    if (empty(Fst.ibd.list))
    {
      Fst.ibd.list = cbind(as.numeric(compute.Fst.ibd$value))
    }
    else
    {
      Fst.ibd.list = cbind(Fst.ibd.list, as.numeric(compute.Fst.ibd$value))
    }

  }

  #head(data.frame(Fst.ibd.list))
  Fst.ibd.df = data.frame(Fst.ibd.list)
  Fst.ibd.df$mean = apply(Fst.ibd.df, 1, mean, na.rm = T) #does the mean over all the col (simulations) values for each row
  #head(Fst.ibd.df)

  distance.geographic.ibd = Dist.to.df(distance.geographic)
  dist.df = data.frame(distance.geographic.ibd, Fst.ibd.df)
  #head(dist.df)

  if (between_within_small_seperate_plot == FALSE) {
    the_plot = ggplot(dist.df) +
      geom_point(aes(x = dist.df$value, y = dist.df$mean)) +
      geom_smooth(
        aes(x = value, y = mean),
        method = 'lm',
        se = F,
        fullrange = T
      ) +
      scale_x_continuous(
        expand = c(0.0, 0.0),
        limits = c(0, max(dist.df$value)),
        breaks = seq(0, (max(dist.df$value) + 1), 1)
      ) + #TODO change limits to dynamic values
      scale_y_continuous(expand = c(0.0, 0.0),
                         limits = c(0, max(dist.df$mean) + 0.1)) +
      coord_cartesian(xlim = c(0, max(dist.df$value) + 1), ylim = NULL) +
      ylab("Paiwise Fst") +
      xlab("Geographical distance between demes") +
      ggtitle(paste("Generation ", raw_data_ind@other$generation),
              subtitle = name_of_simulation)#+
    #cowplot::theme_cowplot()
  } else {
    within_between.df = dist.df


    value.tag = NULL
    size.of.df = nrow(within_between.df)
    #within_between.df= within_between.df
    for (i in 1:size.of.df) {
      if (any(within_between.df[i, 2] == vector_small_frags)) {
        value.tag = c(value.tag, "small_frag")
      } else if (as.numeric(strsplit(as.character(within_between.df[i, 2]), split = "")[[1]][4]) ==
                 as.numeric(strsplit(as.character(within_between.df[i, 1]), split = "")[[1]][4])) {
        value.tag = c(value.tag, "within_frag")
      } else {
        value.tag = c(value.tag, "between_frag")
      }

    }


    within_between.df$value.tag = value.tag
    #head(within_between.df, n = 100)

    the_plot = ggplot(within_between.df) +
      geom_point(
        aes(
          x = within_between.df$value,
          y = within_between.df$mean,
          color = within_between.df$value.tag,
          shape = within_between.df$value.tag
        )
      ) +
      #geom_smooth(aes(x = value, y = value.1), method = 'lm', se = F, fullrange=T)+
      geom_smooth(
        data = subset(x = within_between.df, value.tag == "between_frag") ,
        aes(x = value, y = mean, color = value.tag),
        method = 'lm',
        se = F,
        fullrange = T
      ) +
      geom_smooth(
        data = subset(x = within_between.df, value.tag == "small_frag") ,
        aes(x = value, y = mean, color = value.tag),
        method = 'lm',
        se = F,
        fullrange = T
      ) +
      geom_smooth(
        data = subset(x = within_between.df, value.tag == "within_frag") ,
        aes(x = value, y = mean, color = value.tag),
        method = 'lm',
        se = F,
        fullrange = T
      ) +
      scale_x_continuous(
        expand = c(0.0, 0.0),
        limits = c(0, max(dist.df$value)),
        breaks = seq(0, (max(dist.df$value) + 1), 1)
      ) + #TODO change limits to dynamic values
      scale_y_continuous(expand = c(0.0, 0.0),
                         limits = c(0, max(dist.df$mean) + 0.11)) +
      coord_cartesian(xlim = c(0, max(dist.df$value) + 1), ylim = c(0, 0.2)) +
      ylab("Paiwise Fst") +
      xlab("Geographical distance between demes") +
      theme(legend.position = c(0.8, 0.8)) +
      ggtitle(paste("Generation ", raw_data_ind@other$generation),
              subtitle = name_of_simulation)
    the_plot

  }

  print("IBD plot done.")
  return(the_plot)
}

pairFst_dist_plots_allSims_fourFragClasses = function(array_of_genind_array,
                                                      generation,
                                                      name_of_simulation,
                                                      distanceMethod = "manhattan",
                                                      vector_small_frags = c(2, 4, 5, 7, 14, 15, 16, 17, 23, 26, 27, 28, 33, 34, 35, 41, 46, 47, 48, 49),
                                                      vector_mediumSmall_frags = c(1, 3, 6, 8, 12, 13, 18, 19, 21, 22, 29, 30, 31, 32, 37, 39, 42, 43, 44, 45),
                                                      vector_mediumLarge_frags = c(9, 10, 11, 24, 25, 36, 38, 40),
                                                      vector_large_frags = c(20)) {
  ##test
  #array_of_genind_array=raw_data_multi_sim[1:2]
  #str(array_of_genind_array,give.head = T,max.level = 2)
  #dim(raw_data_multi_sim)
  #dimnames(raw_data_multi_sim)
  #class(raw_data_multi_sim[1,2][[1]])
  #class(raw_data_single_sim[[1]]@other$generation)
  #raw_data_multi_sim[1,1]
  #raw_data_multi_sim[[1]][[2]]@other$generation

  #generation=20500
  #between_within_small_seperate_plot = T
  #vector_small_frags = c(2,4,5,7,14,15,16,17,23,26,27,28,33,34,35,41,46,47,48,49)
  #vector_mediumSmall_frags = c(1,3,6,8,12,13,18,19,21,22,29,30,31,32,37,39,42,43,44,45)
  #vector_mediumLarge_frags= c(9,10,11,24,25,36,38,40)
  #vector_large_frags = c(20)
  #distanceMethod = "manhattan"
  ##

  if (between_within_small_seperate_plot == TRUE &&
      is.null(vector_small_frags)) {
    stop("between_within_small_seperate_plot is TRUE but vector_small_frags is NULL")
  }



  raw_data_ind = Get_genind_from_gen(a_genind_array = array_of_genind_array[[1]], generation = generation)
  raw_data_pop = genind2genpop(raw_data_ind, process.other = T, quiet = T)

  #raw_data_ind@other$xy
  distance.geographic = dist(raw_data_pop@other$xy, method = distanceMethod)

  #Fst.ibd.list = data.frame()

  Dgen.list = NULL

  raw_data_pop = NULL

  pairwiseFstMatrixList = vector(mode = "list", length = length(array_of_genind_array))

  for (k in 1:length(array_of_genind_array))
  {
    # we create two instances of data because the pairwise.fst function takes genind objects instead of genepop
    raw_data_ind = Get_genind_from_gen(a_genind_array = array_of_genind_array[[k]], generation = generation)


    compute.Fst.ibd_newMethod = pairwise.fstb(gsp = raw_data_ind)
    compute.Fst.ibd_newMethod = as.dist(compute.Fst.ibd_newMethod,
                                        diag = F,
                                        upper = F)
    pairwiseFstMatrixList[[k]] = (compute.Fst.ibd_newMethod)

    ####compute.Fst.ibd = pairwise.fst(x = raw_data_ind, res.type = "matrix") # if we use the res.type="dist" we lose line/col names
    ####compute.Fst.ibd = as.dist(compute.Fst.ibd)
    #is.euclid(compute.Fst.ibd)
    #cailliez(compute.Fst.ibd)
    ####compute.Fst.ibd = Dist.to.df(compute.Fst.ibd)
    #class(compute.Fst.ibd)


  }

  # head(data.frame(Fst.ibd.list))
  # Mean over each of the pairwise values for all the simulations
  pairwiseFstReduce = Reduce('+', pairwiseFstMatrixList)
  pairwiseFstMean = pairwiseFstReduce / length(pairwiseFstMatrixList)
  # head(pairwiseFstMean)
  pairwiseFstMean_df = Dist.to.df(pairwiseFstMean)

  #Fst.ibd.df = data.frame(Fst.ibd.list)
  #Fst.ibd.df$mean = apply(Fst.ibd.df,1,mean,na.rm=T) #does the mean over all the col (simulations) values for each row
  #head(Fst.ibd.df)

  distance.geographic.ibd = Dist.to.df(distance.geographic)
  #dist.df = data.frame(distance.geographic.ibd,Fst.ibd.df)

  #geographic and genetic data into one dataframe. Mean values of pairwiseFst renamed to "mean"
  dist.df = data.frame(distance.geographic.ibd, mean = pairwiseFstMean_df$value)


  # the following line separates the subpops into pops, for example pop13A will be "13" and pop1B will be "1"
  #strsplit(as.character(within_between.df[1432,2]),split = "[aA-zZ]+")[[1]][2]
  #strsplit(as.character(within_between.df[1432,1]),split = "[aA-zZ]+")[[1]][2]

  dist.df= Classify_Pop_Comparison(dist.df,
                                   vector_small_frags = vector_small_frags,
                                   vector_mediumSmall_frags = vector_mediumSmall_frags,
                                   vector_mediumLarge_frags = vector_mediumLarge_frags,
                                   vector_large_frags = vector_large_frags)

  # CORRELATION LEVEL PEARSON CORR (MANTEL TEST)
  #names of population comparisons
  mypops = unique(dist.df$popVS)

  #pearson correlation for each population comparison
  perPopCorr = lapply(X= mypops, function(X){
    a = cor.test(dist.df[dist.df$popVS == X,]$value,dist.df[dist.df$popVS == X,]$mean,exact =TRUE,alternative = "greater",method = "pearson")
    a$pop = X
    return(a)
  })
  #get relevant values
  perPopResults = sapply(perPopCorr, function(X){c(X$p.value,(X$p.value)<0.05,X$pop,X$estimate )})
  #add correlation that uses all the populations
  allPopCorr = cor.test(dist.df$value, dist.df$mean, exact =TRUE,alternative = "greater",method = "pearson")
  allMyResults=cbind(perPopResults, c(allPopCorr$p.value,(allPopCorr$p.value<0.05),"Overall",allPopCorr$estimate))
  #transform data into dataframe format
  allMyResults = t(allMyResults)

  allMyResults=data.frame(allMyResults)
  colnames(allMyResults) = c("pvalue","p<0.05","popComparison","correlation")
  allMyResults$correlation = as.numeric(levels(allMyResults$correlation))[allMyResults$correlation]
  allMyResults$pvalue = as.numeric(levels(allMyResults$pvalue))[allMyResults$pvalue]
  allMyResults$pvalue = as.numeric(format(round(allMyResults$pvalue, 5), nsmall = 5))

  # CORRELATION INFORMATION INTO MAIN DF
  dist.df$pvalue = allMyResults[match(dist.df$popVS,allMyResults$popComparison),1]
  dist.df$p0.05 = allMyResults[match(dist.df$popVS,allMyResults$popComparison),2]
  dist.df$corr = allMyResults[match(dist.df$popVS,allMyResults$popComparison),4]



  the_plot_all = ggplot(dist.df) +
    #geom_point(aes(x = value,y = mean,color = popVS), alpha = 0.5) +
    geom_smooth(
      color = "black",
      aes(x = value, y = mean, color = "Overall"),
      method = 'lm',
      se = F,
      fullrange = T
    )+
    geom_boxplot(
      color = "purple",
      alpha = 0.1,
      # we can do boxplots on continuous vars as long as we group them
      aes(x = value, y = mean, color = popVS, group = value)
    )



  the_plot_between = ggplot() +
    geom_smooth(data = dist.df,
      #color = "black",
      aes(x = value, y = mean, color = "Overall"),
      method = 'lm',
      se = F,
      fullrange = T
    ) +
    stat_smooth(
      data = subset(x = dist.df, popVS == "Within") ,
      aes(x = value, y = mean, color = popVS),
      method = 'lm',
      se = F,
      fullrange = T
    ) +
    geom_smooth(
      data = subset(x = dist.df, popVS == "SvsS") ,
      aes(x = value, y = mean, color = popVS),
      method = 'lm',
      se = F,
      fullrange = T
    ) +
    geom_smooth(
      data = subset(x = dist.df, popVS == "SvsMS") ,
      aes(x = value, y = mean, color = popVS),
      method = 'lm',
      se = F,
      fullrange = T
    ) +
    geom_smooth(
      data = subset(x = dist.df, popVS == "SvsML") ,
      aes(x = value, y = mean, color = popVS),
      method = 'lm',
      se = F,
      fullrange = T
    ) +
    geom_smooth(
      data = subset(x = dist.df, popVS == "SvsL") ,
      aes(x = value, y = mean, color = popVS),
      method = 'lm',
      se = F,
      fullrange = T
    ) +
    geom_smooth(
      data = subset(x = dist.df, popVS == "MSvsMS") ,
      aes(x = value, y = mean, color = popVS),
      method = 'lm',
      se = F,
      fullrange = T
    ) +
    geom_smooth(
      data = subset(x = dist.df, popVS == "MSvsML") ,
      aes(x = value, y = mean, color = popVS),
      method = 'lm',
      se = F,
      fullrange = T
    ) +
    geom_smooth(
      data = subset(x = dist.df, popVS == "MSvsL") ,
      aes(x = value, y = mean, color = popVS),
      method = 'lm',
      se = F,
      fullrange = T
    ) +
    geom_smooth(
      data = subset(x = dist.df, popVS == "MLvsML") ,
      aes(x = value, y = mean, color = popVS),
      method = 'lm',
      se = F,
      fullrange = T
    ) +
    geom_smooth(
      data = subset(x = dist.df, popVS == "MLvsL") ,
      aes(x = value, y = mean, color = popVS),
      method = 'lm',
      se = F,
      fullrange = T
    ) +
    geom_smooth(
      data = subset(x = dist.df, popVS == "LvsL") ,
      aes(x = value, y = mean, color = popVS),
      method = 'lm',
      se = F,
      fullrange = T
    )


  the_plot_box = ggplot(dist.df) +
    geom_boxplot(
    data = subset(x = dist.df, popVS == c("Within")) ,
    aes(x = value, y = mean, color = popVS, group = value)
  ) +
    geom_boxplot(
      data = subset(x = dist.df, popVS == c("SvsS")) ,
      aes(x = value, y = mean, color = popVS, group = value)
    ) +
    geom_boxplot(
      data = subset(x = dist.df, popVS == c("SvsMS")) ,
      aes(x = value, y = mean, color = popVS, group = value)
    ) +
    geom_boxplot(
      data = subset(x = dist.df, popVS == c("SvsML")) ,
      aes(x = value, y = mean, color = popVS, group = value)
    ) +
    geom_boxplot(
      data = subset(x = dist.df, popVS == c("SvsL")) ,
      aes(x = value, y = mean, color = popVS, group = value)
    ) +
    geom_boxplot(
      data = subset(x = dist.df, popVS == c("MSvsMS")) ,
      aes(x = value, y = mean, color = popVS, group = value)
    ) +
    geom_boxplot(
      data = subset(x = dist.df, popVS == c("MSvsML")) ,
      aes(x = value, y = mean, color = popVS, group = value)
    ) +
    geom_boxplot(
      data = subset(x = dist.df, popVS == c("MSvsL")) ,
      aes(x = value, y = mean, color = popVS, group = value)
    ) +
    geom_boxplot(
      data = subset(x = dist.df, popVS == c("MLvsML")) ,
      aes(x = value, y = mean, color = popVS, group = value)
    ) +
    geom_boxplot(
      data = subset(x = dist.df, popVS == c("MLvsL")) ,
      aes(x = value, y = mean, color = popVS, group = value)
    ) +
    geom_boxplot(
      data = subset(x = dist.df, popVS == c("LvsL")) ,
      aes(x = value, y = mean, color = popVS, group = value)
    )


  if (FALSE) {
    new_data_test = subset(x = within_between.df, value.tag.between == "within_frag")## subset data
    new_model = lm(data = new_data_test,
                   formula = new_data_test$mean ~ new_data_test$value) ## create linear regression model gen ~ geo



    ## method to build geom_smooth of the data without actually using geom_smooth

    smooth_values = data.frame(predict(new_model, se = T)) ## builds smooth values along with standard error



    df <- data.frame(
      cbind(
        value = new_data_test$value,
        mean = new_data_test$mean,
        fit = smooth_values$fit,
        upperBound = smooth_values$fit + 2 * smooth_values$se.fit,
        lowerBound = smooth_values$fit - 2 * smooth_values$se.fit
      )
    )

    g <- ggplot(df, aes(value, mean))
    #g <- g + geom_point()
    g <-
      g + geom_linerange(aes(ymin = lowerBound, ymax = upperBound))
    g <- g + geom_point(aes(value, fit))
    g <- g + geom_smooth(method = "lm",
                         color = "green",
                         size = 2)
    g <-
      g + geom_abline(
        slope = new_model$coefficients[[2]],
        intercept = new_model$coefficients[[1]],
        color = "red"
      )
    g
  }


  the_plots_list = list(the_plot_all,
                        the_plot_between,
                        the_plot_box
                        #the_plot_within,
                        #the_plot_withinBetween
                        )

  #color pallete to be used in the plots
  #TODO: compare and make consistent with palette from Hs plots
  forestComparisonPallete = c("SvsS" = "#FF0700",
                              "SvsMS" = "#FF4F00",
                              "SvsML" = "#EC4900",
                              "SvsL" = "#FF7C00",
                              "MSvsMS" = "#FFA700",
                              "MSvsML" = "#FFBF00",
                              "MSvsL" = "#FFD600",
                              "MLvsML" = "#A9F500",
                              "MLvsL" = "#00DA1B",
                              "LvsL" = "#0B7AC3",
                              "Within" = "#585F5E",
                              "Overall" = "#000000")

  for (i in 1:length(the_plots_list)) {

    the_plots_list[[i]] =
      the_plots_list[[i]] +
      scale_x_continuous(
        expand = c(0.0, 0.0),
        limits = c(0, max(dist.df$value) + 1),
        breaks = seq(0, (max(dist.df$value) + 1), 2)
      ) +
      scale_y_continuous(expand = c(0.0, 0.0),
                         limits = c(0, max(dist.df$mean) + 0.11)) +
      #TODO: possible make a dynamic ylim
      coord_cartesian(xlim = c(0, max(dist.df$value) + 1), ylim = c(0, 1.01)) +
      ylab("Paiwise Fst") +
      xlab("Geographical distance between demes") +
      #scale_size(range=c(5,20))+
      theme(#legend.position = c(0.83, 0.7),
            legend.position = "right",
            #legend.background = element_rect(size=0.5, linetype="solid",colour ="grey"),
            legend.title = element_text(colour="black", size=8,
                                       face="plain"),
            legend.text = element_text(colour="black", size=8,
                                       face="plain")) +
      ggtitle(paste("IBD - Generation ", raw_data_ind@other$generation),
              subtitle = name_of_simulation)+
      scale_color_manual(name = "Fragment\ncomparison",
                         values = forestComparisonPallete
                         #labels = c("<= 17", "17 < qsec <= 19", "> 19")
      )
  }


  sortedAllMyResults = allMyResults[order(allMyResults$correlation),]
  #resultsSummaryTable = tableGrob(sortedAllMyResults,rows = NULL)

  tableAsPlot = ggplot(sortedAllMyResults,aes(x=reorder(x = popComparison,X = correlation),y=correlation))+
    geom_col(aes(fill=popComparison))+
    geom_text(aes(label=ifelse(pvalue<0.05 ,"*",'')),colour="black",
               position = position_nudge(y = 0.04))#+
    #annotate("text", label = "* = (pvalue<0.05)", x = 7, y = 0.6, size = 5, colour = "black")

  tableAsPlot = tableAsPlot+
    ylab("Pearson correlation") +
    xlab("Fragment comparison") +
    ggtitle(paste("IBD Statistics - Generation ", raw_data_ind@other$generation),
            subtitle = name_of_simulation)+
    scale_fill_manual(name = "Fragment\ncomparison",
                      values = forestComparisonPallete)+
    theme(legend.position = "none")


  the_plots_list[[4]]=tableAsPlot

  print("IBD plot done.")

  return(the_plots_list)
}


plotFourFragClassesIBD_toFile = function(ibd.plots.list, the_path_to_plot_folder, name_of_simulation){

  lapply(X = ibd.plots.list,function(X){

    lapply(seq_along(X), function(i){
      if(i == 1){
        gen = strsplit(x = X[[i]]$labels$title, split = "  ")[[1]][2]
        save_the_plot(
          plotObject = X[[i]],
          pathToPlot = the_path_to_plot_folder,
          nameOfPlotFile = paste(name_of_simulation, "_IBDplot_All_", gen,
                                 sep = "")
        )
      }else if(i == 2){
        gen = strsplit(x = X[[i]]$labels$title, split = "  ")[[1]][2]
        save_the_plot(
          plotObject = X[[i]],
          pathToPlot = the_path_to_plot_folder,
          nameOfPlotFile = paste(name_of_simulation, "_IBDplot_between_", gen,
                                 sep = "")
        )
      }else if(i == 3){
        gen = strsplit(x = X[[i]]$labels$title, split = "  ")[[1]][2]
        save_the_plot(
          plotObject = X[[i]],
          pathToPlot = the_path_to_plot_folder,
          nameOfPlotFile = paste(name_of_simulation, "_IBDplot_boxplot_", gen,
                                 sep = "")
        )
      }else if(i == 4){
        gen = strsplit(x = X[[i]]$labels$title, split = "  ")[[1]][2]
        save_the_plot(
          plotObject = X[[i]],
          pathToPlot = the_path_to_plot_folder,
          nameOfPlotFile = paste(name_of_simulation, "_IBD_correlation_", gen,
                                 sep = "")
        )
      }
    })

  })


}

################
#COMPUTE PAIRWISE_FST/DISTANCE (IBD) PLOT SLOPE OVER TIME
################
#+++++++++++++++
Compute.PairwiseFst.IBD.Slope.OverTime = function(array_of_genind_array,
                                                  name_of_simulation,
                                                  distanceMethod = "manhattan",
                                                  isPairwiseFstdist = TRUE,
                                                  env.change = 0,
                                                  vector_small_frags = c(2, 4, 5, 7, 14, 15, 16, 17, 23, 26, 27, 28, 33, 34, 35, 41, 46, 47, 48, 49),
                                                  vector_mediumSmall_frags = c(1, 3, 6, 8, 12, 13, 18, 19, 21, 22, 29, 30, 31, 32, 37, 39, 42, 43, 44, 45),
                                                  vector_mediumLarge_frags = c(9, 10, 11, 24, 25, 36, 38, 40),
                                                  vector_large_frags = c(20)) {
  #array_of_genind_array=raw_data_multi_sim[1:2]
  #distanceMethod = "manhattan"
  #isPairwiseFstdist = TRUE
  #name_of_simulation=name_of_simulation
  #env.change=20000
  #*******

  generation_list = NULL
  Dgen.mean = NULL

  #length(levels(array_of_genind_array[[4]][[1]]@pop))
  #length(levels(selPopSize(array_of_genind_array[[4]][[1]], nMin = 2)@pop))


  for (j in 1:length(array_of_genind_array[[1]])) {
    # 1:length(array_of_genind_array[[1]])
    Dgen.list = NULL
    for (i in 1:length(array_of_genind_array)) {

      #drop POPs that have less than 2 individuals (related to problems with pairwise.fstb)
      data_with_dropPop = selPopSize(array_of_genind_array[[i]][[j]], nMin = 2)

      # choose if we want dist vs pairwiseFst method
      if (isPairwiseFstdist) {
        compute.gen.dist = pairwise.fstb(gsp = data_with_dropPop)
        compute.gen.dist = as.dist(compute.gen.dist , diag = F, upper = F)
      } else{
        raw_data_pop = genind2genpop(data_with_dropPop,
                                     process.other = T,
                                     quiet = T)
        compute.gen.dist = dist(raw_data_pop)
      }
      Dgen.list[[i]] = (compute.gen.dist)

      #str(Dgen.list)
      #apply(simplify2array(Dgen.list), c(1,2),mean)

    }

    Dgen.reduce = Reduce('+', Dgen.list)
    Dgen.mean[[j]] = Dgen.reduce / length(Dgen.list)
    generation_list = c(generation_list,
                        as.numeric(array_of_genind_array[[1]][[j]]@other$generation))

  }


  ## Calc geographic distance FROM THE LAST GENERATION GENIND OBJECT
  raw_data_pop = genind2genpop(array_of_genind_array[[1]][[length(array_of_genind_array[[1]])]],
                               process.other = T,
                               quiet = T)

  distance.geographic = dist(raw_data_pop@other$xy, method = distanceMethod)
  distance.geographic = Dist.to.df(distance.geographic)

  Dgen.to.df.allSims = NULL

  for (i in 1:length(Dgen.mean)) {
    Dgen.mean[[i]] = as.dist(Dgen.mean[[i]])
    if (empty(Dgen.to.df.allSims)) {
      Dgen.mean.df = Dist.to.df(Dgen.mean[[i]])
      Dgen.to.df.allSims = data.frame(
        row = Dgen.mean.df$row,
        col = Dgen.mean.df$col,
        generationMean = Dgen.mean.df$value
      )
    } else{
      Dgen.mean.df = Dist.to.df(Dgen.mean[[i]])
      Dgen.to.df.allSims = data.frame(Dgen.to.df.allSims, generationMean =
                                        Dgen.mean.df$value)
    }
  }


  new_data_test = data.frame(distance.geographic)
  new_data_test = data.frame(new_data_test, Dgen.to.df.allSims[, 3:length(Dgen.to.df.allSims)])


  new_data_test = Classify_Pop_Comparison(new_data_test,
                                          vector_small_frags = vector_small_frags,
                                          vector_mediumSmall_frags = vector_mediumSmall_frags,
                                          vector_mediumLarge_frags = vector_mediumLarge_frags,
                                          vector_large_frags = vector_large_frags)

  #plot data pFst vs geodist as in normal IBD plot
  if (FALSE) {
    #just testing stuff
    ggplot(new_data_test) +
      geom_smooth(
        aes(
          x = new_data_test$value,
          y = new_data_test$generationMean.19,
          col = "19"
        ),
        method = "lm",
        se = T
      ) +
      geom_point(aes(
        x = new_data_test$value,
        y = new_data_test$generationMean.19,
        col = "19"
      )) +
      geom_smooth(
        aes(
          x = new_data_test$value,
          y = new_data_test$generationMean,
          col = "1",
          size = 2
        ),
        method = "lm",
        se = T
      ) +
      geom_point(aes(
        x = new_data_test$value,
        y = new_data_test$generationMean,
        col = "1"
      )) +
      geom_abline(slope = 0.0006182,
                  intercept = 0.0351354,
                  color = "black")
  }

  if(FALSE){
  slope.df = NULL
  intercept.df = NULL
  all_model_information = NULL


  # excluding first 4 cols because they correspond to "row", "col", "value" and the first generation
  # first gen full of NaNs
  # exclude last col since its the popVS col.
  for (i in 5:(length(new_data_test)-1)) {
    new_model_all = lm(data = new_data_test,
                   formula = new_data_test[, i] ~ new_data_test$value) ## create linear regression model gen ~ geo

    #all_model_information[[i - 3]] = new_model
    slope.df[i - 4] = new_model_all$coefficients[[2]]
    intercept.df[i - 4] = new_model_all$coefficients[[1]]

  }


  #head(slope.df[1:20,])
  slope.df = data.frame(mean = slope.df, generation_list = generation_list[2:length(generation_list)])
  intercept.df = data.frame(mean = intercept.df, generation_list = generation_list[2:length(generation_list)])
  #head(slope.df)
  #length(all_model_information)
  }



  list_slope_intercept_df = Coefficients_per_comparison_class(new_data_test,generation_list)

  slope_intercept_df = ldply(list_slope_intercept_df, data.frame)

  forestComparisonPallete = c("SvsS" = "#FF0700",
                              "SvsMS" = "#FF4F00",
                              "SvsML" = "#EC4900",
                              "SvsL" = "#FF7C00",
                              "MSvsMS" = "#FFA700",
                              "MSvsML" = "#FFBF00",
                              "MSvsL" = "#FFD600",
                              "MLvsML" = "#A9F500",
                              "MLvsL" = "#00DA1B",
                              "LvsL" = "#0B7AC3",
                              "Within" = "#585F5E",
                              "All" = "black")

  p.slope = ggplot(slope_intercept_df) +
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "SvsS"), aes(x=generation,y=slope,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "SvsMS"), aes(x=generation,y=slope,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "SvsML"), aes(x=generation,y=slope,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "SvsL"), aes(x=generation,y=slope,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "MSvsMS"), aes(x=generation,y=slope,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "MSvsML"), aes(x=generation,y=slope,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "MSvsL"), aes(x=generation,y=slope,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "MLvsML"), aes(x=generation,y=slope,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "MLvsL"), aes(x=generation,y=slope,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "LvsL"), aes(x=generation,y=slope,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "Within"), aes(x=generation,y=slope,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "All"), aes(x=generation,y=slope,color=.id))+
    geom_vline(xintercept = env.change, linetype = "longdash") +
    coord_cartesian(xlim = NULL,
                  ylim = NULL,
                  expand = T) + ggtitle(label = "IBD (Fst/Geographical distance) slope over time", subtitle = name_of_simulation)+
    scale_color_manual(name = "Fragment\ncomparison",
                       values = forestComparisonPallete)

  p.intercept = ggplot(slope_intercept_df) +
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "SvsS"), aes(x=generation,y=intercept,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "SvsMS"), aes(x=generation,y=intercept,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "SvsML"), aes(x=generation,y=intercept,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "SvsL"), aes(x=generation,y=intercept,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "MSvsMS"), aes(x=generation,y=intercept,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "MSvsML"), aes(x=generation,y=intercept,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "MSvsL"), aes(x=generation,y=intercept,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "MLvsML"), aes(x=generation,y=intercept,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "MLvsL"), aes(x=generation,y=intercept,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "LvsL"), aes(x=generation,y=intercept,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "Within"), aes(x=generation,y=intercept,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "All"), aes(x=generation,y=intercept,color=.id))+
    geom_vline(xintercept = env.change, linetype = "longdash") +
    coord_cartesian(xlim = NULL,
                    ylim = NULL,
                    expand = T) + ggtitle(label = "IBD (Fst/Geographical distance) intercept over time", subtitle = name_of_simulation)+
    scale_color_manual(name = "Fragment\ncomparison",
                       values = forestComparisonPallete)

  p.slope.mod = ggplot(slope_intercept_df) +
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "SvsS"), aes(x=generation,y=slope_mod,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "SvsMS"), aes(x=generation,y=slope_mod,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "SvsML"), aes(x=generation,y=slope_mod,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "SvsL"), aes(x=generation,y=slope_mod,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "MSvsMS"), aes(x=generation,y=slope_mod,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "MSvsML"), aes(x=generation,y=slope_mod,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "MSvsL"), aes(x=generation,y=slope_mod,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "MLvsML"), aes(x=generation,y=slope_mod,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "MLvsL"), aes(x=generation,y=slope_mod,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "LvsL"), aes(x=generation,y=slope_mod,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "Within"), aes(x=generation,y=slope_mod,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "All"), aes(x=generation,y=slope_mod,color=.id))+
    geom_vline(xintercept = env.change, linetype = "longdash") +
    coord_cartesian(xlim = NULL,
                    ylim = NULL,
                    expand = T) + ggtitle(label = "IBD (Fst/Geographical distance) slope over time", subtitle = name_of_simulation)+
    scale_color_manual(name = "Fragment\ncomparison",
                       values = forestComparisonPallete)

  p.intercept.mod = ggplot(slope_intercept_df) +
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "SvsS"), aes(x=generation,y=intercept_mod,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "SvsMS"), aes(x=generation,y=intercept_mod,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "SvsML"), aes(x=generation,y=intercept_mod,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "SvsL"), aes(x=generation,y=intercept_mod,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "MSvsMS"), aes(x=generation,y=intercept_mod,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "MSvsML"), aes(x=generation,y=intercept_mod,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "MSvsL"), aes(x=generation,y=intercept_mod,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "MLvsML"), aes(x=generation,y=intercept_mod,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "MLvsL"), aes(x=generation,y=intercept_mod,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "LvsL"), aes(x=generation,y=intercept_mod,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "Within"), aes(x=generation,y=intercept_mod,color=.id))+
    geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "All"), aes(x=generation,y=intercept_mod,color=.id))+
    geom_vline(xintercept = env.change, linetype = "longdash") +
    coord_cartesian(xlim = NULL,
                    ylim = NULL,
                    expand = T) + ggtitle(label = "IBD (Fst/Geographical distance) intercept over time", subtitle = name_of_simulation)+
    scale_color_manual(name = "Fragment\ncomparison",
                       values = forestComparisonPallete)


  melt.df = melt(slope_intercept_df, id.vars=c(1,6))

  p.melt = ggplot(melt.df)+
    geom_point(aes(x=generation,y=value,color=.id,group=variable))+
    geom_vline(xintercept = env.change, linetype = "longdash") +
    coord_cartesian(xlim = NULL,
                    ylim = NULL,
                    expand = T) + ggtitle(label = "IBD coefficients over time", subtitle = name_of_simulation)+
    scale_color_manual(name = "Fragment\ncomparison",
                       values = forestComparisonPallete)+
    facet_grid(variable~.,scales="free")


  p.slope.overall = ggplot(slope_intercept_df) + geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "All"), aes(x=generation,y=slope,color=.id))+
    geom_vline(xintercept = env.change, linetype = "longdash") +
    coord_cartesian(xlim = NULL,
                    ylim = NULL,
                    expand = T) + ggtitle(label = "IBD slope over time", subtitle = name_of_simulation)+
    scale_color_manual(name = "Fragment\ncomparison",
                       values = forestComparisonPallete)

  p.intercept.overall = ggplot(slope_intercept_df) +  geom_point(data = subset(x = slope_intercept_df,slope_intercept_df$.id == "All"), aes(x=generation,y=intercept,color=.id))+
    geom_vline(xintercept = env.change, linetype = "longdash") +
    coord_cartesian(xlim = NULL,
                    ylim = NULL,
                    expand = T) + ggtitle(label = "IBD intercept over time", subtitle = name_of_simulation)+
    scale_color_manual(name = "Fragment\ncomparison",
                       values = forestComparisonPallete)

  if(FALSE){p.slope = ggplot(slope.df) +
    #geom_line(aes(x=slope.df$generation_list, y=slope.df$mean))+
    geom_point(aes(x = slope.df$generation_list, y = slope.df$mean)) +
    geom_vline(xintercept = env.change, linetype = "longdash") +
    #geom_abline(data = all_model_information[[2]],slope = all_model_information[[2]]$coefficients[[2]],intercept = all_model_information[[2]]$coefficients[[1]], color="red")+
    #geom_abline(data = all_model_information[[1201]],slope = all_model_information[[1201]]$coefficients[[2]],intercept = all_model_information[[1201]]$coefficients[[1]], color="blue")+
    #geom_smooth(aes(x = slope.df$generation_list, y = slope.df$mean)) +
    coord_cartesian(xlim = NULL,
                    ylim = NULL,
                    expand = T) + ggtitle(label = "IBD slope over time", subtitle = name_of_simulation)

  p.intercept = ggplot(intercept.df) +
    #geom_line(aes(x=intercept.df$generation_list, y=intercept.df$mean))+
    geom_point(aes(x = intercept.df$generation_list, y = intercept.df$mean)) +
    geom_vline(xintercept = env.change, linetype = "longdash") +
    #geom_abline(data = all_model_information[[2]],slope = all_model_information[[2]]$coefficients[[2]],intercept = all_model_information[[2]]$coefficients[[1]], color="red")+
    #geom_abline(data = all_model_information[[1201]],slope = all_model_information[[1201]]$coefficients[[2]],intercept = all_model_information[[1201]]$coefficients[[1]], color="blue")+
    #geom_smooth(aes(x = intercept.df$generation_list, y = intercept.df$mean)) +
    coord_cartesian(xlim = NULL,
                    ylim = NULL,
                    expand = T) + ggtitle(label = "IBD intercept over time", subtitle = name_of_simulation)
  }

  if (FALSE) {
    ggplot() + geom_abline(
      data = all_model_information[[2]],
      slope = all_model_information[[2]]$coefficients[[2]],
      intercept = all_model_information[[2]]$coefficients[[1]],
      color = "red"
    ) +
      geom_abline(
        data = all_model_information[[1197]],
        slope = all_model_information[[1197]]$coefficients[[2]],
        intercept = all_model_information[[1197]]$coefficients[[1]],
        color = "blue"
      ) +
      coord_cartesian(xlim = c(-1, 1),
                      ylim = c(-1, 1),
                      expand = F)
  }

  if (FALSE) {
    ## method to build geom_smooth of the data without actually using geom_smooth
    smooth_values = data.frame(predict(new_model, se = T)) ## builds smooth values along with standard error



    df <- data.frame(
      cbind(
        value = new_data_test$value,
        mean = new_data_test[[53]],
        fit = smooth_values$fit,
        upperBound = smooth_values$fit + 2 * smooth_values$se.fit,
        lowerBound = smooth_values$fit - 2 * smooth_values$se.fit
      )
    )

    g <- ggplot(df, aes(value, mean))
    #g <- g + geom_point()
    g <-
      g + geom_linerange(aes(ymin = lowerBound, ymax = upperBound))
    g <- g + geom_point(aes(value, fit))
    g <- g + geom_smooth(method = "lm",
                         color = "green",
                         size = 2)
    g <-
      g + geom_abline(
        slope = new_model$coefficients[[2]],
        intercept = new_model$coefficients[[1]],
        color = "red"
      )
    g
  }

  return(list(p.slope, p.intercept, p.slope.mod, p.intercept.mod, p.melt, p.slople.overall, p.intercept.overall))

}

Classify_Pop_Comparison = function(dist.df,#df with population names on the first 2 columns
                                   vector_small_frags = c(2, 4, 5, 7, 14, 15, 16, 17, 23, 26, 27, 28, 33, 34, 35, 41, 46, 47, 48, 49),
                                   vector_mediumSmall_frags = c(1, 3, 6, 8, 12, 13, 18, 19, 21, 22, 29, 30, 31, 32, 37, 39, 42, 43, 44, 45),
                                   vector_mediumLarge_frags = c(9, 10, 11, 24, 25, 36, 38, 40),
                                   vector_large_frags = c(20)
){

  popVS = NULL

  SvsS = vector_small_frags
  SvsMS = c(vector_small_frags,vector_mediumSmall_frags)
  SvsML = c(vector_small_frags,vector_mediumLarge_frags)
  SvsL = c(vector_small_frags,vector_large_frags)
  MSvsMS = vector_mediumSmall_frags
  MSvsML = c(vector_mediumSmall_frags,vector_mediumLarge_frags)
  MSvsL = c(vector_mediumSmall_frags,vector_large_frags)
  MLvsML = vector_mediumLarge_frags
  MLvsL = c(vector_mediumLarge_frags,vector_large_frags)
  LvsL = vector_large_frags


  for (i in 1:nrow(dist.df)) {

    if(as.numeric(strsplit(
      as.character(dist.df[i, 2]), split = "[aA-zZ]+"
    )[[1]][2]) == as.numeric(strsplit(
      as.character(dist.df[i, 1]), split = "[aA-zZ]+"
    )[[1]][2])){
      popVS = c(popVS, "Within")
    } else {

      if(any(as.numeric(strsplit(
        as.character(dist.df[i, 1]), split = "[aA-zZ]+"
      )[[1]][2]) == SvsS) &&
      any(as.numeric(strsplit(
        as.character(dist.df[i, 2]), split = "[aA-zZ]+"
      )[[1]][2]) == SvsS)){
        popVS = c(popVS, "SvsS")
      }else if(any(as.numeric(strsplit(
        as.character(dist.df[i, 1]), split = "[aA-zZ]+"
      )[[1]][2]) == MSvsMS) &&
      any(as.numeric(strsplit(
        as.character(dist.df[i, 2]), split = "[aA-zZ]+"
      )[[1]][2]) == MSvsMS)){
        popVS = c(popVS, "MSvsMS")
      }else if(any(as.numeric(strsplit(
        as.character(dist.df[i, 1]), split = "[aA-zZ]+"
      )[[1]][2]) == MLvsML) &&
      any(as.numeric(strsplit(
        as.character(dist.df[i, 2]), split = "[aA-zZ]+"
      )[[1]][2]) == MLvsML)){
        popVS = c(popVS, "MLvsML")
      }else if(any(as.numeric(strsplit(
        as.character(dist.df[i, 1]), split = "[aA-zZ]+"
      )[[1]][2]) == LvsL) &&
      any(as.numeric(strsplit(
        as.character(dist.df[i, 2]), split = "[aA-zZ]+"
      )[[1]][2]) == LvsL)){
        popVS = c(popVS, "LvsL")
      }else if(any(as.numeric(strsplit(
        as.character(dist.df[i, 1]), split = "[aA-zZ]+"
      )[[1]][2]) == SvsMS) &&
      any(as.numeric(strsplit(
        as.character(dist.df[i, 2]), split = "[aA-zZ]+"
      )[[1]][2]) == SvsMS)){
        popVS = c(popVS, "SvsMS")
      }else if(any(as.numeric(strsplit(
        as.character(dist.df[i, 1]), split = "[aA-zZ]+"
      )[[1]][2]) == SvsML) &&
      any(as.numeric(strsplit(
        as.character(dist.df[i, 2]), split = "[aA-zZ]+"
      )[[1]][2]) == SvsML)){
        popVS = c(popVS, "SvsML")
      }else if(any(as.numeric(strsplit(
        as.character(dist.df[i, 1]), split = "[aA-zZ]+"
      )[[1]][2]) == SvsL) &&
      any(as.numeric(strsplit(
        as.character(dist.df[i, 2]), split = "[aA-zZ]+"
      )[[1]][2]) == SvsL)){
        popVS = c(popVS, "SvsL")
      }else if(any(as.numeric(strsplit(
        as.character(dist.df[i, 1]), split = "[aA-zZ]+"
      )[[1]][2]) == MSvsML) &&
      any(as.numeric(strsplit(
        as.character(dist.df[i, 2]), split = "[aA-zZ]+"
      )[[1]][2]) == MSvsML)){
        popVS = c(popVS, "MSvsML")
      }else if(any(as.numeric(strsplit(
        as.character(dist.df[i, 1]), split = "[aA-zZ]+"
      )[[1]][2]) == MSvsL) &&
      any(as.numeric(strsplit(
        as.character(dist.df[i, 2]), split = "[aA-zZ]+"
      )[[1]][2]) == MSvsL)){
        popVS = c(popVS, "MSvsL")
      }else if(any(as.numeric(strsplit(
        as.character(dist.df[i, 1]), split = "[aA-zZ]+"
      )[[1]][2]) == MLvsL) &&
      any(as.numeric(strsplit(
        as.character(dist.df[i, 2]), split = "[aA-zZ]+"
      )[[1]][2]) == MLvsL)){
        popVS = c(popVS, "MLvsL")
      }

    }

  }
  dist.df$popVS = popVS

  return(dist.df)

}



Coefficients_per_comparison_class = function(
  dist_data_all_generations,
  generation_list
){



  # dist_data_all_generations = new_data_test
  list_slope_intercept_data = as.data.frame(matrix(NA,nrow=1,ncol=5))
  list_slope_intercept_all =  as.data.frame(matrix(NA,nrow=1,ncol=5))
  ##cols: slope intercept slope_mod intercept_mod generation
  list_coefficients_per_comparison_class = NULL

  bool_switch = T


  ## for each comparison between populations
  for(j in levels(factor(dist_data_all_generations$popVS))){

    for (i in 5:(length(dist_data_all_generations)-1)) {
      # exclude first four columns and last columns
      # firstGeneration + 1 <= i <= lastGeneration

      #compute linear regression genDist ~ geoDist
      if(bool_switch){ # switch because we only want to compute the stats for "all" once
        dummy.var.all.mod = lm(data = dist_data_all_generations,
                               formula = (dist_data_all_generations[, i]/(1-dist_data_all_generations[, i])) ~ log(dist_data_all_generations$value,exp(1)))
        dummy.var.all = lm(data = dist_data_all_generations,
                           formula = dist_data_all_generations[, i] ~ dist_data_all_generations$value)
        list_slope_intercept_all[i-4,] = rbind(dummy.var.all$coefficients[[2]], #slope
                                               dummy.var.all$coefficients[[1]], #intercept
                                               dummy.var.all.mod$coefficients[[2]], #slope_mod
                                               dummy.var.all.mod$coefficients[[1]], #intercept_mod
                                               generation_list[[i-3]])  #generation
      }
      dummy.var = lm(data = dist_data_all_generations,
                     formula = dist_data_all_generations[, i] ~ dist_data_all_generations$value,
                     subset = dist_data_all_generations$popVS == j)
      dummy.var.mod = lm(data = dist_data_all_generations,
                         formula = (dist_data_all_generations[, i]/(1-dist_data_all_generations[, i])) ~ log(dist_data_all_generations$value,exp(1)),
                         subset = dist_data_all_generations$popVS == j)
      list_slope_intercept_data[i-4,] = rbind(dummy.var$coefficients[[2]], #slope
                                              dummy.var$coefficients[[1]], #intercept
                                              dummy.var.mod$coefficients[[2]], #slope_mod
                                              dummy.var.mod$coefficients[[1]], #intercept_mod
                                              generation_list[[i-3]])  #generation
    }

    bool_switch=F
    names(list_slope_intercept_data) = c("slope","intercept","slope_mod","intercept_mod","generation")
    list_coefficients_per_comparison_class[[j]] = list_slope_intercept_data

  }

  names(list_slope_intercept_all) = c("slope","intercept","slope_mod","intercept_mod","generation")
  list_coefficients_per_comparison_class[["All"]] = list_slope_intercept_all

  return(list_coefficients_per_comparison_class)
}




################
#TEST AND PLOT IBD HYPOTHESIS
################
#+++++++++++++++
Test_IBD_hypothesis = function(genind_object) {
  # The original value of the correlation between the distance matrices is represented by the dot,
  # while histograms represent permuted values (i.e., under the absence of spatial structure).
  # Significant spatial structure would therefore result in the original value being out of the
  # reference distribution.
  Dgen = dist(genind_object)
  Dgeo = dist(genind_object$other$xy)
  ibd = mantel.randtest(Dgen, Dgeo)


  Plot_IBD_hypothesis = function(mantelRandTestResult) {
    plot(mantelRandTestResult)
  }
  Plot_IBD_hypothesis(ibd)
  return(ibd)
}

Test_IBD_hypothesis_ggplot = function(genind_object, distanceMethod = "manhattan") {
  # The original value of the correlation between the distance matrices is represented by the dot,
  # while histograms represent permuted values (i.e., under the absence of spatial structure).
  # Significant spatial structure would therefore result in the original value being out of the
  # reference distribution.

  genind_object = raw_data_multi_sim[[1]][[10]]
  genind_to_pop = genind2genpop(genind_object, process.other = T, quiet = T)
  genind_object@other$xy
  Dgen = dist(genind_to_pop)
  Dgeo = dist(genind_to_pop$other$xy, method = distanceMethod)
  df.Dgen = Dist.to.df(Dgen)
  df.Dgeo = Dist.to.df(Dgeo)
  ibd = mantel.randtest(Dgen, Dgeo)
  as.matrix(as.data.frame(lapply(df.Dgeo, as.numeric)))
  mantel.randtest(df.Dgen, df.Dgeo)
  ibd

  df.ibd = data.frame(ibd$sim)
  head(df.ibd)
  ggplot(data = df.ibd) + geom_histogram(
    aes(x = ibd.sim),
    boundary = 0,
    binwidth = 0.005,
    colour = "black",
    fill = "white"
  ) +
    geom_density(aes(x = ibd.sim), alpha = .2, fill = "#FF6666") +
    geom_vline(xintercept = ibd$obs, color = "red")

  Plot_IBD_hypothesis = function(mantelRandTestResult) {
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

summarySE <-
  function(data = NULL,
           measurevar,
           groupvars = NULL,
           na.rm = FALSE,
           conf.interval = .95,
           .drop = TRUE) {
    ## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
    ##   data: a data frame.
    ##   measurevar: the name of a column that contains the variable to be summariezed
    ##   groupvars: a vector containing names of columns that contain grouping variables
    ##   na.rm: a boolean that indicates whether to ignore NA's
    ##   conf.interval: the percent range of the confidence interval (default is 95%)
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm = FALSE) {
      if (na.rm)
        sum(!is.na(x))
      else
        length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(
      data,
      groupvars,
      .drop = .drop,
      .fun = function(xx, col) {
        c(
          N    = length2(xx[[col]], na.rm = na.rm),
          mean = mean   (xx[[col]], na.rm = na.rm),
          sd   = sd     (xx[[col]], na.rm = na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <-
      datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval / 2 + .5, datac$N - 1)
    datac$ci <- datac$se * ciMult

    return(datac)
  }

## make dist obj into dataframe (so that we can plot it)
Dist.to.df <- function(inDist) {
  if (class(inDist) != "dist")
    stop("wrong input type")
  A <- attr(inDist, "Size")
  B <-
    if (is.null(attr(inDist, "Labels")))
      sequence(A)
  else
    attr(inDist, "Labels")
  if (isTRUE(attr(inDist, "Diag")))
    attr(inDist, "Diag") <- FALSE
  if (isTRUE(attr(inDist, "Upper")))
    attr(inDist, "Upper") <- FALSE
  data.frame(row = B[unlist(lapply(sequence(A)[-1], function(x)
    x:A))],
    col = rep(B[-length(B)], (length(B) - 1):1),
    value = as.vector(inDist))
}


## save plots with set size
save_the_plot = function(plotObject, pathToPlot = "~/Documents/", nameOfPlotFile) {
  ggplot2::ggsave(
    plot = plotObject,
    file = paste(pathToPlot, "/", nameOfPlotFile, ".pdf", sep = ""),
    units = "cm",
    width = 15,
    height = 11,
    dpi = 600
  )

}


PopNameDf_to_PopSizeTypeDf = function(vectorPopNames,
                                      vector_small_frags,
                                      vector_mediumSmall_frags,
                                      vector_mediumLarge_frags,
                                      vector_large_frags){
  ## Converts a given vector of population names with the structure "pop[NUMBER]" and turns it into
  ## one of four given "classes" of populations from a type of forest fragment.
  ## The classes are Large, Small, MedSmall and MedLarge (-forest fragments)


  ##############
  #vectorPopNames = summaryPerPop$ind
  #vector_small_frags = c(2, 4, 5, 7, 14, 15, 16, 17, 23, 26, 27, 28, 33, 34, 35, 41, 46, 47, 48, 49)
  #vector_mediumSmall_frags = c(1, 3, 6, 8, 12, 13, 18, 19, 21, 22, 29, 30, 31, 32, 37, 39, 42, 43, 44, 45)
  #vector_mediumLarge_frags = c(9, 10, 11, 24, 25, 36, 38, 40)
  #vector_large_frags = c(20)
  ##############
  ##############


  popName_to_popSizeType = lapply(X = vectorPopNames, FUN = function(X){

    X[any(unlist(strsplit(X, split = "[aA-zZ]+"))[2]==as.numeric(vector_mediumSmall_frags))] = "MedSmall"
    X[any(unlist(strsplit(X, split = "[aA-zZ]+"))[2]==as.numeric(vector_small_frags))] = "Small"
    X[any(unlist(strsplit(X, split = "[aA-zZ]+"))[2]==as.numeric(vector_large_frags))] = "Large"
    X[any(unlist(strsplit(X, split = "[aA-zZ]+"))[2]==as.numeric(vector_mediumLarge_frags))] = "MedLarge"

    return(X)

  })

  popName_to_popSizeType = unlist(popName_to_popSizeType)
  return(popName_to_popSizeType)
}


###########################################################################
###########################################################################
#++++SECTION 4 - SNPs DATA ANALYSIS - ADD TO THE MAIN SCRIPT
###########################################################################
###########################################################################





Get_seqSNP_data = function(number_of_simulations,
                           path_to_data,
                           generations_list,
                           layer_ID
                           ){

  #number_of_simulations = 2
  #path_to_data = asd
  #generations_list = seq(10,50,10)
  #layer_ID="layer0"

  raw_snp_data = vector(mode = "list", length = number_of_simulations)

  for(simID in 1:number_of_simulations){

    setwd(paste(path_to_data, "/Adegenet_sim_", simID, "/", sep = ""))

    raw_data_per_sim = vector(mode = "list", length = length(generations_list))

      for(generationID in 1:length(generations_list)){

        currentGeneration = generations_list[generationID]

        name_of_file = paste(currentGeneration,layer_ID,"SNP",sep = "_")
        print(name_of_file)
        raw_data_per_sim[[generationID]] = read.PLINK(file = paste(name_of_file,".raw",sep = ""),
                   map.file = paste(name_of_file,".map",sep = ""),
                   saveNbAlleles = T,
                   quiet = T)

        # this block processes the xy coordinates from each individuals name
        # since they do not come explicitly in the data file
        xy_coords = NULL
        for(individual in raw_data_per_sim[[generationID]]@ind.names){
          xy_coords= rbind(xy_coords,
                           cbind(unlist(strsplit(individual,split = "_"))[4], #x position
                                 unlist(strsplit(individual,split = "_"))[5])) #y position
        }
        xy_coords = apply(xy_coords, 2, as.numeric) #xy coords as numerics instead of characters
        raw_data_per_sim[[generationID]]@other$xy = xy_coords

        # adding generation to genlight obj
        raw_data_per_sim[[generationID]]@other$generation = currentGeneration


      }

    raw_snp_data[[simID]] = raw_data_per_sim

  }

  return(raw_snp_data)

}

mean_geo_position_per_genlight_pop = function(raw_snp_data, simulationID, generationID){

  #raw_snp_data = test_snp_data
  #simulationID = 1
  #generationID = 5

  separated_pops = seppop(raw_snp_data[[simulationID]][[generationID]])

  pop_geo_xy = NULL

  for(pops in separated_pops){ # for each pop (genlight obj)

    pop_geo_xy = rbind(pop_geo_xy, pops=colMeans(pops@other$xy)) #mean x and y positions

  }

  row.names(pop_geo_xy) = popNames(raw_snp_data[[simulationID]][[generationID]])


  return(pop_geo_xy)

}

calc_dist_geo_genlight = function(mean_geo_position_per_genlight_pop, dist_method = "manhattan"){

  #calc dist obj directly from mean_geo_position_per_genlight_pop with the given method
  #and transform it into a dataframe
  return(Dist.to.df(dist(mean_geo_position_per_genlight_pop,method = dist_method)))

}



calc_mean_Fst = function(raw_snp_data){

  #raw_snp_data = test_snp_data

  fst_per_gen = NULL
  fst_per_sim = NULL
  fst_per_gen_names = NULL

  for(genIdx in 1:length(raw_snp_data[[1]])){
    for(simIdx in 1:length(raw_snp_data)){

      print(paste("Calculating Fst for Sim:",simIdx,
                  "generation:", raw_snp_data[[simIdx]][[genIdx]]@other$generation,
                  sep = " "))

      # stamppFst gets Fsts, pvalues, bootstrap values, CI values
      stamppFstPvalueObj = stamppFst(raw_snp_data[[simIdx]][[genIdx]])
      # currently only getting Fsts
      fst_per_sim[[simIdx]] = stamppFstPvalueObj$Fsts



    }

    fst_reduce = Reduce('+',fst_per_sim) # sum over all elements of list
    fst_per_gen[[genIdx]] = fst_reduce / length(fst_per_sim) #mean

    fst_per_gen[[genIdx]] = as.dist(fst_per_gen[[genIdx]])
    fst_per_gen[[genIdx]] = Dist.to.df(fst_per_gen[[genIdx]])

    #name each list obj after their respective generation
    fst_per_gen_names = c(fst_per_gen_names,raw_snp_data[[1]][[genIdx]]@other$generation)
    names(fst_per_gen) = fst_per_gen_names
  }

  return(fst_per_gen)

}


plot_IBD = function(mean_fst_per_gen_dflist, geo_dist_df, generation_to_plot){

  #generation_to_plot=50
  #mean_fst_per_gen_dflist=mean_fst_per_gen
  #geo_dist_df = test_snp_geo_data

  data_to_plot = mean_fst_per_gen_dflist[as.character(generation_to_plot)]
  data_to_plot = data_to_plot[[1]] #unlisting

  data_to_plot$geoDist = geo_dist_df$value

  ggplot(data = data_to_plot)+
    geom_point(aes(x=geoDist,y=value))+
    geom_smooth(aes(x=geoDist,y=value),method = "lm")


}




if(FALSE){
asd = "~/Documents/Documents/NetBeansProjects/SINS2_oct_20132/SINS2_oct_2013/SINStoArlequin/output/test_outputFormat/STANDARD_5x5_test"
#get seqSNP data
test_snp_data = Get_seqSNP_data(number_of_simulations = 2,path_to_data = asd,generations_list = seq(130,150,10),layer_ID = "layer0")
#get geo data per pop
test_snp_geo_data = calc_dist_geo_genlight(mean_geo_position_per_genlight_pop(test_snp_data,1,2))
#calc fst
mean_fst_per_gen = calc_mean_Fst(test_snp_data)

#plot IBD
plot_IBD(mean_fst_per_gen_dflist = mean_fst_per_gen,
         geo_dist_df = test_snp_geo_data,
         generation_to_plot = 150)






old_stuff_erase_after_revising = function(){



setwd("~/server_folderN/SINS_Sampler/dist/output/sim_SNPtest_10x10_k100_r04_m01_t100_seqSNP100len_NEWFORMATTEST/sim_SNPtest_10x10_k100_r04_m01_t100_seqSNP100len_NEWFORMATTEST/Adegenet_sim_1")
str(pli500)
pli500 = read.PLINK("500_layer0_AutosomeSNP.raw",map.file = "500_layer0_AutosomeSNP.map",saveNbAlleles=T)
glPlot(pli500)
pli50 = extract.PLINKmap("50_layer0_AutosomeSNP.map",x = pli50)
pli100 = read.PLINK("100_layer0_AutosomeSNP.raw")
pli100 = extract.PLINKmap("100_layer0_AutosomeSNP.map",x = pli100)
plot(pli50)

plot(pli100)

snpposi.plot(pli500@other$position, genome.size = 220,codon = F)
snpposi.plot(position(pli100), genome.size = 220,codon = F)
snpposi.test(position(pli100), genome.size = 220,codon = F)
locNames(pli500)
plot(hclust(dist(pli500)),labels=F)
pli500.dapc = dapc(pli500)
scatter(pli500.dapc)
pli500_pca = glPca(pli500)
scatter(pli500_pca,xlim =c(-6,6) ,ylim = c(-6,6))
myCol <- colorplot(pli100_pca$scores,pli100_pca$scores, transp=TRUE, cex=4)
add.scatter.eig(pli100_pca$eig[1:40],2,1,2, posi="topright", inset=.05, ratio=.3)


grp <- find.clusters(pli100)
plot(grp$Kstat, type="b", col="blue")


pli100@other$chromosome = NULL
chromosome(pli100)

pli100_pca = glPca(pli100)
scatter(pli100_pca,xax = 1,yax = 2,xlim =c(-6,6) ,ylim = c(-6,6))
myCol <- colorplot(pli100_pca$scores,pli100_pca$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")
add.scatter.eig(pli100_pca$eig[1:40],2,1,2, posi="topright", inset=.05, ratio=.3)
library(ape)
tre <- nj(dist(as.matrix(pli100)))
plot(tre, typ="fan", show.tip=FALSE)
tiplabels(pch=20, col=myCol, cex=4)


loadingplot.glPca(pli100_pca)
?loadingplot
loadingplot(as.matrix(pli100)[,2],srt = 70)

pli100_dapc = dapc(pli100)
scatter(pli100_dapc,scree.da=FALSE, bg="white", posi.pca="topright", legend=TRUE,
        txt.leg=paste("group", 1:2), col=c("red","blue"))
compoplot(pli100_dapc, col=c("red","blue"),lab="", txt.leg=paste("group", 1:2), ncol=2)


pli100_glMean <- glMean(pli100)
pli100_glMean <- c(pli100_glMean,1-pli100_glMean)
hist(pli100_glMean, proba=TRUE, col="gold", xlab="Allele frequencies",
     main="Distribution of (second) allele frequencies", nclass=20,xlim = c(0,1))
temp <- density(pli100_glMean,bw=.05)
lines(temp$x, temp$y*1.8,lwd=3)

pli100_glDot = glDotProd(pli100)


lD <- lapply(pli100, function(e) dist(as.matrix(e)))
D <- Reduce("+",lD)
library(ape)
plot(nj(D), type="fan")


library(adegenet)
setwd("~/Documents/Documents/NetBeansProjects/SINS2_oct_20132/SINS2_oct_2013/SINStoArlequin/output/test_outputFormat/STANDARD_5x5_test/Adegenet_sim_1/")
pli10 = read.PLINK("50_layer0_AutosomeSNP.raw")
pli10WithMap = extract.PLINKmap("50_layer0_AutosomeSNP.map",x = pli10)
as.matrix(pli10WithMap)

pli10WithMap@other$chromosome
plot(pli10WithMap)
snpposi.plot(pli10WithMap@other$position, genome.size = 220,codon = F)
glMean(pli10WithMap)

#dist

plot(hclust(dist(pli10WithMap)),labels=F)
summary(dapc(pli10WithMap))
summary()

test.a = as.matrix(pli10)
test.a=test.a[,-1]
colnames(test.a) = c(1:10)
test.a.g = df2genind(test.a, ploidy=2,sep = " ")

?extract.PLINKmap




myPath <- system.file("files/usflu.fasta",package="adegenet")
flu <- fasta2genlight(myPath, chunk=10, parallel=FALSE)
flu@gen[[10]]@ploidy


test.fasta = fasta2genlight("test.fasta")
as.matrix(test.fasta)





a=genind(as.matrix(pli50))

a=glPca(pli50,nf = 2)
scatter(a)

ab=dapc(pli500,n.da = 2,n.pca = 10)
scatter(ab,scree.da=FALSE, bg="white", posi.pca="topright", legend=TRUE,
        txt.leg=paste("group", 1:2), col=c("red","blue"))
compoplot(ab, col=c("red","blue"),lab="", txt.leg=paste("group", 1:2), ncol=2)
  }

}
