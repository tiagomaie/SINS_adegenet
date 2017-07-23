### Analysis run for SINS paper

args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("At least one argument must be supplied (name of simulation).",
       call. = FALSE)
}


name_of_simulation = args[1]
currentDirectory = args[2]
dataDirectory = args[3]

initalGen = as.integer(args[4])
lastGen = as.integer(args[5])
intervalGen = as.integer(args[6])

setwd(currentDirectory)

if (FALSE) {
  name_of_simulation = "test_case_5x5_t100"
  name_of_simulation = "r05_m01_k100_10x10_t6000_f3000"
  name_of_simulation = "r05_m01_k100_21x_t25k_f20k_tw19-25_loTM"
  name_of_simulation = "r05_m01_k100_20x_t25k_f20k_tw19-25"
  name_of_simulation = "STANDARD_5x5_test"
  currentDirectory = "/home/tmaie/R_Projects/SINS_adegenet"
  dataDirectory = "~/ServerFolders/Elephant_ServerFolder/SINS_Sampler/dist/output"
  dataDirectory = ""
  initalGen = 100
  lastGen = 0
  intervalGen = 10


  setwd(currentDirectory)

}

source("SINS_data_analysis.r")


#++++++++++++++++++++++
######## INPUT ########
#++++++++++++++++++++++

# set parameters
the_number_of_simulations = Set_number_of_simulations(5)
the_generation_list = Set_generations_list(initalGen, lastGen, intervalGen)
the_layers_list = Set_layers_list("layer0")
the_markers_list = Set_markers_list( "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10")

the_path_to_data = paste(dataDirectory,
                         "/",
                         name_of_simulation,
                         "/",
                         name_of_simulation,
                         "/",
                         sep = "")
the_path_to_plot_folder = file.path("~", "R_Projects","R_Plots","SINS_R_Plots", name_of_simulation)

# Create dir for plots
dir.create(the_path_to_plot_folder, recursive = T)
#setwd(the_path_to_plot_folder)


#+++++++++++

# get names of the files and then get data from these files
the_file_names = Get_file_names(seq(10,100,10), the_layers_list)

raw_data_multi_sim = Get_data_from_files_multi_sim(
  the_file_names = the_file_names,
  set_markers = the_markers_list,
  generations_list = seq(10,100,10),
  path_to_data = the_path_to_data,
  number_of_simulations = the_number_of_simulations,
  is_marker_diploid = TRUE
)




#ibdQuemere=Make_pairFst_distance_plot_QUEMERE_Paper()
#save_the_plot(plotObject = ibdQuemere,pathToPlot = the_path_to_plot_folder,nameOfPlotFile = "IBDplot_Quemere_corrected")

if (is.element("t25k", unlist(strsplit(
  name_of_simulation, split = "_", fixed = T
)))) {
  print("IBD over time")
  IBD.over.time = Compute.PairwiseFst.IBD.Slope.OverTime(array_of_genind_array = raw_data_multi_sim,
                                                         name_of_simulation = name_of_simulation,
                                                         env.change = 20000)

  save_the_plot(
    plotObject = IBD.over.time[[1]],
    pathToPlot = the_path_to_plot_folder,
    nameOfPlotFile = paste(name_of_simulation, "_IBDOverTime_slope",
                           sep = "")
  )

  save_the_plot(
    plotObject = IBD.over.time[[2]],
    pathToPlot = the_path_to_plot_folder,
    nameOfPlotFile = paste(name_of_simulation, "_IBDOverTime_intercept",
                           sep = "")
  )



  if (is.element("loTM", unlist(strsplit(
    name_of_simulation, split = "_", fixed = T
  )))) {

    print("Compute.Hs.over.all.sims.per.locus")
    hs.allsim.per.locus = Compute.Hs.over.all.sims.per.locus(raw_data_multi_sim)
    hs_mean_plot = Make.HS.mean.plots(
      compute.Hs.all.sims = hs.allsim.per.locus,
      initial.gen = 19000,
      env.change = 20000,
      last.gen = 25000,
      subtitle = name_of_simulation,
      number_of_simulations = the_number_of_simulations,
      showPerDeme = T,
      vector_small_frags = c(2, 4, 5, 7, 14, 15, 16, 17, 23, 26, 27, 28, 33, 34, 35, 41, 46, 47, 48, 49),
      vector_mediumSmall_frags = c(1, 3, 6, 8, 12, 13, 18, 19, 21, 22, 29, 30, 31, 32, 37, 39, 42, 43, 44, 45),
      vector_mediumLarge_frags = c(9, 10, 11, 24, 25, 36, 38, 40),
      vector_large_frags = c(20)
    )
    save_the_plot(
      plotObject = hs_mean_plot,
      pathToPlot = the_path_to_plot_folder,
      nameOfPlotFile = paste(name_of_simulation, "_HSmeanPlot",
                             sep = "")
    )



  print("IBD plots")
  ibd.plots.list = lapply(X = seq(19000, 25000, 1000), function(X) {
    pairFst_dist_plots_allSims_fourFragClasses(
      array_of_genind_array = raw_data_multi_sim,
      generation = X,
      name_of_simulation = name_of_simulation
    )#[[2]]
  })

  }else{


    print("Compute.Hs.over.all.sims.per.locus")
    hs.allsim.per.locus = Compute.Hs.over.all.sims.per.locus(raw_data_multi_sim)
    hs_mean_plot = Make.HS.mean.plots(
      compute.Hs.all.sims = hs.allsim.per.locus,
      initial.gen = 19000,
      env.change = 20000,
      last.gen = 25000,
      subtitle = name_of_simulation,
      number_of_simulations = the_number_of_simulations,
      showPerDeme = T,
      vector_small_frags = c(1,2),
      vector_mediumSmall_frags = c(3,4,5),
      vector_mediumLarge_frags = c(6,7),
      vector_large_frags = c(8,9)
    )
    save_the_plot(
      plotObject = hs_mean_plot,
      pathToPlot = the_path_to_plot_folder,
      nameOfPlotFile = paste(name_of_simulation, "_HSmeanPlot",
                             sep = "")
    )


    print("IBD plots")
    ibd.plots.list = lapply(X = seq(19000, 25000, 1000), function(X) {
      pairFst_dist_plots_allSims_fourFragClasses(
        array_of_genind_array = raw_data_multi_sim,
        generation = X,
        name_of_simulation = name_of_simulation,
        vector_small_frags = c(1,2),
        vector_mediumSmall_frags = c(3,4,5),
        vector_mediumLarge_frags = c(6,7),
        vector_large_frags = c(8,9)
      )#[[2]]
    })


  }

  plotFourFragClassesIBD_toFile(ibd.plots.list,
                                the_path_to_plot_folder,
                                name_of_simulation)

} else{

  print("IBD over time")
  IBD.over.time = Compute.PairwiseFst.IBD.Slope.OverTime(array_of_genind_array = raw_data_multi_sim,
                                                         name_of_simulation = name_of_simulation,
                                                         env.change = 3000)
  save_the_plot(
    plotObject = IBD.over.time[[1]],
    pathToPlot = the_path_to_plot_folder,
    nameOfPlotFile = paste(name_of_simulation, "_IBDOverTime_slope",
                           sep = "")
  )

  save_the_plot(
    plotObject = IBD.over.time[[2]],
    pathToPlot = the_path_to_plot_folder,
    nameOfPlotFile = paste(name_of_simulation, "_IBDOverTime_intercept",
                           sep = "")
  )

  print("Compute.Hs.over.all.sims.per.locus")
  hs.allsim.per.locus = Compute.Hs.over.all.sims.per.locus(raw_data_multi_sim)
  hs_mean_plot = Make.HS.mean.plots(
    compute.Hs.all.sims = hs.allsim.per.locus,
    initial.gen = 0,
    env.change = 3000,
    last.gen = 6000,
    subtitle = name_of_simulation,
    number_of_simulations = the_number_of_simulations,
    showPerDeme = T,
    vector_small_frags = c(1,2),
    vector_mediumSmall_frags = c(3,4),
    vector_mediumLarge_frags = c(5),
    vector_large_frags = c(6)
  )

  save_the_plot(
    plotObject = hs_mean_plot,
    pathToPlot = the_path_to_plot_folder,
    nameOfPlotFile = paste(name_of_simulation, "_HSmeanPlot",
                           sep = "")
  )


  print("IBD plots")
  ibd.plots.list = lapply(X = seq(1000, 6000, 1000), function(X) {
    pairFst_dist_plots_allSims_fourFragClasses(
      array_of_genind_array = raw_data_multi_sim,
      generation = X,
      name_of_simulation = name_of_simulation,
      vector_small_frags = c(1,2),
      vector_mediumSmall_frags = c(3,4),
      vector_mediumLarge_frags = c(5),
      vector_large_frags = c(6)
    )#[[2]]
  })


  plotFourFragClassesIBD_toFile(ibd.plots.list,
                                the_path_to_plot_folder,
                                name_of_simulation)


}
