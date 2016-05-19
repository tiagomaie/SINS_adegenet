library("adegenet")
library("ade4")
library(hierfstat)
library(pegas)
library(ggplot2)

setwd("/home/tiago/Documents/SINS_VALIDATION_AND_TESTING/SINS_Pipeline/SINS/dist/results/SINS_NatalieScenario_Simplified_4/simulation_1/")


#################################################################################################
#################################################################################################
#################################################################################################
######################################### HAPLOID EXAMPLE #######################################
#################################################################################################
#################################################################################################
#################################################################################################



#read first 5 rows of the table
tab5rows <- read.table("layerOne_Y.txt", header = FALSE, nrows = 5)
#get the classes of the columns of the table
classes <- sapply(tab5rows, class)
#read the rest (or a bigger part) of the table knowing the class of the columns speeds up the reading process



if(FALSE){
  Yzome_10000rows = NULL
  Yzome_10000rows_simple=NULL
  Yzome_10000rows <- read.table("layerOne_Y.txt", header = FALSE, nrows=10000, colClasses = classes)
  
  #convert factor to string. dont do it before because apparently its faster to read the table with the right classes
  #already known. might have to test this
  factor_cols=sapply(Yzome_10000rows, is.factor) #get boolean list of cols that are factors
  Yzome_10000rows[factor_cols]=sapply(Yzome_10000rows[factor_cols], as.character) #convert factor cols to char cols
  
  
  x_coords= NULL
  y_coords= NULL
  for(i in 1:length(Yzome_10000rows[,4])){
    
    x_coords = c(x_coords, strsplit(Yzome_10000rows[i,4],split="_")[[1]][4])
    y_coords = c(y_coords, strsplit(Yzome_10000rows[i,4],split="_")[[1]][5])
  }
  x_coords=as.list(as.integer(x_coords))
  y_coords=as.list(as.integer(y_coords))
  
  
  Yzome_import_cols = c("V7")
  Yzome_10000rows_simple = Yzome_10000rows[Yzome_import_cols]

  #Yzome_10000rows[,4]=as.character(Yzome_10000rows[,4])
  for(i in 1:length(Yzome_10000rows[,4])){
    
    Yzome_10000rows[,4][i] = paste(i,"_",Yzome_10000rows[,4][i],sep="")
    
  }
  

  
  the_coords=NULL
  the_coords = cbind(x_coords,y_coords)
  str(as.list(the_coords))
  head(the_coords)

  length(x_coords)
  length(y_coords)
  
  row.names(Yzome_10000rows_simple) = Yzome_10000rows[,4]  
  rownames(Yzome_10000rows_simple)
  
  str(Yzome_10000rows_simple)
  Yzome_10000rows_simple
  
  names(Yzome_10000rows_simple) = c("genotype")
  head(Yzome_10000rows_simple)
  
  Yzome_genind10000=NULL
  Yzome_genind10000 = df2genind(X = Yzome_10000rows_simple, ploidy = 1, ncode = 3)
  
  #the right way to add stuff to @other
  Yzome_genind10000$other$xy = the_coords
  
  ##########################^^#############################
  
  Yzome_genind10000@tab[500:1000,]
  
  head(Yzome_genind10000@tab)
  head(Yzome_genind10000[1:5]@other$xy)

  fstat(Yzome_genind10000)
  
  
  
  
}else{
  #Yzome_10rows <- read.table("layerOne_Y.txt", header = FALSE, nrows=10, colClasses = classes)
  #Yzome_import_cols = c("V7")
  #Yzome_10rows_simple= Yzome_10rows[Yzome_import_cols]  
  #row.names(Yzome_10rows_simple) = Yzome_10rows[,4]
  #Yzome_genind10 = df2genind(X = Yzome_10rows_simple, ploidy = 1, ncode = 3)
}

#################################################################################################
#################################################################################################
#################################################################################################
####################################### HAPLOID EXAMPLE END #####################################
#################################################################################################
#################################################################################################
#################################################################################################



#################################################################################################
#################################################################################################
#################################################################################################
######################################### BIPLOID EXAMPLE #######################################
#################################################################################################
#################################################################################################
#################################################################################################


setwd("/home/tiago/Documents/Documents/NetBeansProjects/SINS2_oct_20132/SINS2_oct_2013/SINStoArlequin/output/test_outputFormat/STANDARD_5x5_test/Adegenet_sim1/")

Read_data_func = function(name_of_file, markers_list){
  
  #A1_marker_data = NULL
  #A1_marker_simple=NULL
  #col_names = c("ID","A1a","A1b","x","y","pop") defunct
  #A1_marker_data <- read.table("90_layer0_A1", header = FALSE, colClasses = classes, row.names = 1, col.names = col_names) defunct
  
  #A1_marker_data <- read.table("90_layer0_A1", header = FALSE, colClasses = classes, row.names = 1)
  #A1_marker_data <- read.table(name_of_file, header = FALSE, colClasses = classes, row.names = 1)
  #head(A1_marker_data)
  current_generation=NULL
  current_generation =strsplit(x=name_of_file,split = "_")[[1]][1]

  data_file = NULL
  marker_data = NULL
  bind_markers=NULL
  a_file_name = NULL
  for(marker in markers_list){
    a_file_name =paste(name_of_file,marker,sep = "")
    #read first 5 rows of the table
    #tab5rows <- read.table("90_layer0_A1", header = FALSE, nrows = 5,row.names =1 )
    tab5rows <- read.table(file = a_file_name, header = FALSE, nrows = 5,row.names =1)
    
    #get the classes of the columns of the table
    classes <- sapply(tab5rows, class)
    #read the rest (or a bigger part) of the table knowing the class of the columns speeds up the reading process
    
      
    data_file <- read.table(a_file_name, header = FALSE, colClasses = classes, row.names = 1)
    marker_data = as.matrix(paste(data_file$V2,data_file$V3,sep = "/"))
    bind_markers = cbind(bind_markers,marker_data)
      
  }
  markers_df2genind = df2genind(X = bind_markers, ploidy = 2, ncode = 3,ind.names = rownames(data_file),sep = "/")
  markers_df2genind@pop = data_file$V6
  xy_coords = cbind(data_file$V4,data_file$V5)
  markers_df2genind@other$xy = xy_coords
  markers_df2genind@other$generation = current_generation

  if(FALSE){
    #allele_cols = c("A1a","A1b")
    #A1_marker_simple = A1_marker_data[allele_cols]
    
    A1_marker_simple = as.matrix(paste(A1_marker_data$V2,A1_marker_data$V3,sep = "/"))
    #head(A1_marker_simple)
    
    A1df2genind = df2genind(X = A1_marker_simple, ploidy = 2, ncode = 3,ind.names = rownames(A1_marker_data),sep = "/")
    #head(A1df2genind@tab)
    
    
    #A1df2genind
    A1df2genind@pop = A1_marker_data$V6
    xy_coords = cbind(A1_marker_data$V4,A1_marker_data$V5)
    A1df2genind@other$xy = xy_coords
    #head(A1df2genind@tab)
    #head(A1df2genind@pop)
  }
  return(markers_df2genind)
}


Set_markers_list = function(...){
  the_markers_list=NULL
  the_markers_list = c(...)
  return(the_markers_list)
}
  
Set_generations_list = function(initialTimeStep, finalTimeStep, timeStepIncrement){
  the_generations_list=NULL
  the_generations_list = seq(initialTimeStep,finalTimeStep,timeStepIncrement)
  return(the_generations_list)
}


Set_layers_list = function(...){
  the_layers_list=NULL
  the_layers_list = c(...)
  return(the_layers_list)
}
  
Get_file_names = function(generations_list, layers_list){
  the_file_names = NULL

  for(j in generations_list){
   
    the_file_names = rbind(the_file_names, paste(j,"_",layers_list,"_",sep = ""))
     
  }
  return(the_file_names)
}

Get_data_from_files = function(the_file_names, set_markers,generations_list){
  the_data_list=NULL
  for(i in the_file_names){
    the_data_list = c(the_data_list,Read_data_func(i, set_markers))
  }
  #set names of elements in list to the generations that they belong to
  names(the_data_list) = generations_list
  return(the_data_list)
}


####################################################################################################
####################################################################################################
####################################################################################################
#Using summaries
####################################################################################################
####################################################################################################
####################################################################################################

#layer0_gen80 = Read_data_func("80_layer0_",set_markers)
#str(layer0_gen80)


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
  
  
  bartlett_test = bartlett.test(list(summary_genind$Hexp,summary_genind$Hobs))
  
  the_ttest = t.test(summary_genind$Hexp,summary_genind$Hobs,pair=T,var.equal=TRUE,alter="greater")
  
  results = list(summary_genind,bartlett_test,the_ttest)

  
  #return(summary_genind)
  return(results)
}


##TODO: attempt to show summaries with ggplot
Do_summary_calcs_ggplot = function(the_genind_obj){
  the_genind_obj=layer0_gen80
  summary_genind = summary(the_genind_obj)
  names(summary_genind)
  
  df_summary = data.frame(n=summary_genind$n,
                          n.by.pop=summary_genind$n.by.pop,
                          #summary_genind$loc.n.all,
                          pop.n.all=summary_genind$pop.n.all#,
                          #summary_genind$NA.perc,
                          #summary_genind$Hobs,
                          #summary_genind$Hexp
                          )
  
  
  par(mfrow=c(2,2))
  
  ggplot()+geom_jitter(aes(df_summary$n.by.pop, df_summary$pop.n.all), xlab="Colonies sample size",
       ylab="Number of alleles",main="Alleles numbers and sample sizes",
       type="n")
  text(summary_genind$n.by.pop,summary_genind$pop.n.all,lab=names(summary_genind$n.by.pop))
  barplot(summary_genind$loc.n.all, ylab="Number of alleles",
          main="Number of alleles per locus")
  barplot(summary_genind$Hexp-summary_genind$Hobs, main="Heterozygosity: expected-observed",
          ylab="Hexp - Hobs")
  barplot(summary_genind$n.by.pop, main="Sample sizes per population",
          ylab="Number of genotypes",las=3)
  
  bartlett.test(list(summary_genind$Hexp,summary_genind$Hobs))
  
  t.test(summary_genind$Hexp,summary_genind$Hobs,pair=T,var.equal=TRUE,alter="greater")
  
  
  return(summary_genind)
  
}


if(FALSE){
summary_calc_list = list()
pair_dist_plots=list()
a_counter = 1

for(i in the_data_list[]){
  
  
  if(a_counter>1){
    pair_dist_plots[[a_counter]] = pairDistPlot(i,boxplot = FALSE,data=TRUE,jitter = TRUE,violin = FALSE,within = TRUE) 
  summary_calc_list[[a_counter]] = do_summary_calcs(i)
  
  }
  #print(a_counter)
  a_counter = a_counter+1
  
}

pairDistPlot(the_data_list[[1]],grp =the_data_list[[1]]$pop, boxplot = TRUE,data=TRUE,jitter = TRUE,violin = TRUE,within = TRUE)



propShared(obj = the_data_list[[10]])
the_data_list[[2]]@tab
summary_calc_list[[2]]
pair_dist_plots[[2]]$data

#set names of elements in list to the generations that they belong to
names(summary_calc_list) = set_generations

as.data.frame(summary_calc_list[[5]])
a =data.frame(n.by.pop = summary_calc_list[[5]]$n.by.pop, pop.n.all=summary_calc_list[[5]]$pop.n.all)
a[[2]]
str(summary_calc_list[[2]])


pairDistPlot(layer0_gen80)

ggplot()+geom_jitter(data = a,aes(x=npop, y=popall))
+geom_jitter(aes(color=pop))


placeholder_var = do_summary_calcs(layer0_gen80)
test_var = NULL
test_var = list()
test_var$b = placeholder_var
test_var = list(test_var, do_summary_calcs(layer0_gen80))
str(test_var)
}

#DEFUNCT
if(FALSE){
summary_A1=summary(A1df2genind)
summary_A1$pop.n.all
summary_A1$n.by.pop
names(summary_A1)

A1genind2genpop = genind2genpop(A1df2genind,process.other = TRUE)
summary_A1pop = summary(A1genind2genpop)
names(summary_A1pop)

par(mfrow=c(2,2))
plot(x = summary_A1$n.by.pop, y = summary_A1$pop.n.all, xlab="Colonies sample size",
     ylab="Number of alleles",main="Alleles numbers and sample sizes",
     type="n")
text(summary_A1$n.by.pop,summary_A1$pop.n.all,lab=names(summary_A1$n.by.pop))
barplot(summary_A1$loc.n.all, ylab="Number of alleles",
        main="Number of alleles per locus")
barplot(summary_A1$Hexp-summary_A1$Hobs, main="Heterozygosity: expected-observed",
        ylab="Hexp - Hobs")
barplot(summary_A1$n.by.pop, main="Sample sizes per population",
        ylab="Number of genotypes",las=3)
}
#/DEFUNCT

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

fstat(layer0_gen80)
Fst(as.loci(layer0_gen80))

hw.test(layer0_gen80, B=0)

Gtest <- gstat.randtest(layer0_gen80,nsim=500)
plot(Gtest)

do_pairwise_fst = pairwise.fst(layer0_gen80)
is.euclid(do_pairwise_fst)

  (layer0_gen80)

#####################
#####################
#####################
pca1 <- dudi.pca(A1df2genind,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

s.class(pca1$li, pop(A1df2genind))
title("PCA of A1df2genind dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

the_col <- funky(15)
s.class(pca1$li, pop(A1df2genind),xax=1,yax=3, col=transp(the_col,.6), axesell=FALSE,
        cstar=0, cpoint=3, grid=FALSE)
#####################
#####################
#####################




A1genind2genpop = genind2genpop(A1df2genind)
Dgen <- dist.genpop(A1genind2genpop,method=2)
Dgeo <- dist(A1genind2genpop$other$xy)
ibd <- mantel.randtest(Dgen,Dgeo)
ibd
plot(ibd)


Dgen <- dist(A1df2genind$tab)
Dgeo <- dist(other(A1df2genind)$xy)
ibd <- mantel.randtest(Dgen,Dgeo)
ibd
plot(ibd)


plot(Dgeo, Dgen)

library(MASS)
dens <- kde2d(Dgeo,Dgen, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo, Dgen, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(Dgen~Dgeo))
title("Isolation by distance plot")


summary(A1df2genind$pop)

levels(A1df2genind$pop)

temp <- A1df2genind$pop
levels(temp) <- c(1,2,3,4)
temp <- as.numeric(as.character(temp))
plot(A1df2genind$other$xy, pch=temp, cex=1.5, xlab='x', ylab='y')
legend("bottom",leg=c("Pop A", "Pop B","Pop C","Pop D"),pch=c(1,2,3,4))


A1pairfst = pairwise.fst(A1df2genind)
is.euclid(A1pairfst)
D <- dist(A1df2genind$tab)
gab <- chooseCN(A1df2genind$other$xy,ask=FALSE,type=2)





#convert factor to string. dont do it before because apparently its faster to read the table with the right classes
#already known. might have to test this
factor_cols=sapply(A1zome_10000rows, is.factor) #get boolean list of cols that are factors
A1zome_10000rows[factor_cols]=sapply(A1zome_10000rows[factor_cols], as.character) #convert factor cols to char cols


x_coords= NULL
y_coords= NULL
for(i in 1:length(A1zome_10000rows[,4])){
  
  x_coords = c(x_coords, strsplit(A1zome_10000rows[i,4],split="_")[[1]][4])
  y_coords = c(y_coords, strsplit(A1zome_10000rows[i,4],split="_")[[1]][5])
}
x_coords=as.list(as.integer(x_coords))
y_coords=as.list(as.integer(y_coords))


A1zome_import_cols = c("V7","V8")
A1zome_10000rows_simple = A1zome_10000rows[A1zome_import_cols]

#Yzome_10000rows[,4]=as.character(Yzome_10000rows[,4])
for(i in 1:length(A1zome_10000rows[,4])){
  
  A1zome_10000rows[,4][i] = paste(i,"_",A1zome_10000rows[,4][i],sep="")
  
}



the_coords=NULL
the_coords = cbind(x_coords,y_coords)
str(as.list(the_coords))
head(the_coords)

length(x_coords)
length(y_coords)

row.names(A1zome_10000rows_simple) = A1zome_10000rows[,4]  
rownames(A1zome_10000rows_simple)

str(A1zome_10000rows_simple)
A1zome_10000rows_simple

names(A1zome_10000rows_simple) = c("genotype_1","genotype_2")
head(A1zome_10000rows_simple)
tail(A1zome_10000rows_simple)
A1zome_genind10000=NULL
A1zome_genind10000 = df2genind(X = A1zome_10000rows_simple, ploidy = 2, ncode = 3)
tail(A1zome_genind10000@tab)
#the right way to add stuff to @other
A1zome_genind10000$other$xy = the_coords



##defining populations##
the_pop = NULL
for (i in length(the_coords[,1])){
  
  if(the_coords[i,]$x_coords == 9 && the_coords[i,]$y_coords == 9)
    the_pop = c(the_pop,"cluster_3")
  else if(the_coords[i,]$x_coords == 7 && the_coords[i,]$y_coords == 7)
    the_pop = c(the_pop,"cluster_1")
  else if(the_coords[i,]$x_coords == 7 && the_coords[i,]$y_coords == 7)
    the_pop = c(the_pop,"cluster_1")
    
    
}



##########################^^#############################

A1zome_genind10000@tab[500:1000,]

head(A1zome_genind10000@tab)
head(A1zome_genind10000[1:5]@other$xy)

fstat(A1zome_genind10000)










#################################################################################################
#################################################################################################
#################################################################################################
########################################## PLINK TESTS ##########################################
#################################################################################################
#################################################################################################
#################################################################################################

df10 = genind2df(Yzome_genind10, sep="|")
df10
#TEST
some_obj = read.PLINK(file = "/home/tiago/PycharmProjects/ToTestOnly/test_input.raw")
some_obj
?read.plink
as.matrix(some_obj)
some_obj_mydata = read.PLINK(file="/home/tiago/PycharmProjects/ToTestOnly/my_data_test_input.raw")
as.matrix(some_obj_mydata)

glMean(some_obj_mydata, FALSE)
glMean(some_obj_mydata)
glSum(some_obj_mydata)
glSum(some_obj_mydata, FALSE)
summary(some_obj_mydata)
glPlot(some_obj_mydata)
glPca(some_obj_mydata)
glDotProd(center = TRUE,some_obj_mydata)
?glDotProd
some_obj_mydata@gen

as.matrix(some_obj_mydata)
ploidy(some_obj_mydata)
temp = as.matrix(some_obj_mydata)/ploidy(some_obj_mydata)
apply(temp,2,mean, na.rm=TRUE)

myFreq <- glMean(some_obj_mydata)
hist(myFreq, proba=TRUE, col="gold", xlab="Allele frequencies",
     main="Distribution of (second) allele frequencies")
temp <- density(myFreq)
lines(temp$x, temp$y*1.8,lwd=3)

some_obj_mydata
?"genlight"
some_obj_mydata@gen[[10]]@snp
some_obj_mydata@gen[[10]]@snp[[1]]
summary(some_obj_mydata@gen[[10]]@snp[[1]])
some_obj_mydata$gen
some_obj_mydata[1]
position(some_obj_mydata[1])


# test


test_snp_2_1 = new("SNPbin", c(1,1,1,2,1,1,2,1))
test_snp_2_2 = new("SNPbin", c(1,1,1,2,1,2,2,1))

test_genlight = new("genlight", list(test_snp_2_2, test_snp_2_1))

temp = as.matrix(test_genlight)/ploidy(test_genlight)
apply(temp,2,mean, na.rm=TRUE)

myFreq <- glMean(test_genlight)
hist(myFreq, proba=TRUE, col="gold", xlab="Allele frequencies",
     main="Distribution of (second) allele frequencies")
temp <- density(myFreq)
lines(temp$x, temp$y*1.8,lwd=3)
