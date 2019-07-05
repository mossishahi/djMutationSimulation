## Here we defined a function, generate_scaffolds, which gets a fragmented draft genome and 
# and generates scaffolds wih mean.size whihch you set! 
# the length of scaffolds is variates around the mean size and by varicance factor you declare
#that how far you want to be different from the mean size of average in scaffolds size


# onelinedata <- readLines("./Data/OneLineString.txt")






base.genome <- read.delim("./Data/simulated_scaffolds.txt", sep = "\n", header = FALSE)

## Here we Mutate the draft genome declared as the input of this function and
## generate n different draft genoms which are Mutated against input draft
#finally we write each draft genome in different txt file
mutate_draft <- function(base.genome, mutation.rate, version.number){
  mutation.count <- mutation.rate * nchar(base.genome)
  base.genome.scaffolds <- as.character(base.genome$V1)
  base.genome.scaffolds.headers <- base.genome.scaffolds[grep(">",base.genome.scaffolds)]
  base.genome.scaffolds.sequences <- base.genome.scaffolds[-grep(">",base.genome.scaffolds)]
  base.genome.scaffolds.length <- nchar(base.genome.scaffolds.sequences)
  print("--------")
  print(base.genome.scaffolds.length)
  contigious.draft <- paste0(base.genome.scaffolds.sequences,collapse = '')
  # print(contigious.draft)
  # print(nchar(contigious.draft))
  mutation.indices <- sample(1:nchar(contigious.draft), mutation.count)
  print("These indices are chosen randomly to be Mutated")
  print(mutation.indices)
  
  for (i in 1:mutation.count) {
    substr(contigious.draft, mutation.indices[i], mutation.indices[i]) <- mutate(substr(contigious.draft, mutation.indices[i], mutation.indices[i]))
  }
  
  generate_scaffolds(contigious.draft, 
                     base.genome.scaffolds.length,
                     "DjMutatedScaffold: ",
                     paste0("Mutated_draf_v", version.number),
                     "./Data/Mutated_drafts/")
}

#Here we mutate one locus at a genome sequence, randomly with 0.25 chance to be mutated to A or T or C or G 
mutate <- function(base){
  Mutation <- ""
  switch(base,
         "A" = {    
           rand <- round(runif(1,1,4))
           switch (rand,
                   "1" =Mutation <- "T",
                   "2" =Mutation <- "C",
                   "3" =Mutation <- "G")
         },
         "T" = {    
           rand <- round(runif(1,1,4))
           switch (rand,
                   "1" =Mutation <- "A",
                   "2" =Mutation <- "C",
                   "3" =Mutation <- "G")
         },
         "C" = {    
           rand <- round(runif(1,1,4))
           switch (rand,
                   "1" =Mutation <- "T",
                   "2" =Mutation <- "A",
                   "3" =Mutation <- "G")
         },
         "G" = {    
           rand <- round(runif(1,1,4))
           switch (rand,
                   "1" =Mutation <- "T",
                   "2" =Mutation <- "C",
                   "3" =Mutation <- "A")
         }
  )
  return(Mutation)
}

#Here we get a draft genome and generate scaffolds as which the size of every one is declared by get_scaffolds_size
generate_scaffolds <- function(draft.genome, scaffolds.length, scaffold.header, draft.name, draft.path){
  scaffolds.number <- length(scaffolds.length)
  start.index <- 0
  print(paste0("number of scaffolds: ",scaffolds.number))
  for(i in 1:scaffolds.number){
    end.index <- start.index + scaffolds.length[i]
    print(paste0(scaffold.header,i," > ",start.index,":",end.index))
    scaff <- substr(draft.genome, start.index, end.index)
    index.str <- paste0(">",scaffold.header,i)
    if(i==1)
      write.table(paste(index.str,scaff,sep = '\n'), paste0(draft.path, draft.name,".txt"), append = FALSE,sep = '/n',col.names=FALSE,row.names=FALSE,quote = FALSE)
    else
      write.table(paste(index.str,scaff,sep = '\n'), paste0(draft.path, draft.name,".txt"), append = TRUE,sep = '/n',col.names=FALSE,row.names=FALSE, quote = FALSE)
    start.index <- (end.index+1)
  }
  write.table(scaffolds.length, paste0(draft.path, draft.name,"_scaff_lenght",".txt"), append = TRUE,sep = '/n',col.names=FALSE,row.names=FALSE, quote = FALSE)
} 

#Here we get an integer as the total size of a draft genome, an integer as the average size of scaffolds and generate 
# some numbers from a normal distribution with mean= scaff.mean.size and sd= variance factor
get_scaffolds_size <- function(total.size, scaff.mean.size, variance.factor){
  scaffolds.number <- floor(total.size/scaff.mean.size)
  scaffolds.length <- floor(rnorm(scaffolds.number, scaff.mean.size,sd = variance.factor*scaff.mean.size))  
  return(scaffolds.length)
}


mutate_draft(base.genome, 0.1, 1)




# start.indexV <- c()
# 
# onelinedata <- readLines("./Data/OneLineString.txt")
# # generate_scaffolds(onelinedata, 450, 0.1)
# 
# 
# get_scaffolds_size <- function(total.size, scaff.mean.size, variance.factor){
#   scaffolds.number <- floor(total.size/scaff.mean.size)
#   print(scaffolds.number)
#   scaffolds.length <- floor(rnorm(scaffolds.number, scaff.mean.size,sd = variance.factor*scaff.mean.size))  
#   return(scaffolds.length)
# }
# 
# 
# generate_scaffolds <- function(draft.genome, scaffolds.length, scaffold.header, draft.name, draft.path){
#   scaffolds.number <- length(scaffolds.length)
#   start.index <- 0
#   print(paste0("number of scaffolds: ",scaffolds.number))
#   for(i in 1:scaffolds.number){
#     end.index <- start.index + scaffolds.length[i]
#     print(paste0(scaffold.header,i," > ",start.index,":",end.index))
#     scaff <- substr(draft.genome, start.index, end.index)
#     index.str <- paste0(">",scaffold.header,i)
#     if(i==1)
#       write.table(paste(index.str,scaff,sep = '\n'), paste0(draft.path, draft.name,".txt"), append = FALSE,sep = '/n',col.names=FALSE,row.names=FALSE,quote = FALSE)
#     else
#       write.table(paste(index.str,scaff,sep = '\n'), paste0(draft.path, draft.name,".txt"), append = TRUE,sep = '/n',col.names=FALSE,row.names=FALSE, quote = FALSE)
#     start.index <- (end.index+1)
#   }
#   write.table(scaffolds.length, paste0(draft.path,draft.name,"_scaff_lenght",".txt"), append = TRUE,sep = '/n',col.names=FALSE,row.names=FALSE, quote = FALSE)
# } 



#generate_scaffolds(onelinedata, get_scaffolds_size(nchar(onelinedata), 450, 0.15), "DjsimulatedScaffolds:", "second_test", "./Data/" )


#get_scaffolds_size(nchar(onelinedata), 450, 0.15)



# generate_scaffolds <- function(draft.genome,scaff.mean.size,variance.factor){
#   total_size <- nchar(draft.genome)
#   scaffolds.number <- floor(total_size/scaff.mean.size)
#   scaffolds.length <- floor(rnorm(scaffolds.number,scaff.mean.size,sd = variance.factor*scaff.mean.size))  
#   start.index <- 0
#   print(paste0("number of scaffolds: ",scaffolds.number))
#   
#   for(i in 1:scaffolds.number){
#     end.index <- start.index + scaffolds.length[i]
#     print(paste0("> Dj:simulatedScaffold:",i," > ",start.index,":",end.index))
#     scaff <- substr(draft.genome, start.index, end.index)
#     index.str <- paste(">Dj:simulatedScaffold:",i)
#     if(i==1)
#       write.table(paste(index.str,scaff,sep = '\n'),'./simulated_scaffolds.txt',append = FALSE,sep = '/n',col.names=FALSE,row.names=FALSE,quote = FALSE)
#     else
#       write.table(paste(index.str,scaff,sep = '\n'),'./simulated_scaffolds.txt',append = TRUE,sep = '/n',col.names=FALSE,row.names=FALSE, quote = FALSE)
#     start.index <- (end.index+1)
#   }
#   write.table(scaffolds.length,"./simulated_scaffolds_length.txt",append = TRUE,sep = '/n',col.names=FALSE,row.names=FALSE, quote = FALSE)
# } 

getwd()

