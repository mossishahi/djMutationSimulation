olinedata <- readLines("./Data/OneLineString.txt")

# data <- readLines("./Data/draft_genome.txt")
# data_length <- nchar(data)
# barplot(data_length)
# plot(data_length)
# concat <- paste(data,collapse = '')
# concat
# write.table(concat, ,file="./Data/OneLineString.txt",sep = '\n',row.names = FALSE,col.names = FALSE)
# write.table(data,file="./Data/commaseperated.txt",sep = ',,,',row.names = FALSE,col.names = FALSE)




start.indexV <- c()

generate_scaffolds <- function(draft.genome,scaff.mean.size,variance.factor){
  total_size <- nchar(draft.genome)
  scaffolds.number <- floor(total_size/scaff.mean.size)
  scaffolds.length <- floor(rnorm(scaffolds.number,scaff.mean.size,sd = variance.factor*scaff.mean.size))  
  start.index <- 0
  print(paste0("number of scaffolds: ",scaffolds.number))
 
   for(i in 1:scaffolds.number){
    end.index <- start.index + scaffolds.length[i]
    print(paste0("> Dj:simulatedScaffold:",i," > ",start.index,":",end.index))
    scaff <- substr(draft.genome, start.index, end.index)
    index.str <- paste(">Dj:simulatedScaffold:",i)
    if(i==1)
     write.table(paste(index.str,scaff,sep = '\n'),'./simulated_scaffolds.txt',append = FALSE,sep = '/n',col.names=FALSE,row.names=FALSE,quote = FALSE)
    else
      write.table(paste(index.str,scaff,sep = '\n'),'./simulated_scaffolds.txt',append = TRUE,sep = '/n',col.names=FALSE,row.names=FALSE, quote = FALSE)
    start.index <- (end.index+1)
  }
  
  write.table(scaffolds.length,"./simulated_scaffolds_length.txt",append = TRUE,sep = '/n',col.names=FALSE,row.names=FALSE, quote = FALSE)
} 

getwd()       

slength <- readLines("./Data/simulated_scaffolds_length.txt")
slength <- as.numeric(slength)

#slength <- sort(slegnth)
plot(slength)

ggplot(data = slength)+geom_density2d(data = slength)

geom_density_2d(data=slength)

draft.scaff.length <-readLines("./Data/scaffolds_length.txt")
draft.scaff.length <- as.numeric(draft.scaff.length)
plot(density(draft.scaff.length))
length(draft.scaff.length)
hist(draft.scaff.length)

draft.scaff.length <- sort(draft.scaff.length,decreasing = TRUE)
draft.scaff.length[1]
draft.scaff.length <- draft.scaff.length[1:80000]
plot(density(draft.scaff.length))





## Here we Mutate the draft genome declared as the input of this function and
## generate n different draft genoms which are Mutated against input draft
#finally we write each draft genome in different txt file

onelinedata <- readLines("./Data/OneLineString.txt")
length(onelinedata)

typeof(onelinedata)

base.genome <- read.delim("./Data/simulated_scaffolds.txt", sep = "\n", header = FALSE)


mutate_draft <- function(base.genome, .mutation.rate, version.number){
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
  print("These indice are chosen randomly to be Mutated")
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
generate_scaffolds <- function(draft.genome, scaffolds.length, scaffold.header, draft.name, draft.path){
  scaffolds.number <- length(scaffolds.length)
  start.index <- 0
  print(paste0("number of scaffolds: ",scaffolds.number))
  for(i in 1:scaffolds.number){
    end.index <- start.index + scaffolds.length[i]
    print(paste0(scaffold.header,i," > ",start.index,":",end.index))
    scaff <- substr(draft.genome, start.index, end.index)
    index.str <- paste0(">",scaffold.header,i)
    print("***")
    print(nchar(paste(index.str,scaff,sep = '\n')))
    if(i==1)
      write.table(paste(index.str,scaff,sep = '\n'), paste0(draft.path, draft.name,".txt"), append = FALSE,sep = '/n',col.names=FALSE,row.names=FALSE,quote = FALSE)
    else
      write.table(paste(index.str,scaff,sep = '\n'), paste0(draft.path, draft.name,".txt"), append = TRUE,sep = '/n',col.names=FALSE,row.names=FALSE, quote = FALSE)
    start.index <- (end.index+1)
  }
  write.table(scaffolds.length, paste0(draft.path, draft.name,"_scaff_lenght",".txt"), append = TRUE,sep = '/n',col.names=FALSE,row.names=FALSE, quote = FALSE)
} 
get_scaffolds_size <- function(total.size, scaff.mean.size, variance.factor){
  scaffolds.number <- floor(total.size/scaff.mean.size)
  scaffolds.length <- floor(rnorm(scaffolds.number, scaff.mean.size,sd = variance.factor*scaff.mean.size))  
  return(scaffolds.length)
}



# generate_scaffolds <- function(draft.genome, scaffolds.length, scaffold.header, draf.name, draft.path){
#   scaffolds.number <- length(scaffolds.length)
#   start.index <- 0
#   print(paste0("number of scaffolds: ",scaffolds.number))
#   for(i in 1:scaffolds.number){
#     end.index <- start.index + scaffolds.length[i]
#     print(paste0(scaffold.header,i," > ",start.index,":",end.index))
#     scaff <- substr(draft.genome, start.index, end.index)
#     index.str <- paste(scaffold.header,i)
#     if(i==1)
#       write.table(paste(index.str,scaff,sep = '\n'), paste0(draft.path, draft.name,".txt"), append = FALSE,sep = '/n',col.names=FALSE,row.names=FALSE,quote = FALSE)
#     else
#       write.table(paste(index.str,scaff,sep = '\n'), paste0(draf.draft.path, draft.name,".txt"), append = TRUE,sep = '/n',col.names=FALSE,row.names=FALSE, quote = FALSE)
#     start.index <- (end.index+1)
#   }
#   write.table(scaffolds.length, paste0(draft.path,"_scaff_lenght",".txt"), append = TRUE,sep = '/n',col.names=FALSE,row.names=FALSE, quote = FALSE)
# } 




