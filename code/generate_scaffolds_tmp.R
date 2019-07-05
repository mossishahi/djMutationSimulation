


## Here we Mutate the draft genome declared as the input of this function and
## generate a different draft genome which is Mutated! we mutate the locus with indices which are generated randomly
#finally we write the new draft genome in new txt file
mutate_draft <- function(base.genome, mutation.rate, version.number){
  base.genome.scaffolds <- as.character(base.genome$V1)
  base.genome.scaffolds.headers <- base.genome.scaffolds[grep(">",base.genome.scaffolds)]
  base.genome.scaffolds.sequences <- base.genome.scaffolds[-grep(">",base.genome.scaffolds)]
  base.genome.scaffolds.length <- nchar(base.genome.scaffolds.sequences)
  #print("--------")
  print(base.genome.scaffolds.length)
  contigious.draft <- paste0(base.genome.scaffolds.sequences,collapse = '')
  mutation.count <- round(mutation.rate * nchar(contigious.draft))
  print(paste0("Number of mutations:",mutation.count))
  mutation.indices <- sample(1:nchar(contigious.draft), mutation.count)
  # print("These indices are chosen randomly to be Mutated")
  # print(mutation.indices)

  for (i in 1:mutation.count) {
    substr(contigious.draft, mutation.indices[i], mutation.indices[i]) <- mutate(substr(contigious.draft, mutation.indices[i], mutation.indices[i]))
  }
  # call this function to write mutated genome in file
  generate_scaffolds(contigious.draft, 
                     base.genome.scaffolds.length,
                     paste0("DjMutatedGenome_V:",version.number," scaffoldNumber :"),
                     paste0("Mutated_draf_v", version.number),
                     "./Data/Mutated_drafts/")
}

#Here we mutate one locus at a genome sequence, randomly with 0.25 chance to be mutated to A or T or C or G, by the chance of 0.25
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

## Here we defined a function named generate_scaffolds, which gets a fragmented draft genome and 
# and generates scaffolds wih sizes you specify in the vector, scaffolds.length
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
# some numbers from a normal distribution with mean=scaff.mean.size and sd= variance factor
get_scaffolds_size <- function(total.size, scaff.mean.size, variance.factor){
  scaffolds.number <- floor(total.size/scaff.mean.size)
  scaffolds.length <- floor(rnorm(scaffolds.number, scaff.mean.size,sd = variance.factor*scaff.mean.size))  
  return(scaffolds.length)
}


mutate_draft(base.genome, 0.1, 1)


