
#### load matter

library(matter)
# nchunks <- 4L * max(1L, detectCores() - 1L)
matter_defaults(verbose = FALSE, nchunks = 16L)
raw <- read.csv("Selevsek2015_DIA_Spectronaut_input-001.xls", sep="\t")
annotation <- read.csv("Selevsek2015_DIA_Spectronaut_annotation.csv")
output = SpectronauttoMSstatsFormat(raw, annotation = annotation)

 
method = "TMP"
cens = "NA"
impute = TRUE
MSstatsConvert::MSstatsLogsSettings(FALSE)
input = MSstatsPrepareForDataProcess(output, 2, NULL)
input = MSstatsNormalize(input, "EQUALIZEMEDIANS")
input = MSstatsMergeFractions(input)
input = MSstatsHandleMissing(input, "TMP", TRUE, "NA", 0.999)
input = MSstatsSelectFeatures(input, "all")
processed = getProcessed(input)
input = MSstatsPrepareForSummarization(input, method, impute, cens, FALSE)
protein_indices = split(seq_len(nrow(input)), list(input$PROTEIN))
num_proteins = length(protein_indices)
summarized_results = vector("list", num_proteins)
single_protein = input[protein_indices[[1]],]

bro = function(i) {
    single_protein = input[protein_indices[[i]],]
    MSstatsSummarizeSingleTMP(
        single_protein, impute, censored_symbol, remove50missing)
}

# input object not found... This is where we would add the input to the cluster, but 
# we should probably just use it as the input.  
summarized_results = chunkLapply(seq_len(num_proteins), bro, BPPARAM=SnowParam(workers=4))


summarized_results = chunkLapply(seq_len(num_proteins), function(i) {
    single_protein = input[protein_indices[[i]],]
    MSstatsSummarizeSingleTMP(
        single_protein, impute, censored_symbol, remove50missing)
}, BPPARAM=SnowfastParam(workers=4))
# Turn input into matter list.... Refactor to a matrix in each list.
summarized = MSstatsSummarizeWithSingleCore(input, method, impute, cens, FALSE, TRUE)
length(summarized) # list of summarization outputs for each protein
head(summarized[[1]][[1]])

protein_ids = matter_vec(seq_len(protein))
Y = chunk_lapply(protein_ids, FUN = function(i) { i + 5
}, BPPARAM=SnowParam(workers=4))
# as.list.default(X) no method for coercing this S4 class to a vector - BiocParallel
# Chunks making it not work for some reason....
# definitely default setting wrong here....  Github version don't work.
fast
Y = lapply(seq_len(num_proteins), FUN = function(i) { i + 5
})


matter_defaults(verbose=TRUE, nchunks=20L)

## time cluster startup

# 4 workers
sp4 <- SnowParam(workers=4)
sfp4 <- SnowfastParam(workers=4)

system.time(bpstart(sp4))
  #  user  system elapsed 
  # 0.052   0.019   0.472 
system.time(bpstart(sfp4))
  #  user  system elapsed 
  # 0.003   0.002   0.089 

bpstop(sp4)
bpstop(sfp4)

# 8 workers
sp8 <- SnowParam(workers=8)
sfp8 <- SnowfastParam(workers=8)

system.time(bpstart(sp8))
  #  user  system elapsed 
  # 0.053   0.051   0.933 
system.time(bpstart(sfp8))
  #  user  system elapsed 
  # 0.002   0.003   0.091 

bpstop(sp8)
bpstop(sfp8)

# 12 workers
sp12 <- SnowParam(workers=12)
sfp12 <- SnowfastParam(workers=12)

system.time(bpstart(sp12))
  #  user  system elapsed 
  # 0.057   0.072   1.320 
system.time(bpstart(sfp12))
  #  user  system elapsed 
  # 0.002   0.004   0.101 

bpstop(sp12)
bpstop(sfp12)

# 16 workers
sp16 <- SnowParam(workers=16)
sfp16 <- SnowfastParam(workers=16)

system.time(bpstart(sp16))
  #  user  system elapsed 
  # 0.059   0.089   1.738 
system.time(bpstart(sfp16))
  #  user  system elapsed 
  # 0.004   0.005   0.123 

bpstop(sp16)
bpstop(sfp16)

#### load data

dir <- "~/workspace/MSstats/"
file <- "50um_sample1_2_3_pos.imzML"
path <- file.path(dir, file)

p <- CardinalIO::parseImzML(path, ibd=TRUE, check=FALSE)
x <- p$ibd$intensity
print(x)
# <65133 length> matter_list :: out-of-memory list
# (5.49 MB real | 0 bytes shared | 9.67 GB virtual)

#### compare SnowParam versus SnowfastParam

matter_defaults(serialize=TRUE)

# using file-based matter 'x'
# reading data on manager

system.time(chunkLapply(x, summary,
	BPPARAM=SnowParam(workers=4)))
  #  user  system elapsed 
  # 15.003   5.921  37.243 

system.time(chunkLapply(x, summary,
	BPPARAM=SnowfastParam(workers=4)))
  #  user  system elapsed 
  # 3.711   5.382  22.833

matter_defaults(serialize=FALSE)

# using file-based matter 'x'
# reading data on workers

system.time(chunkLapply(x, summary,
	BPPARAM=SnowParam(workers=12)))
  #  user  system elapsed 
  # 0.427   0.454  16.967 

system.time(chunkLapply(x, summary,
	BPPARAM=SnowfastParam(workers=12)))
  #  user  system elapsed 
  # 0.239   0.025  15.381 

# using file-based matter 'x'
# comparison to sequential processing

system.time(chunkLapply(x, summary,
	BPPARAM=NULL))
 #   user  system elapsed 
 # 53.765   6.327  60.265 

#### compare file-based versus R matrix

matter_defaults(serialize=NA)

x2 <- as.matrix(x)
  #  user  system elapsed 
  # 1.769   2.076   3.851 
mem(x2)
#     real   shared  virtual 
# 19.34 GB     0 GB     0 GB 

# using R matrix 'x2'
# serializing data in memory

system.time(chunkLapply(x2, summary,
	BPPARAM=SnowParam(workers=12)))
 #   user  system elapsed 
 # 11.416   6.411  35.800 

system.time(chunkLapply(x2, summary,
	BPPARAM=SnowfastParam(workers=12)))
  #  user  system elapsed 
  # 1.247   5.928  20.188 

system.time(chunkLapply(x2, summary,
	BPPARAM=MulticoreParam(workers=12)))
 #   user  system elapsed 
 # 47.229  12.600  10.791

# mcp12 <- MulticoreParam(workers=12)
# memtime(chunkLapply(x2, summary, BPPARAM=mcp12),
# 	BPPARAM=mcp12)
# crashes....

rm(x2)
mem(reset=TRUE)

#### compare shared memory versus alternatives

x3 <- fetch(x, BPPARAM=SnowfastParam(workers=12))
print(x3)
# <65133 length> matter_list :: out-of-memory list
# (5.49 MB real | 19.34 GB shared | 19.34 GB virtual)

# using file-based matter 'x'
# using shared memory-based matter 'x3'

system.time(chunkLapply(x, summary,
	BPPARAM=SnowfastParam(workers=12)))
  #  user  system elapsed 
  # 0.235   0.022  15.351

system.time(chunkLapply(x3, summary,
	BPPARAM=SnowfastParam(workers=12)))
  #  user  system elapsed 
  # 0.242   0.023  15.424 

sfp12 <- SnowfastParam(workers=12)
memtime(chunkLapply(x3, summary, BPPARAM=sfp12),
	BPPARAM=sfp12)
# $cluster$end
#         node      real  shared max real max shared    temp
# 1  localhost 158.19 MB 0 bytes  0.16 GB    0 bytes 0 bytes
# 2  localhost 158.19 MB 0 bytes  0.16 GB    0 bytes 0 bytes
# 3  localhost 190.28 MB 0 bytes  1.04 GB    0 bytes 0 bytes
# 4  localhost 190.28 MB 0 bytes  1.04 GB    0 bytes 0 bytes
# 5  localhost 190.28 MB 0 bytes  1.04 GB    0 bytes 0 bytes
# 6  localhost 190.28 MB 0 bytes  1.04 GB    0 bytes 0 bytes
# 7  localhost 190.28 MB 0 bytes  1.04 GB    0 bytes 0 bytes
# 8  localhost 190.28 MB 0 bytes  1.04 GB    0 bytes 0 bytes
# 9  localhost 190.28 MB 0 bytes  1.04 GB    0 bytes 0 bytes
# 10 localhost 190.28 MB 0 bytes  1.04 GB    0 bytes 0 bytes
# 11 localhost 190.28 MB 0 bytes  1.04 GB    0 bytes 0 bytes
# 12 localhost 190.28 MB 0 bytes  1.04 GB    0 bytes 0 bytes

# $cluster$total
# [1] 10.7 GB

# $overhead
#     real   shared 
# 10.83 GB     0 GB 

# $total
# [1] 30.39 GB

# $time
#    user  system elapsed 
#   0.239   0.018   9.605 
