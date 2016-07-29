	targets <- read.table(file="targets",header=TRUE)
	mainDir <-"/homes/rsingh/Ensemble"
	for(i in seq(along=targets[,1])){
	outDir <- paste("QC_",targets[,1][i],sep="")
        dir.create(file.path(mainDir,outDir))
        }

        for(i in 1:nrow(targets)){
	out <- paste("QC_",targets[,1][i],sep="")
        command <- paste("fastqc -o",out,"[-f fastq]",sep="")
	command1 <- paste(command,targets[,2][i],targets[,3][i])
	system(command1)
	}


