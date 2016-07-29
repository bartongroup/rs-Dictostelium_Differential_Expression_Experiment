	targets <- read.table(file="targets",header=TRUE)
	gtf <- "DictyEnsemble.gtf"	
	for(i in 1:nrow(targets)){
	Inp <- paste(targets[,1][i],"_sort.sam",sep="")
	count <- paste(targets[,1][i],".count",sep="")
	command <- paste("htseq-count -s no -a 10 ",Inp, gtf," > ",count)
	system(command)
	}
		
