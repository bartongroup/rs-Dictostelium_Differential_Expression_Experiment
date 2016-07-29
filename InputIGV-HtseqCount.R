	targets <- read.table(file="targets",header=TRUE)

        for(i in 1:nrow(targets)){
	files <- paste(targets[,1][i],"/","accepted_hits.bam",sep="")

	## sort by name and generate SAM files for htseq-count
	output <- paste(targets[,1][i],"_sort",sep="")
	command <- paste("samtools sort -n",files,output)
	system(command)
	Inp <- paste(targets[,1][i],"_sort.bam",sep="")
	out <- paste(targets[,1][i],"_sort.sam",sep="")
	command1 <- paste("samtools view -o",out,Inp)
	system(command1)
	
	## sort by position and index for IGV/IGB
	out <- paste(targets[,1][i],"-sortIGV",sep="")
	command2 <- paste("samtools sort",files,out)
	system(command2)
	Inp1 <- paste(targets[,1][i],"-sortIGV.bam",sep="")
	command3 <- paste("samtools index",Inp1)
	system(command3)
	}


	

	

