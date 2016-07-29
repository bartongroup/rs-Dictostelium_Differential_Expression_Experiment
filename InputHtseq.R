	files <- list.files(pattern="*-sortIGV.bam_duplRem.bam$")
	for(i in 1:length(files)){
	output <- paste(files[i],"_sort",sep="")
	command <- paste("samtools sort -n",files[i],output)
	system(command)
	Inp <- paste(files[i],"_sort.bam",sep="")
	out <- paste(files[i],"_sort.sam",sep="")
        command1 <- paste("samtools view -o",out,Inp)
	system(command1)
	}

