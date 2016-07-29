	a <- read.table(file="EnrichedGeneCategory.txt",sep="\t",header=TRUE)
	cat <- read.table(file="Ontology/OntologyUsed",header=TRUE)

	common <- cat[which(cat$GO %in% a$category),]
	common <- common[,-3]
	df1 <- aggregate(common[2],common[-2], FUN = function(X) paste(unique(X), collapse=", "))
	DEG <- read.table(file="DEGPvalue.txt",sep="\t",header=TRUE)
	b <- DEG[which(rownames(DEG) %in% df1$ID),]
	write.table(b,file="DEGListFC-EnrichedCatWF.txt",sep="\t")
	c <- common[which(common$ID %in% rownames(DEG)), ] #### DEG and GO categories
	write.table(c,file="DEGList-EnrichedCatWF.txt",sep="\t")	

	con <- file("AllEnrichedGoCategory.txt","r")
	term <- readLines(con)
#	goid <- data.frame()
#	 for(i in 1:length(term)){
#	 id <- term[grep(tomatch[1],term)]
#	 def <- term[grep(tomatch[2],term)]
#	 d <- cbind(id,def)
#	 goid <- rbind(goid,d)
#	 }
	tomatch <- c("GOID","Term")
	id <- term[grep(tomatch[1],term)]
	def <- term[grep(tomatch[2],term)]
	goid <- data.frame(id,def)
	goid1 <- as.data.frame(sapply(goid,gsub,pattern="GOID:",replacement=""))
	goid1 <- as.data.frame(sapply(goid1,gsub,pattern="Term:",replacement=""))
	colnames(goid1) <- c("GO","Term")
	write.table(goid1,file="DEGListTerms-EnrichedCatWF.txt",sep="\t")

	
	


	#df2 <- 
	#combine <- cbind(b,common[which(rownames(b) %in% df1$ID),])
	
