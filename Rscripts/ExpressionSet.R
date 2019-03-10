exampleSet <- ExpressionSet(assayData=exprs,phenoData=phenoData,
	experimentData=experimentData,annotation="hgu95av2")  
minimalSet <- ExpressionSet(assayData=exprs)
	1.assayData
	exprs <- as.matrix(read.table(exprsFile, header=TRUE, sep="\t",row.names=1,as.is=TRUE))
	2.phenoData <- new("AnnotatedDataFrame",data=pData, varMetadata=metadata)
				pData <- read.table(pDataFile,row.names=1, header=TRUE, sep="\t")
				all(rownames(pData)==colnames(exprs))
				metadata <- data.frame(labelDescription=
					c("Patient gender",
					"Case/control status",
					"Tumor progress on XYZ scale"),
					row.names=c("gender", "type", "score"))
	3.annotation <- "hgu95av2"
	4.experimentData <- new("MIAME",name="Pierre Fermat",lab="Francis Galton Lab",
			contact="pfermat@lab.not.exist",title="Smoking-Cancer Experiment",
			abstract="An example ExpressionSet",url="www.lab.not.exist",
			other=list(notes="Created from text files"))

ExpressionSet Basics
exampleSet$type
featureNames()
sampleNames()
varLabels()
exprs()
phenoData()
pData(phenoData)