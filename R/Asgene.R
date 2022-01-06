#' Title
#' AsgeneDB: A functional gene database for metagenomic profiling of arsenic metabolism
#' Description
#' A manually curated arsenic functional gene database/tool (AsgeneDB) was developed for rapid and accurate metagenomic analysis.
#' @param analysis Choose the target you want to analyze,abundance or taxonomy; type = "character"; default = "abundance"
#' @param workdir Specify directory for sequence file location;type = "character"; default = "./"
#' @param method Specify the database searching tool you plan to use, currently diamond, usearch and blast are supported; type = "character"; default = "diamond"
#' @param toolpath Specify directory for searching tool location; type = "character"; default = "./"
#' @param seqtype Specify the extensions of your sequence files, e.g. fastq, fastq.gz, fasta,fasta.gz, fq, fq.gz, fa, fa.gz; type = "character"; default = "fasta"
#' @param filetype Specify your sequence type, nucl or prot; type = "character"; default = "nucl"
#' @param out Specify the directory for the final result file; type = "character"; default = "./"
#'
#' @return
#' @export
#'
#' @examples
#' library(Asgene)
#' Asgene("abundance","./","diamond","./","nucl","fasta","./")
#' Asgene("taxonomy","./","diamond","./","nucl","fasta","./")
Asgene <- function(analysis="abundance",workdir="./",method="diamond",toolpath="./",search_parameters = "-e 1e-10 -p 28 --query-cover 80 --id 50",seqtype="fasta",filetype="nucl",out="./"){
  #Install dependent packages
  if (!requireNamespace("dplyr", quietly = TRUE))
    install.packages("dplyr")
  library("dplyr")
  #user-define metagenomic comparison parameters
  search_parameters <- "-e 1e-4 -p 28 --query-cover 80 --id 50 "
  system("mkdir sample_file")

  #Call comparison tool
  if ( method == "diamond" ) {
    system(paste(toolpath,"diamond makedb --in ./AsgeneDB.fa --db ./AsgeneDB",sep = ""));
    file <- list.files(path=workdir, pattern=filetype);
    for (i in file) {
      file_1 <- paste(workdir,i,sep="");
      out <- gsub(filetype,"diamond",file);
      out <- paste("./sample_file/",out,sep = "");
      if (seqtype == "nucl") {
        system(paste(toolpath,"diamond blastx ",search_parameters," -d ./AsgeneDB.dmnd -q ",file_1," -o ",out,sep = ""))}
      if (seqtype == "prot"){
        system(paste(toolpath,"diamond blastp ",search_parameters," -d ./AsgeneDB.dmnd -q ",file_1," -o ",out,sep = ""))}
    }
  }

  if ( method == "usearch") {
    if ( grepl("gz",filetype) == "TRUE" ) {
      print("Only fastq and fasta files are supported by usearch!\n");
      q();
    }
    file <- list.files(path=workdir, pattern= filetype);
    for (i in file) {
      file_1 <- paste(workdir,i,sep="");
      out <- gsub(filetype,"usearch",file);
      out <- paste("./sample_file/",out,sep = "");
      system(paste(toolpath,"usearch -usearch_global ",file_1," -db ./AsgeneDB.fa ",search_parameters," -blast6out ",out,sep = ""))}
  }

  if ( method == "blast" ) {
    if ( grepl("gz|fastq|fq",filetype,perl = TRUE) == "TRUE" ) {
      print("Only fasta files are supported by blast program!");
      q();
    }
    file <- list.files(path=workdir, pattern= filetype);
    system(toolpath,"makeblastdb -dbtype prot -input_type fasta -in ./total_100% -out ./total_100%");
    for (i in file) {
      file_1 <- paste(workdir,i,sep="");
      out <- gsub(filetype,"blast",file);
      out <- paste("./sample_file/",out,sep = "");
      if (seqtype == "nucl") {
        system(paste(toolpath,"blastx -db ./AsgeneDB.fa ",search_parameters," -query ",file_1," -out ",out,sep = ""))}
      if (seqtype == "prot"){
        system(paste(toolpath,"blastp -db ./AsgeneDB.fa ",search_parameters," -query ",file_1," -oout ",out,sep = ""))}}
  }

  #Merge metagenomic PE files (take union)
  list <-read.table("./sampleinfo.txt",sep = "\t",header =F)
  system("mkdir sample_merge")
  for (i in list[,1]){
    a <- read.table(file=paste("./sample_file/",i,"_R1.",method,sep=""),sep = "\t")
    a1 <- data.frame(a)
    b <- read.table(file=paste("./sample_file/",i,"_R2.",method,sep=""),sep = "\t")
    b1 <- data.frame(b)
    e <- rbind(a1,b1)
    e1 <- aggregate(e[,3]~e[,1],data=e,FUN="max")
    colnames(e1)[1] <- "V1"
    colnames(e1)[2] <- "V3"
    e1 <- merge(e1,e,by=c("V1","V3"))
    e1 <- e1 %>% filter(!duplicated(V1))
    write.table(e1,file=paste("sample_merge/",i,".",method,sep=""), sep= " ",quote = FALSE,row.names = FALSE,col.names = F)
  }

  #count the abundance values of each sample and standardize
  if ( analysis == "abundance" ) {
    file.name <- list.files(path="./sample_merge/",pattern = method)
    for (i in file.name){
      a <- read.table(file=paste("./sample_merge/",i,sep = ""),sep = " ",header =F)
      i<- gsub(paste(".",method,sep = ""),"",i)
      a1<- mutate(a,v=i)
      #add total reads
      list<-data.frame(list)
      colnames(list)[1]<-"v"
      sample <- list$v
      t<- list [grep(i,list$v),]
      a2<- merge(t,a1,by="v")
      a2<- a2[,c(1,2,3,5)]
      b <- read.table("./asgene.map",sep = "")
      colnames(a2)[4] <- "pi"
      colnames(b)[1] <- "pi"
      colnames(b)[2] <- "gene"
      g<-merge(a2,b,by="pi")
      g<-g[,-6]
      colnames(g)[3]<-"totalreads"
      u <- read.table("./length.txt",sep = "\t",header = F)
      colnames(u)[1] <- "pi"
      result <- merge(g,u,by="pi")
      colnames(result)[6]<-"length"
      x <- result
      x <-x[order(x$gene),]

      y<-count(x,x$pi)
      colnames(y)[1]<-"pi"
      v<-merge(x,y,by="pi")
      v1 <- v[!duplicated(v$pi),]

      #Formula=SUM protein{[mapped reads/(total reads*protein length)]*10^9}
      v1$totalreads <- as.numeric(v1$totalreads)
      df<- data.frame(v1)
      df<- mutate(df,c=v1$n*10^9/(v1$length*3))
      df<- mutate(df,d=c/v1$totalreads)
      #group_by_sum
      df1<-aggregate(df$d,by=list(df$gene),sum)
      colnames(df1) <- c("gene",as.character(i))
      gsub(".diamond","",df1)
      write.table(df1, file = paste("./sample_merge/",i, ".csv",sep = ""),sep=",",quote = FALSE,row.names = FALSE)
    }

    a <- list.files(path="./sample_merge/",pattern=".csv")
    n <- length(a)
    merge.data <- read.csv(file = paste("./sample_merge/",a[1],sep = ""),header=T,sep=",")
    merge.data1 <- data.frame(merge.data)
    for (i in 2:n){
      new.data =read.csv(file=paste("./sample_merge/",a[i],sep = ""),sep = ",",header =T)
      merge.data1 = merge(merge.data1,new.data,by="gene",all=T)
    }
    write.csv(merge.data1,file=paste(out,"sample_abundance.csv",sep=""),row.names=F)
  }

  #count functional gene drive species of all samples
  if ( analysis == "taxonomy" ) {
    list <- list.files(path="./sample_merge/",pattern = method)
    system("mkdir gene_tax")
    for (i in list){
      a <- read.table(file=paste("./sample_merge/",i,sep = ""),sep = " ",header =F)
      colnames(a)[1] = "reads_id"
      colnames(a)[3] = "protein_id"
      a <- a[,c(1,3)] %>% filter(!duplicated(a[,1]))
      i<- gsub(".diamond","",i)
      a1<- mutate(a,"sample"=i)
      id_gene_tax_pathway <- read.csv("./id_gene_tax_pathway_total.csv",sep=",",header = T)
      sample_gene_tax_pathway <- merge(a1,id_gene_tax_pathway,by="protein_id" )
      write.table(sample_gene_tax_pathway, file = paste("./gene_tax/",i,".csv",sep = ""),sep=",",quote = FALSE,row.names = FALSE)
    }

    a <- list.files(path="./gene_tax/",pattern = ".csv")
    n <- length(a)
    merge.data <- read.csv(file = paste("./gene_tax/",a[1],sep = ""),header=T,ssep=",")
    merge.data1 <- data.frame(merge.data)
    for (i in 2:n){
      new.data =read.csv(file=paste("./gene_tax/",a[i],sep = ""),sep = ",",header =T)
      merge.data1 = merge(merge.data1,new.data,all=T)
    }
    write.csv(merge.data1,file=paste(out,"sample_gene_tax_pathway.csv",sep=""),row.names=F)
  }
}

