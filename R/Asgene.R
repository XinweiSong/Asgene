#' AsgeneDB: A functional gene database for metagenomic profiling of arsenic metabolism#' AsgeneDB: A functional gene database for metagenomic profiling of arsenic metabolism
#' @description A manually curated arsenic functional gene database (AsgeneDB) and R package (Asgene package) are developed for rapid and accurate metagenomic analysis.
#' @param analysis Choose the target you want to analyze,abundance or taxonomy; type = "character"; default = "abundance"
#' @param workdir Specify directory for sequence file location;type = "character"; default = "./"
#' @param method Specify the database searching tool you plan to use, currently diamond, usearch and blast are supported; type = "character"; default = "diamond"
#' @param toolpath Specify directory for searching tool location; type = "character"; default = "./"
#' @param search_parameters Define metagenomic comparison parameters; type = "character"; default = "-e 1e-4 -p 28 --query-cover 80 --id 50"
#' @param seqtype Specify your sequence type, nucl or prot; type = "character"; default = "nucl"
#' @param filetype Specify the extensions of your sequence files, e.g., fastq, fastq.gz, fasta,fasta.gz, fq, fq.gz, fa, fa.gz; type = "character"; default = "fasta"
#' @param PE Specifies whether your metagenomic data were PE files, noting that PE files are named "_R1 with extension" or "_R2 with extension" (e.g., XINWEI_R1.fasta and XINWEI_R2.fasta); type = "logical"; default = TRUE
#' @param output Specify the directory for the final result file; type = "character"; default = "./"
#'
#' @return
#' @export
#'
#' @examples
#' library(Asgene)
#' Asgene(anlysis = "abundance", workdir = "./", method = "diamond", toolpath = "./", search_parameters = "-e 1e-4 -p 28 --query-cover 80 --id 50", seqtype = "nucl", filetype = "fasta", PE = TRUE, out = "./")
Asgene <- function(analysis = "abundance", workdir = "./", method = "diamond", toolpath = "./", search_parameters = "-e 1e-4 -p 28 --query-cover 80 --id 50", seqtype = "nucl", filetype = "fasta", PE = TRUE, output = "./") {
  setwd("./")
  # Install dependent packages
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr")
  }
  if (!requireNamespace("seqinr", quietly = TRUE)) {
    install.packages("seqinr")
  }
  library(dplyr)
  library(seqinr)
  # Call comparison tool
  if (method == "diamond") {
    system("mkdir sample_file")
    write.fasta(sequences = AsgeneDB, names =names(AsgeneDB), file.out = 'AsgeneDB')
    system(paste(toolpath, "diamond makedb --in AsgeneDB --db ./AsgeneDB", sep = ""))
    file <- list.files(path = workdir, pattern = filetype)
    for (i in file) {
      file_1 <- paste(workdir, i, sep = "")
      out <- gsub(filetype, "diamond", i)
      out <- paste("./sample_file/", out, sep = "")
      if (seqtype == "nucl") {
        system(paste(toolpath, "diamond blastx ", search_parameters, " -d ./AsgeneDB.dmnd -q ", file_1, " -o ", out, sep = ""))
      }
      if (seqtype == "prot") {
        system(paste(toolpath, "diamond blastp ", search_parameters, " -d ./AsgeneDB.dmnd -q ", file_1, " -o ", out, sep = ""))
      }
    }
  }
  else if (method == "usearch") {
    if (grepl("gz", filetype)) {
      stop("Only fastq and fasta files are supported by usearch!")
    } else {
      system("mkdir sample_file")
      file <- list.files(path = workdir, pattern = filetype)
      for (i in file) {
        file_1 <- paste(workdir, i, sep = "")
        out <- gsub(filetype, "usearch", i)
        out <- paste("./sample_file/", out, sep = "")
        system(paste(toolpath, "usearch -usearch_global ", file_1, " -db AsgeneDB ", search_parameters, " -blast6out ", out, sep = ""))
        print(i)
      }
    }
  } else if (method == "blast") {
    if (grepl("gz|fastq|fq", filetype, perl = TRUE)) {
      stop("Only fasta files are supported by blast program!")
    } else {
      system("mkdir sample_file")
      file <- list.files(path = workdir, pattern = filetype)
      system(toolpath, "makeblastdb -dbtype prot -input_type fasta -in AsgeneDB -out ./AsgeneDB")
      for (i in file) {
        file_1 <- paste(workdir, i, sep = "")
        out <- gsub(filetype, "blast", i)
        out <- paste("./sample_file/", out, sep = "")
        if (seqtype == "nucl") {
          system(paste(toolpath, "blastx -db ./AsgeneDB ", search_parameters, " -query ", file_1, " -out ", out, sep = ""))
        }
        if (seqtype == "prot") {
          system(paste(toolpath, "blastp -db ./AsgeneDB ", search_parameters, " -query ", file_1, " -out ", out, sep = ""))
        }
        print(i)
      }
    }
  } else {
    stop("Specify the database searching tool you plan to use!")
  }

  # Merge metagenomic PE files (take union)
  list <- read.table("./sampleinfo.txt", sep = "\t", header = F)
  if (PE == TRUE) {
    system("mkdir sample_merge")
    for (i in list[, 1]) {
      info1 <- file.info(paste("./sample_file/", i, "_R1.", method, sep = ""))
      info2 <- file.info(file = paste("./sample_file/", i, "_R2.", method, sep = ""))
      if (info1$size != 0 && info2$size != 0) {
        file_R1 <- read.table(file = paste("./sample_file/", i, "_R1.", method, sep = ""), sep = "\t")
        file_R1 <- data.frame(file_R1)
        file_R2 <- read.table(file = paste("./sample_file/", i, "_R2.", method, sep = ""), sep = "\t")
        file_R2 <- data.frame(file_R2)
        file_total <- rbind(file_R1, file_R2)
        file_total_2 <- aggregate(file_total[, 3] ~ file_total[, 1], data = file_total, FUN = "max")
        file_total_2 <- as.data.frame(file_total_2)
        names(file_total_2)[1] <- "V1"
        names(file_total_2)[2] <- "V3"
        file_total_2 <- merge(file_total_2, file_total_1, by = c("V1", "V3"))
        file_total_2 <- file_total_2 %>% filter(!duplicated(V1))
        write.table(file_total_2, file = paste("sample_merge/", i, ".", method, sep = ""), sep = " ", quote = FALSE, row.names = FALSE, col.names = F)
        print(i)
      } else {
        next
      }
    }
  }

  # count the abundance values of each sample and standardize
  if (analysis == "abundance") {
    if (PE == TRUE) {
      file.name <- list.files(path = "./sample_merge/", pattern = method)
    } else {
      file.name <- list.files(path = "./sample_file/", pattern = method)
    }
    for (i in file.name) {
      if (PE == TRUE) {
        file <- read.table(file = paste("./sample_merge/", i, sep = ""), sep = " ", header = F)
      } else {
        info <- file.info(paste("./sample_file/", i, sep = ""))
        if (info$size != 0) {
          file <- read.table(file = paste("./sample_file/", i, sep = ""), sep = "\t", header = F)
        } else {
          next
        }
      }
      i <- gsub(paste(".", method, sep = ""), "", i)
      file_1 <- mutate(file, v = i)
      # add total reads
      list <- data.frame(list)
      names(list)[1] <- "v"
      t <- list[grep(i, list$v), ]
      file_2 <- merge(t, file_1, by = "v")
      if (PE == TRUE) {
        file_2 <- file_2[, c(1, 2, 3, 5)]
      } else {
        file_2 <- file_2[, c(1, 2, 3, 4)]
      }
      asgene_map <- asgene.map
      asgene_map <- as.data.frame(asgene_map)
      file_2 <- as.data.frame(file_2)
      names(file_2)[4] <- "pi"
      names(asgene_map)[1] <- "pi"
      names(asgene_map)[2] <- "gene"
      result_g <- merge(file_2, asgene_map, by = "pi")
      result_g <- result_g[, -6]
      names(result_g)[3] <- "totalreads"
      length <- length
      names(length)[1] <- "pi"
      result <- merge(result_g, length, by = "pi")
      result <- as.data.frame(result)
      names(result)[6] <- "length"
      result_x <- result
      result_x <- result_x[order(result_x$gene), ]

      result_y <- dplyr::count(result_x, result_x$pi)
      result_y <- as.data.frame(result_y)
      names(result_y)[1] <- "pi"
      result_v <- merge(result_x, result_y, by = "pi")
      result_v1 <- result_v[!duplicated(result_v$pi), ]

      # Formula=SUM protein{[mapped reads/(total reads*protein length*3)]*10^9}
      result_v1$totalreads <- as.numeric(result_v1$totalreads)
      result_df <- data.frame(result_v1)
      result_df <- mutate(result_df, c = result_v1$n * 10^9 / (result_v1$length * 3))
      result_df <- mutate(result_df, d = c / result_v1$totalreads)
      # group_by_sum
      result_df1 <- aggregate(result_df$d, by = list(result_df$gene), sum)
      result_df1 <- as.data.frame(result_df1)
      names(result_df1) <- c("gene", as.character(i))
      system("mkdir sample_gene_abundance")
      write.table(result_df1, file = paste("./sample_gene_abundance/", i, ".csv", sep = ""), sep = ",", quote = FALSE, row.names = FALSE)
      print(i)
    }

    abundance_csv <- list.files(path = "./sample_gene_abundance/", pattern = ".csv")
    n <- length(abundance_csv)
    merge.data <- read.csv(file = paste("./sample_gene_abundance/", abundance_csv[1], sep = ""), header = T, sep = ",")
    merge.data1 <- data.frame(merge.data)
    for (i in 2:n) {
      new.data <- read.csv(file = paste("./sample_gene_abundance/", abundance_csv[i], sep = ""), sep = ",", header = T)
      merge.data1 <- merge(merge.data1, new.data, by = "gene", all = T)
    }
    write.csv(merge.data1, file = paste(output, "sample_abundance.csv", sep = ""), row.names = F)
  }

  # count functional gene drive species of all samples
  if (analysis == "taxonomy") {
    if (PE == TRUE) {
      list <- list.files(path = "./sample_merge/", pattern = method)
    } else {
      list <- list.files(path = "./sample_file/", pattern = method)
    }
    system("mkdir sample_gene_tax")
    for (i in list) {
      if (PE == TRUE) {
        a <- read.table(file = paste("./sample_merge/", i, sep = ""), sep = " ", header = F)
      } else {
        a <- read.table(file = paste("./sample_file/", i, sep = ""), sep = "\t", header = F)
      }
      a <- as.data.frame(a)
      names(a)[1] <- "reads_id"
      names(a)[3] <- "protein_id"
      a <- a[, c(1, 3)] %>% filter(!duplicated(a[, 1]))
      i <- gsub(".diamond", "", i)
      a1 <- mutate(a, "sample" = i)
      sample_gene_tax_pathway <- merge(a1, id_gene_tax_pathway, by = "protein_id")
      write.table(sample_gene_tax_pathway, file = paste("./sample_gene_tax/", i, ".csv", sep = ""), sep = ",", quote = FALSE, row.names = FALSE)
    }

    a <- list.files(path = "./sample_gene_tax/", pattern = ".csv")
    n <- length(a)
    merge.data <- read.csv(file = paste("./sample_gene_tax/", a[1], sep = ""), header = T, sep = ",")
    merge.data1 <- data.frame(merge.data)
    for (i in 2:n) {
      new.data <- read.csv(file = paste("./sample_gene_tax/", a[i], sep = ""), sep = ",", header = T)
      merge.data1 <- merge(merge.data1, new.data, all = T)
    }
    write.csv(merge.data1, file = paste(output, "sample_gene_tax_pathway.csv", sep = ""), row.names = F)
    print(n)
    }
}
