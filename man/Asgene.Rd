% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Asgene.R
\name{Asgene}
\alias{Asgene}
\title{AsgeneDB: A curated orthology arsenic metabolism gene database and computational tool for metagenome annotation}
\usage{
Asgene(
  analysis = "abundance",
  workdir = "./",
  method = "diamond",
  toolpath = "./",
  search_parameters = "-e 1e-4 -p 28 --query-cover 80 --id 50",
  seqtype = "nucl",
  filetype = "fasta",
  PE = TRUE,
  output = "./",
  test.data = FALSE
)
}
\arguments{
\item{analysis}{Choose the target you want to analyze,abundance or taxonomy; type = "character"; default = "abundance"}

\item{workdir}{Specify directory for sequence file location;type = "character"; default = "./"}

\item{method}{Specify the database searching tool you plan to use, currently diamond, usearch and blast are supported; type = "character"; default = "diamond"}

\item{toolpath}{Specify directory for searching tool location; type = "character"; default = "./"}

\item{search_parameters}{Define metagenomic comparison parameters; type = "character"; default = "-e 1e-4 -p 28 --query-cover 80 --id 50"}

\item{seqtype}{Specify your sequence type, nucl or prot; type = "character"; default = "nucl"}

\item{filetype}{Specify the extensions of your sequence files, e.g., fastq, fastq.gz, fasta,fasta.gz, fq, fq.gz, fa, fa.gz; type = "character"; default = "fasta"}

\item{PE}{Specifies whether your metagenomic data were PE files, noting that PE files are named "_R1 with extension" or "_R2 with extension" (e.g., XINWEI_R1.fasta and XINWEI_R2.fasta); type = "logical"; default = TRUE}

\item{output}{Specify the directory for the final result file; type = "character"; default = "./"}

\item{test.data}{Example datasets are provided as input and output to help users better understand this package; type = "logical"; default = FALSE}
}
\description{
A manually curated arsenic functional gene database (AsgeneDB) and R package (Asgene package) are developed for rapid and accurate metagenomic analysis.
}
\examples{
library(Asgene)
Asgene(anlysis = "abundance", workdir = "./", method = "diamond", toolpath = "./", search_parameters = "-e 1e-4 -p 28 --query-cover 80 --id 50", seqtype = "nucl", filetype = "fasta", PE = TRUE, out = "./", test.data = FALSE)
Asgene(anlysis = "abundance", workdir = "./", method = "diamond", toolpath = "./", search_parameters = "-e 1e-4 -p 28 --query-cover 80 --id 50",seqtype = "prot", output = "./", test.data = TRUE)
}
