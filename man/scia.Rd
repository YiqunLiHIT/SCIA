\name{scia}
\alias{scia}
\title{scia}
\usage{
scia(matr, cla, nam, cor=FALSE, net=HPRD, set1=0.02, set2=0.02,pvalue=0.05)
}
\description{
The establishment of CSSPN (Conditional Specific Shortest Pathway Network) with gene expression profile. Output is a matrix with two columns of character indicating the CSSPN, each row presents node names of an edge in the network. The nodes of CSSPN can be used for functional enrichment analysis.
}
\arguments{
  \item{matr}{a numeric matrix with rows represent genes and columns represent samples}

  \item{cla}{a vetor presents labels of samples with two different values}
 
  \item{cor}{a logical indicating whether the change of correlation (edge) between genes should be considered. If TRUE, the sample size should be sufficient}

  \item{net}{a matrix with two columns of character indicating the background network, each row presents node names of an edge in the network}

  \item{nam}{a vector of characters presents the gene names of the input matrix, with the same length as the rows of input matrix}

  \item{set1}{a vector of characters presents the initial gene set of CSSPN or a numeric value presents the proportion of the most significant differential expressed genes by t-test}

  \item{set1}{a vector of characters presents the terminal gene set of CSSPN or a numeric value presents the proportion of the most significant differential expressed genes by t-test}

\value{
matrix
}

\examples{
  data(data,package="SCIA")
  data(HPRD,package="SCIA")
  cla=c(rep(1,51),rep(2,51))
  scia(matr=data, cla=cla, nam=row.names(data), cor=FALSE, net=HPRD, set1=0.02,set2=0.0001,pvalue=0.05,perm=200)
}
\author{
Yiqun Li
}
