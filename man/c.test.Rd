\name{c.test}
\alias{c.test}
\title{c.test}
\usage{
c.test(matr,cla,cor=FALSE,net=NULL,nam=NULL)
}
\description{
node change and edge change test for a given gene set.
}
\arguments{
  \item{matr}{a numeric matrix with rows represent genes and columns represent samples}

  \item{cla}{a vetor presents labels of samples with two different values}

  \item{cor}{a logical indicating whether the change of correlation (edge) between genes should be considered. If TRUE, the sample size should be sufficient}

  \item{net}{a matrix with two columns of character indicating the background network, each row presents node names of an edge in the network. When net is set as NULL, ever correlation (edge) between genes will be considered if cor is TRUE}

  \item{nam}{a vector of characters presents the gene names of the input matrix, with the same length as the rows of input matrix. If cor is set as false or net is set as NULL, this parameter can be set as NULL}
}

\value{
vector
}

\examples{
##generating a simulated dataset
a<-rnorm(500,1,1)
a<-matrix(a,5,100)
b<-rnorm(500,1.3,1)
b<-matrix(b,5,100)
matr<-cbind(a,b)
cla<-c(rep(1,100),rep(2,100))

c.test(matr,cla,cor=TRUE)
##testing with a background pathway
nam = c("TP53","ZNF1","ZNF2","ZNF3","ZNF7")
net = cbind(c("TP53","ZNF1","ZNF2","ZNF3"),c("ZNF1","ZNF2","ZNF3","ZNF7"))

c.test(matr,cla,cor=TRUE,nam=nam,net=net)
}

\author{
Yiqun Li
}

