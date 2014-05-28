\name{adjust}
\alias{adjust}
\title{Adjust the dimensions of a matrix to the dimensions of another matrix}
\description{Adjust the dimensions of a matrix to the dimensions of another matrix.}
\usage{
adjust(source, target, remove = TRUE, add = TRUE)
}
\details{
An adjacency matrix (the \code{source} matrix) is compared to another adjacency matrix (the \code{target} matrix) by matching the row or column labels. If the target matrix contains rows/columns which are not present in the source matrix, new rows and columns with the corresponding labels and NA values in the cells are inserted into the source matrix. If the source matrix contains rows/columns which are not present in the target matrix, these rows and columns are removed from the source matrix.
}
\arguments{
\item{source}{ A matrix object which should be adjusted. The horizontal and vertical dimensions of the matrix must be equal. }
\item{target}{ A matrix object to which the source matrix is compared in terms of labels. }
\item{remove}{ Should rows and columns that are not present in the target matrix be removed? }
\item{add}{ Should rows and columns that are present in the target matrix but not in the source matrix be added to the source matrix? }
}
\seealso{
\link{xergm-package} \link{handleMissings} \link{preprocess}
}
