% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prediction.R
\name{PredCellTypes}
\alias{PredCellTypes}
\title{Make predictions of cell types based on the single-cell RNA-seq digital expression
profiles using a supervised classifier, SuperCT}
\usage{
PredCellTypes(
  object,
  species = "human",
  model = "CellTypes",
  use.cells = object@use.cells,
  results.dir = getwd()
)
}
\arguments{
\item{object}{CellESet object.}

\item{species}{Either 'human' or 'mouse' for now.}

\item{model}{Choose a supported model. i.e 'CellTypes', 'CellStates'.
Please refer to https://github.com/weilin-genomics/SuperCT for more details.}

\item{use.cells}{Character vector specify cells to make prediction for.}

\item{results.dir}{Specify directory to save downloaded required files for the prediction.}
}
\value{
Predicted cell identities saved in \code{object@meta.data[['pred_types']]}
}
\description{
Make predictions of cell types based on the single-cell RNA-seq digital expression
profiles using a supervised classifier, SuperCT
}
\examples{
\dontrun{
myces <- PredCellTypes(myces, species = 'human', model = 'models', results.dir = '~/TestFiles')
}
}
\references{
Xie Peng and Gao Mingxuan (2019) \emph{SuperCT: a supervised-learning framework for
enhanced characterization of single-cell transcriptomic profiles},
\url{https://doi.org/10.1093/nar/gkz116} \emph{Nucleic Acids Research}
}
