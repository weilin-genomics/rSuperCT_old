#' Make predictions of cell types based on the single-cell RNA-seq digital expression
#' profiles using a supervised classifier, SuperCT
#' @param object CellESet object.
#' @param species Either 'human' or 'mouse' for now.
#' @param model Choose a supported model. i.e 'CellTypes', 'CellStates'.
#' Please refer to https://github.com/weilin-genomics/SuperCT for more details.
#' @param use.cells Character vector specify cells to make prediction for.
#' @param results.dir Specify directory to save downloaded required files for the prediction.
#' @importFrom reticulate import
#' @importFrom utils download.file read.csv
#' @return
#' Predicted cell identities saved in \code{object@meta.data[['pred_types']]}
#' @export
#' @references
#' Xie Peng and Gao Mingxuan (2019) \emph{SuperCT: a supervised-learning framework for
#' enhanced characterization of single-cell transcriptomic profiles},
#' \url{https://doi.org/10.1093/nar/gkz116} \emph{Nucleic Acids Research}
#' @examples
#' \dontrun{
#' myces <- PredCellTypes(myces, species = 'human', model = 'models', results.dir = '~/TestFiles')
#' }
PredCellTypes <- function(object, species = 'human', model = 'CellTypes', use.cells = object@use.cells, results.dir = getwd()){
  classOfObj <- class(object)
  if(classOfObj != 'CellESet')
    stop('A object of class CellESet expected.')
  use.cells <- intersect(x = use.cells, y = colnames(x = object@raw.data))
  if(length(x = use.cells) == 0)
    stop('No available cell names.')
  if(! species %in% c('human', 'mouse'))
    stop('Unrecognized species. Specify as follows: human, mouse.')

  # prepare required files used in prediction
  repo.dir <- paste0(results.dir, '/', model)
  dir.create(path = repo.dir, showWarnings = FALSE)
  all_files <- dir(path = repo.dir, full.names = TRUE)
  model.file <- grep(pattern = 'model.h5', x = all_files, value = TRUE)
  types.file <- grep(pattern = 'type.csv', x = all_files, value = TRUE)
  genes.file <- grep(pattern = 'genes.csv', x = all_files, value = TRUE)
  if(length(x = model.file) == 0){
    model.file <- paste0(repo.dir, '/model.h5')
    download.file(url = paste0('https://raw.githubusercontent.com/weilin-genomics/SuperCT/master/',
                               model, '/v1_model.h5'),
                  destfile = model.file)
  }
  if(length(x = types.file) == 0){
    types.file <- paste0(repo.dir, '/type.csv')
    download.file(url = paste0('https://raw.githubusercontent.com/weilin-genomics/SuperCT/master/',
                               model, '/v1_id2type.csv'),
                  destfile = types.file)
  }
  if(length(x = genes.file) == 0){
    genes.file <- paste0(repo.dir, '/genes.csv')
    download.file(url = 'https://raw.githubusercontent.com/weilin-genomics/SuperCT/master/inthomgenes.csv',
                  destfile = genes.file)
  }
  genes <- read.csv(file = genes.file, header = FALSE, stringsAsFactors = FALSE)
  genes_used <- switch(EXPR = species,
                       human = genes$V1,
                       mouse = genes$V2)
  celltypes <- read.csv(file = types.file, header = TRUE, stringsAsFactors = FALSE)

  # prepare the expression profile and complete it if there's genes not found
  genes_inter <- intersect(x = genes_used, y = rownames(x = object@raw.data))
  genes_supp <- setdiff(x = genes_used, y = genes_inter)
  data1 <- as.matrix(object@raw.data[genes_inter, use.cells, drop = FALSE])
  dimns <- list(genes_supp, use.cells)
  data2 <- matrix(data = 0, nrow = length(x = genes_supp), ncol = length(x = use.cells), dimnames = dimns)
  data <- t(rbind(data1, data2))
  data <- data[, genes_used]
  data[data > 0] <- 1

  # use trained model to predict cell types
  keras.models <- import(module = 'keras.models')
  mm <- keras.models$load_model(model.file)
  pred_res <- mm$predict(x = data)

  # get identity of each cell
  pred_types <- celltypes[max.col(m = pred_res), 'celltype']
  object@meta.data$pred_types <- pred_types

  return(object)
}
