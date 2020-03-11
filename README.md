rSuperCT is the R version of SuperCT, implement the supervised-learning framework for enhanced characterization of single-cell transcriptomic profiles

## Preparation
Before the installation, Python 2.x or 3.x should be available in your system. Type below commands in your R console to check if python is installed. Set initialize = TRUE to initialize Python bindings.
```{r}
install.packages("reticulate")
library(reticulate)
py_available(initialize = TRUE)
```
And, specify your python path if you meant to change python version.
```
use_python(python = 'path/to/your/python', required = TRUE)
py_config() # check configuration of python enviroment.
```
To run rSuperCT smoothly, several modules should be installed. Use either ***pip/pip3 install module_name*** to install if they don't.
```{r eval=FALSE, include=FALSE}
py_module_available(module = 'keras')
py_module_available(module = 'tensorflow')
```
## Installation
To use rSuperCT and load package,
```
devtools::install_github('weilin-genomics/rSuperCT')
library(rSuperCT)
```
## Usage of rSuperCT
Follow below pipelines to predict cell types for scRNA-seq data. For now, only human and mouse are supported.

1) Take pbmc_small dataset in Seurat package as example, the input format for ImportData() function can be Seurat object, (sparse)matrix or data frame in which the rows represent for features and the columns for cells.
```{r}
library(Seurat)
myces <- ImportData(pbmc_small)
myces # show the number of cells and features
```
2) To predict cell types using pre-trained and well-named models published by weilin-genomics/SuperCT/models,
specified model by model parameter and corresponding files will be downloaded and prepared from github to local results.dir.
As prediction done, A column pred_types was saved in meta.data slot of your CellESet object.

```{r}
dir.create('./TestFiles', showWarnings = FALSE)
myces <- PredCellTypes(myces, species = 'human', model = 'models', results.dir = './TestFiles')
```
3) Visualize the predictive results for summarization
```{r}
table(myces@meta.data$pred_types)
plotHist(myces)
```
## References
[1] Xie, P. et.al. (2019). SuperCT: a supervised-learning framework for enhanced characterization of single-cell transcriptomic profiles. https://doi.org/10.1093/nar/gkz116. Nucleic Acids Research.

[2] https://github.com/weilin-genomics/SuperCT
## Contact
If you have any suggestion, questions and bugs report, feel free to contact weilin.baylor@gmail.com.
