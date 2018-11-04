Here, we provide scripts used to retrieve and preprocess the data for the GTEx application.
Gene expression data can be obtained from [recount](https://jhubiostatistics.shinyapps.io/recount/).
The age information, however, requires dbGAP access.
The script `getData.Rmd` retrieves the data from recount and the sample information. Pre-processing of the counts is done in `preprocessData.R`. The design matrix consisting of the principal components on each tissue and the response as input for the regression is prepared in `prepareData4PCRegression.Rmd`.
