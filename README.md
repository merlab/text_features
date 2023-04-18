# Text-Mining Based Feature Selection for Anticancer Drug Response Prediction 

## Abstract

Motivation: Predicting anticancer treatment response from baseline genomic data is a
critical obstacle in personalized medicine. Machine learning methods are
commonly used for predicting drug response from gene expression data. In the
process of constructing these machine learning models, one of the most
significant challenges is identifying appropriate features among a massive
number of genes.

Results: In this study, we utilize features (genes) extracted using the text-mining of
scientific literatures. Using two independent cancer pharmacogenomic datasets,
we demonstrate that text-mining-based features outperform traditional feature
selection techniques in machine learning tasks. In addition, our analysis
reveals that text-mining feature-based machine learning models trained on in
vitro data also perform well when predicting the response of in vivo cancer
models. Our results demonstrate that text-mining based feature selection is an
easy to implement approach that is suitable for building machine learning
models for anticancer drug response prediction.

## Descriptions of some scripts in /R

For each argument, refer to "Reference for valid arguments" table below, for all supported inputs.  

### /R/init_train.R  
**Purpose:** call functions in train_functions.R to produce models  
**Call with arguments:** pSet, method, problem, drugname (optional)  
**Example:** rScript init_train.R CCLE glmnet regression  
**Note:** If drugname is specified, model will only be trained for that drug. Otherwise, model will be trained for all drugs in /R/genes  

### /R/crossvalidate.R  
**Purpose:** use models produced by init_train.R to test on another dataset  
**Call with arguments:** pSet, method, problem, drugname (optional)  
**Example:** rScript crossvalidate.R GDSC glmnet regression  
**Note:** If drugname is specified, model will only be tested for that drug. Otherwise, model will be tested for all drugs in /R/genes  

### /R/plot_data.R  
**Purpose:** produces box and whisker plot from train or test results   
**Call with arguments:** pSet, method, problem, drugname, metric, type  
**Example:** rScript plot_data.R GDSC glmnet regression Adavosertib COR test

## Reference for valid arguments  
| Argument  | Valid inputs |
| ------------- | ------------- |
| pSet  | GDSC, CCLE  |
| method  | glmnet, rf  |
| problem | regression, class  |
| metric | Accuracy, BalancedAccuracy, AUC, COR  |
| type | train, test, all  |

## Contact:

Arvind Mer: [Email](mailto:amer@uottawa.ca)
