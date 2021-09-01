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
