library(caret)
library(survival)
library(survminer)

inFile <- 'data/erlotinib_CCLE_PDXE.rds'
df <- readRDS(inFile)
##---- canter and scale data ----------
df$ccle$x <- scale(df$ccle$x)
df$pdx$x <- scale(df$pdx$x)

##-------------------
mod <- readRDS('data/erlotinib_CCLE-CTRPv2_model_TextMining.rds')
prd <- data.frame(time=df$pdx$OS)
prd$pred <- predict(mod, df$pdx$x)
prd$class <- ifelse(prd$pred>=median(prd$pred), 'Sensitive', 'Resistant')
prd$class <- factor(prd$class)

##------censor at 60 days -------
prd$status <- 1
prd$status[prd$time > 60] <- 0
prd$time[prd$time > 60] <- 60

fit <- survfit(Surv(time, status) ~ class, data = prd)

ggsurvplot(fit, pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, risk.table.col = "strata",
           ggtheme = theme_classic(),
           legend.labs = levels(prd$class),
           legend.title="Class",
           palette = c(#"#A73030FF", "#0073C2FF"
             "#DC0000FF", "#20854EFF" #"#0072B5FF"
             #"#E64B35FF" ,"#00A087FF"
             #"#DC0000FF", "#3B3B3BFF"
             ))

