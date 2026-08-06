library(readr)
library(dplyr)
library(survival)
library(survminer)

expr_annotated <- readRDS(
  "data/processed/expr_annotated.rds"
)

pheno_clean <- readRDS(
  "data/processed/pheno_clean.rds"
)
#Kaplan–Meier survival analysis
#Load the survival data

surv <- read_tsv("data/raw/TCGA-BRCA.survival.tsv.gz")
dim(surv)

head(surv)

colnames(surv)
#Merge OAS3 with survival
oas3_surv <- oas3_long %>%
  left_join(surv, by = "sample")
#removing missing survival values
oas3_surv <- oas3_surv %>%
  filter(!is.na(OS.time),
         !is.na(OS))
dim(oas3_surv)

median_expr <- median(oas3_surv$expression)

median_expr

oas3_surv <- oas3_surv %>%
  mutate(
    OAS3_group = ifelse(
      expression >= median_expr,
      "High",
      "Low"
    )
  )
table(oas3_surv$OAS3_group)

install.packages(c("survival", "survminer"))
library(survival)
library(survminer)
#fit the kaplan-meier  model
fit <- survfit(
  Surv(OS.time, OS) ~ OAS3_group,
  data = oas3_surv
)
#plot the kaplan-meier curve
ggsurvplot(
  fit,
  data = oas3_surv,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  xlab = "Time (days)",
  ylab = "Overall survival probability",
  legend.title = "OAS3",
  legend.labs = c("High", "Low")
)

ggsave(
  "figures/Figure2_OAS3_Expression_Impact_Overall_Survival_BRCA.png",
  width = 6,
  height = 5,
  dpi = 300
)
cox <- coxph(
  Surv(OS.time, OS) ~ OAS3_group,
  data = oas3_surv
)

summary(cox)

saveRDS(oas3_surv,
        "data/processed/oas3_surv.rds")
