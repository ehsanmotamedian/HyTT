# MARS model fitting for sequence-derived translational burden (P x L)
# Requires: earth, readxl, caret

# install.packages(c("earth", "readxl", "caret"))

library(earth)
library(readxl)
library(caret)

# ---- Paths (relative to this script's location) ----
data_dir   <- "mars_model/data"
output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE)

input_file <- file.path(data_dir, "training_data.xlsx")

# ---- Load and clean data ----
df <- read_excel(input_file)

if (!"ProteinExp" %in% names(df)) {
  stop("Target column 'ProteinExp' not found -- check the column name in the input file.")
}

df$ProteinExp <- as.numeric(df$ProteinExp)
df <- na.omit(df)

# Use all numeric columns except the target as predictors
numeric_cols   <- sapply(df, is.numeric)
predictor_cols <- names(df)[numeric_cols & names(df) != "ProteinExp"]

# ---- Train / test split (80 / 20) ----
set.seed(123)
trainIndex <- createDataPartition(df$ProteinExp, p = 0.8, list = FALSE)
train_df   <- df[trainIndex, ]
test_df    <- df[-trainIndex, ]

# ---- Fit MARS model ----
# 10-fold cross-validation repeated 30 times, backward pruning
mars_model <- earth(
  ProteinExp ~ .,
  data    = train_df[, c("ProteinExp", predictor_cols)],
  degree  = 1,
  nk      = 100,
  nprune  = 30,
  penalty = 10,
  nfold   = 10,
  ncross  = 30,
  pmethod = "backward"
)

# ---- Predictions and performance ----
pred_train <- as.vector(predict(mars_model, train_df[, predictor_cols]))
pred_test  <- as.vector(predict(mars_model, test_df[, predictor_cols]))

res_train <- postResample(pred = pred_train, obs = train_df$ProteinExp)
res_test  <- postResample(pred = pred_test, obs = test_df$ProteinExp)

# ---- Key results ----
cat("=== Model formula ===\n")
cat(format(mars_model, style = "pmax"), "\n")

cat("\n=== Variable importance ===\n")
print(evimp(mars_model))

cat("\n=== Training set ===\n")
cat("R2 =", round(res_train["Rsquared"], 3),
    " | RMSE =", round(res_train["RMSE"], 4), "\n")

cat("\n=== Held-out test set ===\n")
cat("R2 =", round(res_test["Rsquared"], 3),
    " | RMSE =", round(res_test["RMSE"], 4), "\n")

# ---- Parity plot ----
png(file.path(output_dir, "ProteinExp_Observed_vs_Predicted_test.png"),
    width = 900, height = 700, res = 120)

plot(test_df$ProteinExp, pred_test,
     xlab = "Observed ProteinExp (test set)",
     ylab = "Predicted ProteinExp (MARS)",
     main = paste("Test set: R2 =", round(res_test["Rsquared"], 3),
                  "| RMSE =", round(res_test["RMSE"], 4)),
     pch = 19, col = "blue", cex = 1.2)
abline(0, 1, col = "red", lwd = 2, lty = 2)
grid()
dev.off()
