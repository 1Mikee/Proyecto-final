
library(DynamicCancerDriverKM)
library(e1071)
library(caret) 
library(dplyr)  
library(pROC)   
library(tidyverse)
library(class) 
library(rpart) 
library(glmnet)


# Cargar datos
view(DynamicCancerDriverKM::BRCA_normal)
view(DynamicCancerDriverKM::BRCA_PT)

load("C:\\Users\\home\\Desktop\\Electiva\\RStudio\\data\\geneScore.rdata")


normal_pt <- rbind(BRCA_normal, BRCA_PT)

df <- normal_pt[, !(names(normal_pt) %in% c("barcode", "bcr_patient_barcode", "bcr_sample_barcode", "vital_status", "days_to_death", "treatments_radiation_treatment_or_therapy"))]

any(is.na(df))

muestras <- as.matrix(df[, -1])

umbral <- 0.0002 * max(muestras)

genes_expresados <- muestras > umbral

verdaderos_por_gen <- colSums(genes_expresados)

umbral_eliminar_columna <- nrow(muestras) * 0.2

columnas_a_conservar <- which(verdaderos_por_gen >= umbral_eliminar_columna)

filtro_genes <- df[, c(1, columnas_a_conservar + 1)] 

geneScore <- prub$features

# Obtener los nombres de genes en filtered_data
genes <- colnames(filtro_genes)[-1]  # Excluir la columna "sample_type"

# Encontrar los genes comunes
genes_comunes <- intersect(geneScore, genes)

genes_comunes <- prub[geneScore %in% genes_comunes, ]

gene_sorted <-  genes_comunes %>% arrange(desc(score))
top_genes <- gene_sorted[1:100, ]


top_100 <- top_genes$features

y <- filtro_genes$sample_type

X <- filtro_genes[, top_100]

y <- as.factor(y)


set.seed(123)  # Semilla para reproducibilidad
trainIndex <- createDataPartition(y, p = 0.8, list = FALSE) # se puede modificar la divicion del modelo para ver otros resultados p = 0.7 
train_data <- X[trainIndex, ]
test_data <- X[-trainIndex, ]
train_labels <- y[trainIndex]
test_labels <- y[-trainIndex]

# Modelo svm con geneScore PIK3R1

model <- svm(train_labels ~ ., data = cbind(train_data, train_labels), kernel = "linear")

predictions <- predict(model, newdata = cbind(test_data, test_labels))

confusionMatrix(predictions, test_labels)

roc_curve <- roc(as.numeric(predictions), as.numeric(test_labels))
roc_curve


# Modelo de Regresión Logistica con geneScore PIK3R1

logistic_model <- cv.glmnet(as.matrix(train_data), train_labels, family = "binomial")

predictions <- predict(logistic_model, newx = as.matrix(test_data), s = "lambda.1se", type = "response")

predicted_labels <- as.factor(ifelse(predictions > 0.5, levels(y)[2], levels(y)[1]))

conf_matrix <- confusionMatrix(predicted_labels, test_labels)
conf_matrix

precision <- posPredValue(predicted_labels, test_labels)
paste("Precisión del modelo de regresión logística:", precision)

roc_curve <- roc(as.numeric(predicted_labels), as.numeric(test_labels))
roc_curve

normalized_train_data <- scale(train_data)
normalized_test_data <- scale(test_data)

knn_model <- knn(train = normalized_train_data, test = normalized_test_data, cl = train_labels, k = 5)

knn_conf_matrix <- confusionMatrix(knn_model, test_labels)
knn_conf_matrix

# Modelo de Árboles de decisiones con geneScore PIK3R1

tree_model <- rpart(train_labels ~ ., data = train_data, method = "class")

plot(tree_model)
text(tree_model, pretty = 0)

tree_predictions <- predict(tree_model, newdata = test_data, type = "class")

tree_conf_matrix <- confusionMatrix(tree_predictions, test_labels)
print(tree_conf_matrix)

