setwd("/Users/bstan/Documents/UW/Courses/BIOST 527")
rm(list = ls())
library(tidyverse)
library(tidyr)
library(tinytex)
library(class)

############
# Question 1
############

max_p = 500
N_list = list("10","100","1000","10000")
d = data.frame(matrix(ncol = length(N_list)+1, nrow = max_p))
colnames(d) <- c("p",N_list)
d$p <- 1:max_p
for(k in 1:length(N_list)){
  med_dist <- NULL
  N=strtoi(N_list[[k]])
  for(i in 1:max_p){
    med_dist<- c(med_dist,(1-0.5^(1/N))^(1/i))
  }
  d[,k+1] <- med_dist
}

plot(
     d[["p"]],d[["10"]], 
     xlab="Number of Dimensions (p)", 
     ylab="Median Distance from Origin to Closest Point", 
     col="blue",
     type="l",
     lty=1,
     main="Median Distance from Origin to Closest Point \n for Uniform Samples in p Dimensions"
     )
lines(d[["p"]],d[["100"]],
      col="black",
      type="l",
      lty=1)
lines(d[["p"]],d[["1000"]],
      col="red",
      type="l",
      lty=1)
lines(d[["p"]],d[["10000"]],
      col="#228B22",
      type="l",
      lty=1)
legend("bottomright", legend=c("N=10","N=100","N=1000","N=10000"),col=c("blue","black","red","#228B22"), lty=c(1,1,1,1), cex=0.8)

############
# Question 2
############

p <- 1:200
L <- .1^(1/p)
plot(
  p,L, 
  xlab="Number of Dimensions (p)", 
  ylab="Length of Side of Hypercube (L(p))", 
  col="blue",
  type="l",
  lty=1,
  main="Length of Side of Smaller Hypercube v. Number of Dimensions (p)"
)

############
# Question 4
############
## Load data
zipcode_train_all <- data.frame(read.table(file="datasets/zip.train",sep=" ", header=FALSE))
zipcode_train <- subset(zipcode_train_all,V1 %in% c(5,8))
zipcode_test_all <- data.frame(read.table(file="datasets/zip.test",sep=" ", header=FALSE))
zipcode_test <- subset(zipcode_test_all,V1 %in% c(5,8))
## Turn label into 1 for 8 and 0 for 5
zipcode_train$V1 <- ifelse(zipcode_train$V1==8,1,0)
zipcode_test$V1 <- ifelse(zipcode_test$V1==8,1,0)

## Check for nulls
count_nulls <- sapply(zipcode_train, function(x) sum(is.na(x)))
count_nulls[count_nulls>0]
drops <- c("V258")
zipcode_train <- zipcode_train[,!(names(zipcode_train) %in% drops)]

## Fit linear regression model and calculate error
lm_fit <- lm(V1 ~ .-V1,data=zipcode_train)
training_predict <- ifelse(predict(lm_fit)<0.5,0,1)
conf_matrix <- table(training_predict,zipcode_train[["V1"]])
ls_train_error = 1-(conf_matrix[1,1]+conf_matrix[2,2])/sum(conf_matrix) # Training accuracy for LS classification
ls_train_error

test_predict <- ifelse(predict(lm_fit,zipcode_test[,-1])<0.5,0,1)
conf_matrix <- table(test_predict,zipcode_test[["V1"]])
ls_test_error = 1-(conf_matrix[1,1]+conf_matrix[2,2])/sum(conf_matrix) # Test accuracy for LS classification
ls_test_error

## Fit KNN and calculate error
set.seed(44)
k_range <- 1:100
accuracy_train <- NULL
for(i in k_range){
  knn_pred <- knn(zipcode_train[,-1],zipcode_train[,-1],zipcode_train[["V1"]],k=i)
  conf_matrix <- table(knn_pred,zipcode_train[["V1"]])
  accuracy_train <- c(accuracy_train,1-(conf_matrix[1,1]+conf_matrix[2,2])/sum(conf_matrix))
}
accuracy_train

accuracy_test <- NULL
for(i in k_range){
  knn_pred <- knn(zipcode_train[,-1],zipcode_test[,-1],zipcode_train[["V1"]],k=i)
  conf_matrix <- table(knn_pred,zipcode_test[["V1"]])
  accuracy_test <- c(accuracy_test,1-(conf_matrix[1,1]+conf_matrix[2,2])/sum(conf_matrix))
}
accuracy_test 

## Plot values of error
plot(
      1/k_range,accuracy_test,
      xlab="Model Complexity in KNN (1/k)",
      ylab="Error",
      col="blue",
      type="l",
      lty=1,
      main="Error of KNN and Least Squared Classifiers v. Model Complexity in KNN",
      ylim = c(0,0.045)
    )
lines(
      1/k_range,accuracy_train,
      col="#228B22",
      type="l",
      lty=1
      )
abline(a=ls_test_error,b=0,col="blue",lty=2)
abline(a=ls_train_error,b=0,col="#228B22",lty=2)
legend("topright",
       legend=c("Training Error (KNN)","Training Error (Least Squares)","Test Error (KNN)","Test Error (Least Squares)"),
       col=c("#228B22","#228B22","blue","blue"),
       lty=c(1,2,1,2), cex=1
       )

## Find the values of k for which KNN outperforms LS
indices_ls_improved = which(accuracy_test<ls_test_error)
accuracy_test[indices_ls_improved]
k_range[indices_ls_improved]

## Find the values of k corresponding to the minimum error
min_indices = which(accuracy_test==min(accuracy_test))
accuracy_test[min_indices]
k_range[min_indices]

############
# Question 5
############
## Define x1, x2 and f
set.seed(10)
x1 <- seq(0, 10, length.out = 100) #runif(1000, min = 0, max = 1)
x2 <- seq(0, 2*pi, length.out = 100) # values used for nonlinear mean model
custom_fxn <- function(x1, x2) {  # Create custom function
  z <- (x1*sin(x2))/(x1+x2) # function used for nonlinear mean model
  #z <- 1.5*x1+4*x2 # function used for linear mean model
  return(z)
}
custom_fxn_y <- function(x1, x2) {  # Create custom function
  z <- (x1*sin(x2))/(x1+x2)+rnorm(1, mean=0, sd=0.2) # function used for nonlinear mean model
  #z <- 1.5*x1+4*x2+rnorm(1, mean=0, sd=1) # function used for linear mean model
  return(z)
}


## Generate matrix of values for use in image()
fxn <- outer(x1, x2, custom_fxn)
image(x1, x2, fxn,main=expression("Plot of f(x"[1]*",x"[2]*")"))

## Generate 100 training and test values from data
x1_train <- sample(x1, 100, replace=TRUE)
x2_train <- sample(x2, 100, replace=TRUE)
x1_test <- sample(x1, 100, replace=TRUE)
x2_test <- sample(x2, 100, replace=TRUE)
train_df <- data_frame(x1_train,x2_train)
colnames(train_df) <- c("x1","x2")
train_df$f <- mapply(custom_fxn_y, train_df$x1, train_df$x2)

test_df <- data_frame(x1_test,x2_test)
colnames(test_df) <- c("x1","x2")
test_df$f <- mapply(custom_fxn_y, test_df$x1, test_df$x2)

## Fit linear regression and compute train/test error
lm_fit <- lm(f ~ x1 + x2,data=train_df)
summary(lm_fit)
ls_train_error <- mean(lm_fit$residuals^2) # MSE for training data

predicted_vals <- predict(lm_fit,test_df)
actuals_preds <- data.frame(cbind(actual=test_df[["f"]], predicted=predicted_vals))
actuals_preds$error <- actuals_preds$actual - actuals_preds$predicted
ls_test_error <- mean(actuals_preds$error^2) # MSE for test data

## Fit KNN and compute train/test error
k_range <- 1:100
accuracy_train <- NULL
for(i in k_range){
  predicted_vals <- FNN::knn.reg(train=select(train_df,c("x1","x2")),test=select(train_df,c("x1","x2")),y=train_df$f,k=i)$pred
  predicted_vals
  actuals_preds <- data.frame(cbind(actual=train_df[["f"]], predicted=predicted_vals))
  actuals_preds$error <- actuals_preds$actual - actuals_preds$predicted
  accuracy_train <- c(accuracy_train,mean(actuals_preds$error^2))
}
accuracy_train

accuracy_test <- NULL
for(i in k_range){
  predicted_vals <- FNN::knn.reg(train=select(train_df,c("x1","x2")),test=select(test_df,c("x1","x2")),y=train_df$f,k=i)$pred
  predicted_vals
  actuals_preds <- data.frame(cbind(actual=test_df[["f"]], predicted=predicted_vals))
  actuals_preds$error <- actuals_preds$actual - actuals_preds$predicted
  accuracy_test <- c(accuracy_test,mean(actuals_preds$error^2))
}
accuracy_test

## Find the values of k for which KNN outperforms LS
indices_ls_improved = which(accuracy_test<ls_test_error)
accuracy_test[indices_ls_improved]
k_range[indices_ls_improved]

## Find the value of k with the lowest MSE
min_indices = which(accuracy_test==min(accuracy_test))
accuracy_test[min_indices]
k_range[min_indices]

## Plot MSE values
plot(
  1/k_range,accuracy_test,
  xlab="Model Complexity in KNN (1/k)",
  ylab="Mean Squared Error (MSE)",
  col="blue",
  type="l",
  lty=1,
  main="MSE of KNN and LS Regression v. Model Complexity in KNN - Linear Mean Model",
  ylim = c(0,0.2)
  #ylim = c(0,10)
)
lines(
  1/k_range,accuracy_train,
  col="#228B22",
  type="l",
  lty=1
)
abline(a=ls_test_error,b=0,col="blue",lty=2)
abline(a=ls_train_error,b=0,col="#228B22",lty=2)
legend("topright",
       legend=c("Training Error (KNN)","Training Error (Least Squares)","Test Error (KNN)","Test Error (Least Squares)"),
       col=c("#228B22","#228B22","blue","blue"),
       lty=c(1,2,1,2), cex=1
)