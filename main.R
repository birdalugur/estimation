setwd("/home/ugur/r-projects/estimation/")
library(pracma)
library(Matrix)
library(RSpectra)
library(rmatio)

source('center.R')
source('kalman_filter_diag.R')
source('kalman_smoother_diag.R')
source('kalman_update_diag.R')
source('ricSW.R')
source('smooth_update.R')
source('FactorExtraction.R')


data<-read.mat('testdata.mat')
x <- data$X
p <- data$p
q <- data$q
r <- data$r



result <- FactorExtraction(x,q,r,p)

F <- result$F
VF <- result$VF
A <- result$A
C <- result$C


