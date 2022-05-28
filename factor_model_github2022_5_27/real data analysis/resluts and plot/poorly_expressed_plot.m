%% mean-std plot: poorly-expressed and well-expressed markers
clear;clc;close all
xe = readtable('E:\RR\SAVE\graduate\pfa\R\Ebench.csv');%% change with your pathway
xe = cell2mat(table2cell(xe));
xv = readtable('E:\RR\SAVE\graduate\pfa\R\Vbench.csv');%% change with your pathway
xv = cell2mat(table2cell(xv));

scatter([mean(xv,2);mean(xe,2)],[sqrt(diag(cov(xv')));sqrt(diag(cov(xe')))])