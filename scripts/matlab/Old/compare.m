%compares two csv files
clear all
close all
clc

filename1 = '/scratch/src/cuIBM-FSI/validation/cylinder/Re3000/output/out1.csv';
filename2 = '/scratch/src/cuIBM-FSI/validation/cylinder/Re3000/output/out2.csv';
file1 = csvread(filename1);
file2 = csvread(filename2);

if (length(file1) ~= length(file2))
    1
else
    0
end

for i = 1:length(file1)
    if abs(file1(i) - file2(i)) > 0.00000000001
        i
        file1(i)
        file2(i)
    end
end

2