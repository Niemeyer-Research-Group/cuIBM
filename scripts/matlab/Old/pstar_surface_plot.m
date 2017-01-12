%plot pressure
clc
clear all
% close all

%change these
number = '0';
view = 'out';
figure
%load data
path = '/scratch/src/cuIBM/validation/luo/test/output/';
pstar_path = strcat(path,number,'pstar.csv');
p_path = strcat(path,number,'p0.csv');
tags_path = strcat(path,number,'hybridp.csv');
ghostTags_path = strcat(path,number,'ghostp.csv');
delim = '\t';
pstar = dlmread(pstar_path,delim,1,0);
p = dlmread(p_path,delim,1,0);
tags = dlmread(tags_path,delim,1,0);
ghost = dlmread(ghostTags_path,delim,1,0);
M = zeros(length(p(:,1)),length(p(1,:)));
%manipulate inside/outside
for i =1:length(M(:,1))
    for j = 1:length(M(1,:))
        if tags(i,j)==-1
            M(i,j) = p(i,j);
        elseif tags(i,j)>0
            M(i,j) = pstar(i,j);
        end
        if ghost(i,j) ~=-1
            M(i,j) = nan;
        end
    end
end
%plot area round body
midy = round(length(M(:,1))/2);
midx = round(length(M(1,:))/2);
surf(M((midy-50):(midy+50),(midx-50):(midx+50)))
%surf(M)
title('P+pstar')
xlabel('x')
ylabel('y')