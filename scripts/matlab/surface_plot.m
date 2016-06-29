%plot pressure
clc
clear
% close all
figure
%change these
number = '199';
type = 'u'; %p or u
suffix = '0'; %u: 0, star, hat, hatfinal, empty. p: 0, star, empty
view = 'outt';

%load data
path = strcat('/scratch/src/cuIBM/validation/luo/test/output/',number,type,suffix,'.csv');
tagspath = strcat('/scratch/src/cuIBM/validation/luo/test/output/',number,'ghost',type,'.csv');
delim = '\t';
M = dlmread(path,delim,1,0);
test = M;
N = dlmread(tagspath,delim,1,0);
%manipulate inside/outside
for i =1:length(M(:,1))
    for j = 1:length(M(1,:))
        if strcmp(view,'out')
            if N(i,j)~=-1
                M(i,j) = nan;
            end
        elseif strcmp(view,'in')
            if N(i,j)==-1
                M(i,j) = nan;
            end
        end
    end
end
%plot area round body
midy = round(length(M(:,1))/2);
midx = round(length(M(1,:))/2);
surf(M((midy-50):(midy+50),(midx-50):(midx+50)))
%surf(M)
title(strcat(type,suffix))
xlabel('x')
ylabel('y')
zlabel('z')