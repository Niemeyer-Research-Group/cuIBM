%plot pressure
clc
clear
close all
% figure
%change these
number = '100';
type = 'u'; %p or u
suffix = ''; %u: 0, star, hat, hatfinal, empty. p: 0, star, empty
view = 'out';

%load data
% caseFolder = '/scratch/src/cuIBM/validation/luo/test/output/'
% caseFolder = '/scratch/src/cuIBM/validation/cylinder/Re40/output/';
caseFolder = '/scratch/src/cuIBM/validation/osc/flow/output/';
path = strcat(caseFolder,number,type,suffix,'.csv');
tagspath = strcat(caseFolder,number,'ghost',type,'.csv');
delim = '\t';
M = dlmread(path,delim,1,0);
test = M;
N = dlmread(tagspath,delim,1,0);

% manipulate inside/outside
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
% surf(M((midy-50):(midy+50),(midx-50):(midx+50)))
surf(M), hold on
title(strcat(type,suffix))
xlabel('x')
ylabel('y')
zlabel('z')

%%
clc
clear
close all

number = '100';
type = 'u'; %p or u
view = 'in';

%load data
caseFolder = '/scratch/src/cuIBM/validation/osc/flow/output/';
path = strcat(caseFolder,number,type,'.csv');
tagspath = strcat(caseFolder,number,'ghost',type,'.csv');
delim = '\t';
M = dlmread(path,delim,1,0);
test = M;
N = dlmread(tagspath,delim,1,0);

body = dlmread(strcat('/scratch/src/cuIBM/validation/osc/flow/midPosition'),delim,1,0);
r=0.5;
teta=-pi:0.01:pi;
midx = body(str2num(number)-1,2);
z = body(str2num(number)-1,4);
x=r*cos(teta) + midx;
y=r*sin(teta);

% manipulate inside/outside
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
%plot
h = 0.03125;
X = linspace(-2.0+h,2.0-h,127);
Y = linspace(-2.0+h,2.0-h,127);
surf(X,Y,M(1:127,1:127)), hold on
fill3( x,y,zeros(1,numel(x))+z,[0 0 0] )
title(type)
xlabel('x')
ylabel('y')
zlabel('z')