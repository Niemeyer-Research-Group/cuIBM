%plot pressure
clc
clear
% close all
figure
%change these
number = '0';
type = 'u'; %p or u
suffix = 'hat'; %u: 0, star, hat, hatfinal, empty. p: 0, star, empty
view = 'in';

%load data
% caseFolder = '/scratch/src/cuIBM/validation/luo/test/output/'
% caseFolder = '/scratch/src/cuIBM/validation/cylinder/Re40/output/';
caseFolder = '/scratch/src/cuIBM/validation/osc/gh/output/';
path = strcat(caseFolder,number,type,suffix,'.csv');
ghostpath = strcat(caseFolder,number,'ghost',type,'.csv');
delim = '\t';
u = dlmread(path,delim,1,0);
test = u;
N = dlmread(ghostpath,delim,1,0);

% manipulate inside/outside
for i =1:length(u(:,1))
    for j = 1:length(u(1,:))
        if strcmp(view,'out')
            if N(i,j)~=-1
                u(i,j) = nan;
            end
        elseif strcmp(view,'in')
            if N(i,j)==0
                u(i,j) = nan;
            end
        end
    end
end
%plot area round body
midy = round(length(u(:,1))/2);
midx = round(length(u(1,:))/2);
surf(u((midy-50):(midy+50),(midx-50):(midx+50)))
% surf(u(60:100,1:127)), hold on
% surf(u)
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
ghostpath = strcat(caseFolder,number,'ghost',type,'.csv');
delim = '\t';
u = dlmread(path,delim,1,0);
test = u;
N = dlmread(ghostpath,delim,1,0);

body = dlmread(strcat('/scratch/src/cuIBM/validation/osc/flow/midPosition'),delim,1,0);
r=0.5;
teta=-pi:0.01:pi;
midx = body(str2num(number)-1,2);
z = body(str2num(number)-1,4);
x=r*cos(teta) + midx;
y=r*sin(teta);

% manipulate inside/outside
for i =1:length(u(:,1))
    for j = 1:length(u(1,:))
        if strcmp(view,'out')
            if N(i,j)~=-1
                u(i,j) = nan;
            end
        elseif strcmp(view,'in')
            if N(i,j)==-1
                u(i,j) = nan;
            end
        end
    end
end
%plot
h = 0.03125;
X = linspace(-2.0+h,2.0-h,127);
Y = linspace(-2.0+h/2,2.0-h/2,128);
surf(X,Y,u(1:128,1:127)), hold on
fill3( x,y,zeros(1,numel(x))+z,[0 0 0] )
title(type)
xlabel('x')
ylabel('y')
zlabel('z')

%% plot hybrid and ghost velocities 
clc
clear
close all

number = '100';
type = 'u'; %p or u

%load data
caseFolder = '/scratch/src/cuIBM/validation/osc/flow/output/';
path = strcat(caseFolder,number,type,'.csv');
ghostpath = strcat(caseFolder,number,'ghostu.csv');
hybridpath = strcat(caseFolder,number,'hybridu.csv');
delim = '\t';
u = dlmread(path,delim,1,0);
ghost = dlmread(ghostpath,delim,1,0);
hybrid = dlmread(hybridpath,delim,1,0);

%setup body edge
body = dlmread(strcat('/scratch/src/cuIBM/validation/osc/flow/midPosition'),delim,1,0);
r=0.5;
teta=-pi:0.01:pi;
midx = body(str2num(number),2);
z = body(str2num(number),4);
x=r*cos(teta) + midx;
y=r*sin(teta);

%setup plot
h = 0.03125;
X = linspace(-2.0+h,2.0-h,127);
Y = linspace(-2.0+h/2,2.0-h/2,128);

%plot
hold on

for i =1:length(u(:,1))
    for j = 1:length(u(1,:))
        if ghost(j,i)>0% || hybrid(i,j)>0
            plot1 = plot3(X(i),Y(j),u(j,i),'ro');
        end
        if hybrid(j,i)>0
            plot2 = plot3(X(i),Y(j),u(j,i),'bo');
        end
    end
end
plot3(x,y,zeros(1,numel(x))+z)
[XX,YY] = meshgrid(X,Y);
for i=44:84
    plot([XX(1,i) XX(end,i)],[YY(1,i) YY(end,i)],'color',[0 0 0] + 0.80) %vert lines
    plot([XX(i,1) XX(i,end)],[YY(i,1) YY(i,end)],'color',[0 0 0] + 0.80) %horz lines
end
% plot(,,'color',[0 0 0] + 0.80)
xlabel('x')
ylabel('y')
zlabel('z')
axis([-0.6 0.6 -0.6 0.6])
grid off
legend([plot1 plot2],'Ghost','Hybrid')

%% plot pressure outside body and first value inside
clc
clear
close all

number = '100';

%load data
caseFolder = '/scratch/src/cuIBM/validation/osc/flow/output/';
% caseFolder = '/scratch/src/cuIBM/validation/osc/static/output/';
path = strcat(caseFolder,number,'p.csv');
ghostpath = strcat(caseFolder,number,'ghostp.csv');
delim = '\t';
p = dlmread(path,delim,1,0);
ghost = dlmread(ghostpath,delim,1,0);

% manipulate inside/outside
for i =1:length(p(:,1))
    for j = 1:length(p(1,:))
        if ghost(i,j)==0
            p(i,j) = nan;
        end
    end
end

% h = 0.005;
% s = 0.4;
% X = linspace(-s+h/2,s-h/2,160);
% s=0.2;
% Y = linspace(-s+h/2,s-h/2,80);
h = 0.03125;
X = linspace(-2.0+h/2,2.0-h/2,128);
Y = linspace(-2.0+h/2,2.0-h/2,128);
%plot area round body
surf(X(44:88),Y(44:88),p(44:88,44:88)), hold on
% surf(X,Y,p(73:152,81:240))
xlabel('x')
ylabel('y')
zlabel('z')