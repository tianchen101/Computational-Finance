clear;
clc;
fid = fopen('crsp_index_nominal.dat','r');
C1 = textscan(fid,'%f%f');
fclose(fid);
date = C1{1,1};
index = C1{1,2};
% so the log return is that 
index_size = size(index,1);
log_returns = log(index(2:index_size)) - log(index(1:(index_size-1)));
log_returns = log_returns - mean(log_returns); % normalize the data
log_returns = log_returns/std(log_returns); % normalize the data


[f,x] =hist(log_returns,200);
g=1/sqrt(2*pi)*exp(-0.5*x.^2);%# pdf of the normal distribution
dx = diff(x(1:2));
bar(x,f/sum(f*dx));
hold on
plot(x,g,'r'); hold off

sigma_1 = 1;
