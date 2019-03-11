%% Preprocess
load('coal_mine_disasters.mat')
histogram(T,50)
xlabel("Year")
ylabel("Number")
title("\tau")
%% Parameters
close all;
psi = 3;
bp = 4;
d = bp + 1;
rho = 0.01;
%% Hybrid MCMC

N = 20000;
burnIn = 5000;
M = N + burnIn;
% t
t = zeros(d+1,M);
t(1,:) = 1658;
t(end,:) = 1980;
% init t
t(:,1) = linspace(1658,1980,d+1);
% theta
theta = zeros(1,M);
% lambda
lambda = zeros(d,M);
% n(tau)
nTau = zeros(d,M);

%Gibbs
for i = 1:M-1
    theta(i) = gamrnd(2,1/psi);
    for j = 1:d
        nTau(j,i) = sum(t(j,i)<=T&T<t(j+1,i));
    end
    lambda(:,i) = gamrnd(nTau(:,i)+2,1./(diff(t(:,i))+theta(i)));
    %MH for t
    for k = 2:d
        R = rho*(t(k+1,i)-t(k-1,i));
        epsilon = -R + 2*R*rand;
        tkNew = t(k,i) + epsilon;
        tNew = t(:,i);
        tkNew;
        nTauNew = [];
        if t(k-1,i)<tkNew && tkNew < t(k+1,i)
            for l = 1:d
                nTauNew(l) = sum(tNew(l)<=T&T<tNew(l+1)) ;
            end
            alpha = min(1,...
            prod(lambda(:,i).^nTauNew(:).*exp(-lambda(:,i).*diff(tNew(:))).*diff(tNew))/...
            prod(lambda(:,i).^nTau(:,i).*exp(-lambda(:,i).*diff(t(:,i))).*diff(t(:,i))));
            alpha;
            if rand < alpha
                t(k,i+1) = tkNew;
            else
                t(k,i+1) = t(k,i);
            end
        else
            t(k,i+1) = t(k,i);
        end
    end
end

tau = ceil(mean(t(:,burnIn:M),2));
%% figures for report. (b)(c)
% lambda
figure
length = size(lambda,1);
for i = 1:length
    subplot(length,1,i)
    histfit(lambda(i,burnIn:end),50,'gamma');
    str = sprintf('\\lambda_%d Distribution', i);
    title(str) 
end
xlabel("\lambda")
savename = sprintf('lambda_Breakpoints%drho%.3fpsi3.png',bp,rho);
saveas(gcf,char(savename))
% theta
figure
histfit(theta(burnIn:end),50,'gamma');
xlabel("\theta")
ylabel("Numbers")
title("\theta Distribution")
savename = sprintf("theta_Breakpoints%drho%.3fpsi3.png",bp,rho);
saveas(gcf,char(savename))
% t
figure
length = size(t,1);
for i = 2:length-1 
    subplot(length-2,1,i-1)
    histfit(t(i,burnIn:end),20)
    ylabel("Numbers")
    str = sprintf('t_%d distribution',i);
    title(str)
end
xlabel("Number of iterations (5000 burn in points)")
savename = sprintf("t_Breakpoints%drho%.3fpsi3.png",bp,rho);
saveas(gcf,char(savename))
% t vary
figure
length = size(t,1);
for i = 2:length-1 
    subplot(length-2,1,i-1)
    plot(t(i,:))
    str = sprintf('t_%d',i);
    ylim([1658 1980])
    yticks(1658:100:1980)
    ylabel("Year")
    title(str)
end
xlabel("Number of iterations (5000 burn in points)")
savename = sprintf("t_Breakpoints%drho%.3fpsi3.png",bp,rho);
saveas(gcf,char(savename))
% t distribution
figure
for i = 2:size(t,1)-1      
    hold on
    histogram(t(i,burnIn:end),40)
    str = sprintf('t distribution with %d breakpoints',bp);
    title(str)
end
xlabel("Year")
ylabel("Number")
savename = sprintf("tdistribution_Breakpoints%drho%.3fpsi3.png",bp,rho);
saveas(gcf,char(savename))
% t in whole
figure
hold on
h = histogram(T,50);
for i=2:size(t,1)-1
    SP = tau(i);
    line([SP SP], [0 max(h.Values)/3], 'Color', [1 0 0],'linewidth',3)
end
xlabel("Year")
ylabel("Numbers")
title('t result')
savename = sprintf("tshow_Breakpoints%drho%.3fpsi3.png",bp,rho);
saveas(gcf,char(savename))
