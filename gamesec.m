clear all;
close all;

psnr = 4;
%T = psnr/2;

alpha = 0.1;
Nnodes = 50;
step = 0.2;
max = 10;
min = 0.1;
%Delta = 6;

%Generate observations
mu = - psnr/2;
var = 1;
T = (mu+(-mu))/2;
%This is to specify the mean and variance under H0



numrep = 100000;
obs = sqrt(var).*randn(numrep,Nnodes) + mu; % generates from N(mu,sigma) not 0,1
%obs = randn(numrep,Nnodes);
% It works 


% Generate corrupted indexes

corr = rand(numrep,Nnodes);
icorr = find(corr < alpha);

i=1;
j=1;

del_len = length(T+min:step:T+max);

for Delta = T+min:step:T+max

%Generate corrupted reports

obs(icorr) = Delta;

%Isolate c

for deltais = min:step:max
    
Tislow = T - deltais;
Tisup = T + deltais;

keptIdx = (obs < Tisup) & (obs > Tislow);
kept = obs .* keptIdx;
numkept = sum(keptIdx,2);
ZZ = find(numkept ~= 0);
means = 2*T*randi([0,1],numrep,1); % this is a problem when T=0, two steps added to solve it
idx0s_mean = find(means==0);
rg = rand(length(idx0s_mean),1);
mean0_tmp = means(idx0s_mean);
mean0_tmp(rg < 0.5) =  T-1;
mean0_tmp(rg >= 0.5) = T+1;

means(idx0s_mean) = mean0_tmp;

means(ZZ) = sum(kept(ZZ,:),2)./numkept(ZZ);



%idx = find(numkept == 0);
%if length(idx) > 0
%    fprintf('Division by zero');
%end

% Count errors

payoff(i,j) = sum(means > T)/numrep;

j = j+1;

end

j = 1;
i = i+1;

fprintf('Strategy %d of %d\n',i-1,del_len);

end

%[STR_A_II,STR_D_II,Pe_II,msg]=zsum(payoff); %Solver 2
%[Pe,STR_A,STR_D] = game_solve(payoff); %Solver 3

%solvers 2,3,4 are exactly the same results with a very tiny diff with
%solver 1 because basically it is not designed for zero sum but Mauro
%edited it

%[STR_A,STR_D,Pe]=zerosum(payoff); %solver 4


eq = LemkeHowson(payoff,-payoff); %solver 1
xax_attack = [T+min:step:T+max];
xax_def = T + [min:step:max];

Perr = eq{1}'*payoff*eq{2};
STR_A = eq{1};
STR_D = eq{2};

save('payoff_matrix','payoff');
save('payoff_game_value.mat','Perr');

figure(1);
bar(xax_attack, STR_A,'red');
xlabel('\Delta');
ylabel('Probability');
set(gca,'fontsize',16);
legend('Attack Mixed Strategy');
grid on;
saveas(gcf,'StrA.jpg');
saveas(gcf,'StrA.fig');


figure(2);
bar(xax_def, STR_D,'blue');
xlabel('\eta');
ylabel('Probability');
set(gca,'fontsize',16);
legend('Defender Mixed Strategy');
grid on;
saveas(gcf,'StrD.jpg');
saveas(gcf,'StrD.fig');


figure(3);
[X,Y] = meshgrid(xax_def, xax_attack);
surf(X,Y,payoff);
xlabel('\eta');
ylabel('\Delta');
zlabel('P_e');
set(gca,'fontsize',16);
view(0,90);
colormap(gray);
colorbar;
grid on;
saveas(gcf,'Game.jpg');
saveas(gcf,'Game.fig');




