clear all;

N = 20;

alfa=0.4;
NA = N*alfa; % the number of attackers on average
sigma_x = 1;
mu_1 = -1;
mu_2 =  1;
Threshold = (mu_1+mu_2)/2;
sigma_x_bar = ((N-NA)/(N^2))*sigma_x;


fun = @(x) exp((-(x-mu_1).^2)/(2*sigma_x_bar^2));

idx = 0;
for i=0:0.01:5% i do not know how to set this boundary for Delta trials but I did so by trial and error.
idx = idx+1;
q(idx) = integral(fun,Threshold - ((NA)/N)*i,Inf);

prob(idx) = (1/sqrt(2*pi*sigma_x_bar^2))*q(idx); 


%prob(idx) = 1 - (1/2)*(1+erf((Threshold - i/N) / sqrt(2*sigma_x_bar)));
%%Correct


Delta(idx) = i;
end

plot(Delta,prob);
grid on;
xlabel('Delta: The attack report value');
ylabel('Attack success probability');
