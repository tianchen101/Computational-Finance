function [ value] = drift_L( sigma,r,T,S0,K,N,call )
%Question 2: Implementing the drift lattice for American Put Option
% sigma: the votality
% r: the interest rate
% T: time to expirory
% S0: the initial stock price
% k: the strike price
% N: number of steps
% call: if call option then call==1,else call=0
% output: the predicted price for the option right now
delta_t = T/N;
u = exp(sigma*sqrt(delta_t)+(r-sigma^2/2)*delta_t);
d = exp(-sigma*sqrt(delta_t)+(r-sigma^2/2)*delta_t);
W = S0 * d.^([N:-1:0]') .* u.^([0:1:N]'); % final stock value
W1 = W; % W1 is a copy of the final stock price 
if call==1
    % means this is the call option
    W = max(W-K,0);
else
    % means this is the put option
    W = max(K-W,0);
end

for i = N:-1:1
    W = exp(-r*delta_t)*((1/2)*W(2:i+1) + (1/2)*W(1:i));% the expected value of option at previous step
    W1 = W1(1:i)/d; % the stock prices at previous step
    % choose the optimal to exercise
    if call==1
        W2 = max(W1-K,0);
    else
        W2 = max(0,K-W1); % value for put option
    end
    W = max(W,W2); % value for call option
end
% after the for loop the final value would be stored in W
value = W;

