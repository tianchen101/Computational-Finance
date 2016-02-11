function [ value  ] = drift_S( sigma,r,T,S0,K,N,call )
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
    %W = max(W-K,0);
   index1 = W*exp(sigma*sqrt(delta_t)) <K;
   index2 = W*exp(-sigma*sqrt(delta_t)) > K;
   index3 = logical((1-index2).*(1-index1)); % this is the interval for the final stock value in between
   W(index1)=0;
   W(index2) =W(index2)* ((exp(sigma*sqrt(delta_t))-exp(-sigma*sqrt(delta_t)))/(2*sigma*sqrt(delta_t)))-K;
   if sum(index3)>=1
   W(index3) = (1/sqrt(2*sigma*sqrt(delta_t)))*(W(index3).*(exp(sigma*sqrt(delta_t))-(K/W(index3)))-K*(sigma*sqrt(delta_t)-log(K/W(index3))));
   end
else
    % means this is the put option
    index1 = W*exp(-sigma*sqrt(-sigma*sqrt(delta_t))) >K;
    index2 = W*exp(sigma*sqrt(delta_t)) <K;
    index3 = logical((1-index2).*(1-index1));
    W(index1)=0;
    W(index2) = K - W(index2)*((exp(sigma*sqrt(delta_t))-exp(-sigma*sqrt(delta_t)))/(2*sigma*sqrt(delta_t)));
    if sum(index3)>=1
    W(index3) = (1/(2*sigma*sqrt(delta_t)))*(K*(log(K/W(index3))+sigma*sqrt(delta_t))-W(index3).*((K/W(index3))-exp(-sigma*sqrt(delta_t))));
    end
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

end

