function [ price_0  ] = ame(T,N,r,S0,sigma,S_price,call )
% arguments:
%  T : time to expire
%  N : possible number of value at the expiary time
%  r : risk-free interest rate
%  u : the probability that the value of option goes up
%  d : the probability that the value of option goes down
%  S0: current price of underlying asset
%  S_price : the strik price for this option
% call : if call == 1,calculate the price for the call option
%       otherwise calculate the price for the put option
% output:
% price_0: the price of the option at time 0

delta_t = T/N; % delta time
u = exp(sigma * sqrt(delta_t));
d = exp(-sigma * sqrt(delta_t));
final_stock_value = zeros(N+1,1);
p = (exp(r*delta_t) - d)/(u-d);% risk neutral possibility
%and then set all the possible values


W =  S0*d.^([N:-1:0]') .* u.^([0:1:N]');
W1 = W; % W1 is an array to keep the stock price data
if call == 1
% call option
    W = max(W-S_price,0);
else
    W = max(S_price-W,0);
end

for i = N:-1:1
    W = exp(-r*delta_t)*((1-p)*W(1:i)+p*W(2:i+1));
    W1 = W1(1:i)/d;
    if call==1
        W2 = max(W1-S_price,0); % W2 is to get the value for exercising the option
    else
        W2 = max(S_price-W1,0);
    end
    W = max(W,W2);
end
price_0 = W;



end

