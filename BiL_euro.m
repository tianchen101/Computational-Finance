function [ price_0 ] = BiL_euro( T,N,r,S0,sigma,S_price,call)
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
for i=1:N+1
    
    num_up = i-1;
    num_down = N-(i-1);
    % final_stock_value[i] means goes up (i-1)times
    % and goes down N-(i-1) times
    final_stock_value(i) = S0 * u^(num_up) * d^(num_down);
end
opt_values = zeros(size(final_stock_value,1),1);
strike_prices = ones(N+1,1)*S_price; % construct the strike price vector
if call == 1
    opt_values = max(final_stock_value - strike_prices,0);
else
    opt_values = max(strike_prices-final_stock_value,0);
end

% so we have got the value for all the possible state at the expiary time
% then just to roll back to the initial time, and we can get the price for
% the option now.

for n=N+1:-1:1
    for i=1:n-1
        opt_values(i) = exp(-r*delta_t)*(p*opt_values(i+1)+(1-p)*opt_values(i));
    end
end
% after this the value for option now is opt_value(1)
price_0 = opt_values(1);

end

