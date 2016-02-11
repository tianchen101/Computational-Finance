function [ price_0 ] = BiL_America(T,N,r,S0,sigma,S_price,call)
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
p = (exp(r*delta_t) - d)/(u-d);

for i=1:N+1
    
    num_up = i-1;
    num_down = N-(i-1);
    % final_stock_value[i] means goes up (i-1)times
    % and goes down N-(i-1) times
    final_stock_value(i) = S0 * u^(num_up) * d^(num_down);
end
final_stock_value_copy = final_stock_value; % a copy of the final stock prices
second_last_stock_value = zeros(N,1); % the stock price at t = (N-1) * delta_t
for i=1:N
    second_last_stock_value(i) = final_stock_value(i)*u; 
end

% flag = 0; % will be used to calculate the initial stock price
% 
% if rem(N+1,2)==0 % means that (N+1) is a even number
%     flag = 1;
% else %otherwise N is even, and (N+1) is odd
%     flag = 0;
% end

% will be used to calculate the option price later
    
opt_values = zeros(size(final_stock_value,1),1);
strike_prices = ones(N+1,1)*S_price; % construct the strike price vector
if call == 1 % call option
    opt_values = max(final_stock_value - strike_prices,0);
    opt_values_second_last = max(second_last_stock_value - strike_prices(1:N),0);
else % put option
    opt_values = max(strike_prices-final_stock_value,0);
    opt_values_second_last =  max(strike_prices(1:N)-second_last_stock_value,0);
end

for n=N+1:-1:1
    if rem(N+1 - n,2) == 0
        flag = 1; % use the array derived from the array second_last_stock_price
    else
        flag = 0; % use the array derived from the array final_stock_value_copy
    end
    % next step is to generate the array at the certain time
    if flag==1
        k = (N+1-n)/2;
        stock_price_now = second_last_stock_value(1+k:N-k);
    else% flag ==0
        k = (N-n)/2+1;
        stock_price_now = final_stock_value(1+k:N+1-k);
    end
    if call==1
        opt_value_now = max( stock_price_now - ones( size( stock_price_now,1),1)*S_price,0);
    else
        opt_value_now = max( ones( size( stock_price_now,1 ),1)*S_price-stock_price_now,0);
    end
    for i=1:n-1
        temp_value1 = exp(-r*delta_t)*(p*opt_values(i+1)+(1-p)*opt_values(i));
        opt_values(i) = max(temp_value1,opt_value_now(i));
    end
end
% after this the value for option now is opt_value(1)
price_0 = opt_values(1);

end

