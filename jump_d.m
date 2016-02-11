function [ value ] = jump_d( sigma,r,T,K,S0,eta1,eta2,Pup,lambda,N_sim,N)
%
%input parameter
%sigma: votality
%r : interest rate
%T : time to expire
%K: strike price
%S0: initial price
%eta1,eta2: parameter of double exponential
% N_sim: number of simulations
% N: number of steps

% output
% value: the value of the european option

delta_t = T/N; % length of each step
kappa = Pup*eta1/(eta1-1)+(1-Pup)*eta2/(eta2+1)-1;
drift = (r- sigma*sigma/2.0 - lambda*kappa);
X_old(1:N_sim,1) = log(S0);
X_new(1:N_sim,1) = zeros(N_sim,1);

jump_chek = zeros(N_sim,1);
jump_mask = zeros(N_sim,1);
jump_size = zeros(N_sim,1);
jump_type = zeros(N_sim,1);

for i = 1:N
    jump_chek(:,1) = rand(N_sim,1);
    jump_mask(:,1) = (jump_chek(:,1) <= lambda*delta_t);%decide whether or not a jump will happen
    jump_type(:,1) = rand(N_sim,1); % to get the type of the jump,if number <Pup then jump up happend, if number >Pdown jump down
    index_up = logical(jump_type(jump_type <=Pup));% for a randomly generated uniformly distributed number in [0,1],if
    % it < Pup it can be seen as a jump up
    index_down = logical(jump_type(jump_type > Pup));% for that random number, > Pup means it is actually jump down
    
    jump_type(index_up) = 1; % mark the jump ups as 1
    jump_type(index_down) = -1;% mark the jump down as -1
    
    jump_size(index_up) = jump_type(index_up) .* exprnd(1/eta1,size(index_up,1),1);
    jump_size(index_down) = jump_type(index_down) .* exprnd(1/eta2,size(index_down,1),1);
    % note here the mu1 = 1/eta1 and mu2 = 1/eta2
    % note here for the down jump need to multiply by (-1) to get the
    % actual change. Because the number we sapled from the exprnd() is
    % actually give us the inverse of the corresponding change!

    jump_size = jump_size .* jump_mask; % remember only jump_mask != 0 can a jump happen,otherwise the change would be zero
    X_new = X_old(:,1) + drift*delta_t + sigma*sqrt(delta_t)*randn(N_sim,1)+jump_size(:,1); 
    X_old(:,1) = X_new(:,1);
end
S(:,1) = exp(X_new(:,1));
% so have get the possible final stock prices
% so need to trace back to get the price of the european option
S = S-K;% because this is the european call option
S(S <0) =0; % the min value of an option is 0!!
S = S*exp(-r*T);
value = sum(S,1)/N_sim;
