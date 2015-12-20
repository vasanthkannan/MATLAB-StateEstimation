% In MATLAB one can produce normally-distributed random variables with an expected value of zero and a standard deviation of 1.0 directly using the function randn. Thus:
% 
%      z = ev + randn(100,10)*sd 
% 
% will produce a {100*10} matrix z of random numbers from a distribution with a mean of ev and a standard deviation of sd.
% 


z = randn(1,1)
