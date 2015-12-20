% test chi2inv


clear



K=50;
max = 200;
alpha = zeros(max,1);
P = zeros(max,1);



max_alpha = 1.0;
alpha(1) = max_alpha;
for i = 1:max
    
    TJ(i) = 1 - chi2inv(alpha(i),K);
%     TJ(i) = 2*K - chi2inv(alpha(i),K);
    if i < max
       alpha(i+1) = alpha(i) - max_alpha/max;
    end
end
% 
% max_alpha = 1.0;
% alpha(1) = max_alpha;
% for i = 1:max
%     
%     TJ(i) = chi2inv(alpha(i),K);
%     if i < max
%        alpha(i+1) = alpha(i) - max_alpha/max;
%     end
% end
plot(TJ,alpha)


% working with cdf function

% TJ(1) = 0.0;
% for i = 1:max
%     
%     alpha(i) = ( 1 - chi2cdf(TJ(i),K) );
%     if i < max
%          TJ(i+1) = TJ(i) + 0.5;
%     end
% end
% 
% TJ(1) = 0.0;
% for i = 1:max
%     P(i) = chi2cdf(TJ(i),K);
%     if i < max
%          TJ(i+1) = TJ(i) + 0.5;
%     end
% end
% plot(TJ, P)

% for i = 1:max
%     fprintf(' %7.4f   %5.4f \n',TJ(i), alpha(i) )
% end

% plot(TJ, alpha)