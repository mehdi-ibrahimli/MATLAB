% % Numerical Mathematics
% % Matlab Sheet 1
% % RSGI
% % WS20/21
% % Mehdi Ibrahimli

clear all
clc
covid = readtable('coviddata.csv');            % load covid data
% select country and find its indices
country_ind = find(strcmp(covid.COUNTRY_SHORT_NAME, 'Ukraine'));
country_data = covid(country_ind,:);           % extract the data
AZ = sortrows(country_data, {'REPORT_DATE'});  % sort the data by date

%% Gauss-Newton method 
syms a b;                                      % create symbolic a and b for function creation and derivation
alpha = 1000;                                  % alpha is the numerical value of a 
beta = 0.05;                                   % betta is the numerical value of b 
t = [1:height(AZ)]';                           % t is the number of days 
R_s = AZ.PEOPLE_POSITIVE_CASES_COUNT - (a .* exp(b .* t));  % create residual function with symbolic expressions
J_s = jacobian(R_s, [a b]);                    % create Jacobian matrix with symbolic expressions 
% J'(x)* J(x)* (increment) = -J(x)* R(x)         Normal equation                                
for i = 1:20                                   % 20 iteration given in the problem
    J = double(subs(J_s, [a b],[alpha beta])); % convert Jacobian into a numerical matrix with alpha & beta 
    R = double(subs(R_s, [a b],[alpha beta])); % convert residual functions into a numerical vector with alpha & beta
    A = J' * J;                                % to utilize SVD later reqrite normal equation in the form Ax=b
    B = -J'* R;                                % A= J'(x)*J(x) and b = -J(x)*R
    K = cond(A);                               % checking the condition number of the matrix that is to be inverted
     if K < 1000                               % max condition number given in the problem 1000
        xGN = A\B;                             % find x by inversion
     else                                      % else use truncated singular value decomposition
       [U,S,V] = svd(A);
         if S(1,1) - S(2,2) > 0                % min function fails on numbers larger than 2^42
            S(2,2) = 0;                        % if else statement used instead to truncate minimum value 
            S(1,1) = 1/S(1,1);                 % and to simultaneously pseudo invert
         else
            S(1,1) = 0;
            S(2,2) = 1/S(2,2);
         end
       xGN = (V * S' * U') * B;                % finding vector X which is our increment 
    end    
    alpha = alpha + xGN(1);                    % add the increment to our initial parameter values and 
    beta = beta + xGN(2);                      % enter the next iteration
end
% levenberg-marquardt method
% set options for levenberg-marquardt method
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
ydata = AZ.PEOPLE_POSITIVE_CASES_COUNT;        % specify Y data-case numbers (again in case of seperate run)
xdata = [1:height(AZ)]';                       % specify X data-days
fun = @(x,xdata)x(1)*exp(x(2)*xdata);          % speify the model function
x0 = [50,0];                                   % specify the start point (50,0)
lb = [];                                       % lower bound
ub = [];                                       % upper bound
% obtain parameters with the optimization toolbox
xLM = lsqcurvefit(fun,x0,xdata,ydata,lb,ub,options); 

% plotting
hold on
% scatter plot of txhe actual data cases on the y axis and days on x axis
scatter(xdata,AZ.PEOPLE_POSITIVE_CASES_COUNT, 'b');
yGN = alpha * exp(beta*t);                     % obtain predictions of the cases with the fitted model GN method
yLM = xLM(1) * exp(xLM(2)*xdata);              % obtain predictions of the cases with the fitted model LM method
plot(xdata, yLM, 'r','LineWidth',2);           % plot levenberg-marquardt curve 'red'
plot(t, yGN, 'g','LineWidth',2);               % plot gauss-newton curve 'green'
ylabel('total cases');                         % y label
xlabel('Days');                                % x label
title({'Days VS total cases '});               % title
legend('Real cases','levenberg-marquardt','gauss-newton'); 
hold off
