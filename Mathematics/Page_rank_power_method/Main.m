% Numerical Mathematics MATLAB excercise 2
% Mehdi Ibrahimli

% a). load data
[H,a,URL] = load_data('math_kit.dat');

%% Starting parameters
alpha = 0.8;                                       
v0 = ones(length(H),1)/norm(ones(length(H),1),1);  % starting vector v0
tol = 0.000001;                                    % tolerance as given in the problem
%% 
I = power_method(H,a,alpha,v0,tol);                % identifying the importance vector with power_method func
%% d).
ranked = table(URL,I);                             % create table with URL and their importance values
sorted = sortrows(ranked,2,{'descend'});           % sort the table by their importance values
OUR_Result = sorted(1:10,:);                       % the result of the 10 most important URLs

