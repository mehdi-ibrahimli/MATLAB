% % Numerical Mathematics
%% Matlab Sheet 3
% % RSGI
% % WS20/21
% % Mehdi Ibrahimli
%% part a
clear all
close all
clc
n = 20;              % number of samples 20 stated in the problem
rng('shuffle');      % random number seed based on time
x=rand(n,1);         % generate 20 random numbers (0 1) which our samples            
X = repmat(x,1,n);   % repeat the values of the x to raise to the power later
power = [0:n-1];     % create the power array 0 to n
A = X.^power;        % elementwise raising to power gives the polynomial equation without coefficients
f=exp(x);            % vector of e^x
coef = A\f;          % f*A^-1 gives the coefficients
coeff = flip(coef);  % flip because poly2sym function displays polynomial starting from highest degree
%% part b
syms t;                                                 % symbolic variable
hold on
fplot(exp(t),[0 1],'-.*c');                             % plot e^x with symbolic variable on [0 1]
fplot(poly2sym(coeff,t),[0 1],'--or','Linewidth',1);    % plot the approximation (phi^)
hold off


