% Numerical Mathematics MATLAB excercise 2
% Mehdi Ibrahimli

function res = page_rank_mult(H,a,alpha,v)
    G = alpha * (H+a) + (1 - alpha) * 1/length(H); % Google matrix
    res = G'* v;                                  
end