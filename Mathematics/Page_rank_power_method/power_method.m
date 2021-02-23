% Numerical Mathematics MATLAB excercise 2
% Mehdi Ibrahimli

function p = power_method(H,a,alpha,v0,tol)                
    G = alpha * (H + 1/length(H) * a) + (1 - alpha) * 1/length(H); % Google matrix
    vk = [];                               % allocate space for vector of a previous iteration
    er = 1;                                % value used to break the loop
    while er == 1                          % the iteration
        vk = v0;                           % assign current Vector to the previous iteration
        v0 = G' * v0;                      % calculate current vector
        v0 = v0/norm(v0,1);                % normalize current vector
        if norm(v0 - vk) < tol             % check the difference between previous iteration
            er = 0;                        % break if necessary
        end
    end  
    p = v0;
end