function [v,e_norm,r_norm] = Richardson_iteration(A,b,alpha,num_steps)
%RICHARDSON_ITERATION - solve a linear system like Ax = b
%   input : A           - matrix
%           b           - vector
%           alpha       - relaxation factor 
%           num_steps   - number of iterative steps
%   output: v    - the numerical solution
%           r    - residual
%           e    - error

    % initialization
    u = A\b;                % the exact solution
    I = eye(2);             % the 2¡Á2 unit matrix
    x1=[0.0;0.0];
    
    % Richardson_iteration
    B=I-alpha*A;
    C=alpha*b;
    e_norm=zeros(1,num_steps);
    r_norm=zeros(1,num_steps);
    for i=1:num_steps
        % iteration
        x2=B*x1+C;
        x1=x2;
        % error and residual
        v=x2;
        e=v-u;      % errors
        r=b-A*v;    % residuals    
        e_norm(i)=norm(e); % 2-norm of errors
        r_norm(i)=norm(r); % 2-norm of residuals
    end
    
end

