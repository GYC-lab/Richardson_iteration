%======================================================================
%   Numerical Solution of PDEs on Parallel Computers
%   Homework 2, Richardson iteration 
%   Fall 2022, PKU
%======================================================================

% initial layout set up (optional)
clear;clc;close all;

A=[ 6.0 3.0;
    3.0 4.0];
b=[-3.0;-9.0];
lambda=eig(A);                          % get eigenvalues
max_eig=max(lambda);                    % get the max eigenvalue
u= A\b;                                 % the exact solution

alpha_0=[0.06,0.1,0.2,0.22,0.24,0.4];   % relaxation factor
I = eye(2);                             % the 2¡Á2 unit matrix
for i=1:1:size(alpha_0,2)
    x1=[0;0];
    alpha=alpha_0(i);
    B=I-alpha*A;
    C=alpha*b;
    for k=1:500
        x2=B*x1+C;
        x1=x2;
    end
end
