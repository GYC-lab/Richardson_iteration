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

num_steps = 10;
figure
for i=1:size(alpha_0,2)
    [v,e,r]=Richardson_iteration(A,b,alpha_0(i),num_steps);
%     plot(e);
%     hold on
end
