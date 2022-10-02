%======================================================================
%   Numerical Solution of PDEs on Parallel Computers
%   Homework 2, Richardson iteration 
%   Fall 2022, PKU
%======================================================================

% initial layout set up (optional)
clear;clc;close all;
set(0,'defaultlinelinewidth',2)
set(0,'defaultaxeslinewidth',2);
set(0,'defaultaxesfontsize',16);
set(0,'defaulttextfontsize',16);
set(0,'DefaultLineMarkerSize',12);
set(0,'defaultAxesTickLabelInterpreter','none');  

A=[ 6.0 3.0;                            % LHS
    3.0 4.0];
b=[-3.0;-9.0];                          % RHS
lambda=eig(A);                          % get eigenvalues
max_eig=max(lambda);                    % get the max eigenvalue
u= A\b;                                 % the exact solution
alpha_0=[0.06,0.1,0.2,0.22,0.24,0.4];   % relaxation factor
I = eye(2);                             % the 2¡Á2 unit matrix
num_steps = 50;                         % total iterative steps
figure

% Richardson_iteration and plot for every relaxation factor
e_matrix=[];
r_matrix=[];
for i=1:size(alpha_0,2)
    [v,e,r]=Richardson_iteration(A,b,alpha_0(i),num_steps);
    plot(e);
    hold on
    e_matrix=[e_matrix;e];
    r_matrix=[r_matrix;r];
end

set(gca,'yscale','log');
set(gca,'XLim',[0 num_steps]);
set(gca,'YLim',[0 1000]);
set(gca,'FontName','Times New Roman');
xlabel('Steps','interpreter','latex'); ylabel('$|u-u_{ex}|$','interpreter','latex');
legend(strcat('\alpha=',split(num2str(alpha_0))),'Location','southwest');
legend('boxoff');

% output errors & residuals to table
filename='errors&residuals_of_Richardson_iteration.xlsx';
title = {'steps','errors','residuals'};
for i=1:size(alpha_0,2)
    sheet=i;
    result=[[1:50]',e_matrix(1,:)',r_matrix(1,:)'];
    xlswrite(filename,result,sheet);
end
