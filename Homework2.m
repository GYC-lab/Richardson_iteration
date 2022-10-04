%======================================================================
%   Numerical Solution of PDEs on Parallel Computers
%   Homework 2, Richardson iteration 
%   Fall 2022, PKU
%======================================================================

%% initial layout set up (optional)
clear;clc;close all;
set(0,'defaultlinelinewidth',2)
set(0,'defaultaxeslinewidth',2);
set(0,'defaultaxesfontsize',16);
set(0,'defaulttextfontsize',16);
set(0,'DefaultLineMarkerSize',12);
set(0,'defaultAxesTickLabelInterpreter','none');  
format long

%% initial values
A=[ 6.0 3.0;                            % LHS
    3.0 4.0];
b=[-3.0;-9.0];                          % RHS
lambda=eig(A);                          % get eigenvalues
max_eig=max(lambda);                    % get the max eigenvalue
min_eig=min(lambda);                    % get the min eigenvalue
u= A\b;                                 % the exact solution
alpha_0=[0.06,0.1,0.2,0.22,0.24,0.4];   % relaxation factor
I = eye(2);                             % the 2¡Á2 unit matrix
num_steps = 50;                         % total iterative steps

%% Richardson_iteration and plot for every relaxation factor
figure
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

% output figure
img=gcf;
filename='../images/1_convergence_curve_for_every_alpha';
print(img,[filename,'.png'], '-dpng', '-r600');  

%% output errors & residuals to table
filename='errors&residuals.xlsx';
% title = {'steps','errors','residuals'};
result=[[1:num_steps]',e_matrix'];
xlswrite(filename,result,'error');
result=[[1:num_steps]',r_matrix'];
xlswrite(filename,result,'residual');

%% calulate optimal alpha & rho
alpha_opt= 2/(min_eig+max_eig);
rho_opt=(max_eig-min_eig)/(max_eig+min_eig);

%% Richardson_iteration and plot for optimal relaxation factor
figure
[v,e,r]=Richardson_iteration(A,b,alpha_opt,num_steps);
plot(e);
set(gca,'yscale','log');
set(gca,'XLim',[0 num_steps]);
set(gca,'YLim',[0 10]);
set(gca,'FontName','Times New Roman');
xlabel('Steps','interpreter','latex'); ylabel('$|u-u_{ex}|$','interpreter','latex');
legend('$\alpha_{opt}=0.200000$','Location','southwest','interpreter','latex');
legend('boxoff');

% output figure
img=gcf;
filename='../images/2_convergence_curve_for_optimal_alpha';
print(img,[filename,'.png'], '-dpng', '-r600');  

% rate of convergence
convergence_rate=zeros(1,num_steps-1);
for i=1:num_steps-1
    convergence_rate(i)=e(i+1)/e(i);
end

%% Richardson_iteration with left preconditioner
M=[ 6.0 0.0;                            % LHS after preconditioning
    0.0 4.0];
A_pre=inv(M)*A;
b_pre=inv(M)*b;
figure
e_matrix_pre=[];
r_matrix_pre=[];
for i=1:size(alpha_0,2)
    [v,e,r]=Richardson_iteration(A_pre,b_pre,alpha_0(i),num_steps);
    plot(e);
    hold on
    e_matrix_pre=[e_matrix_pre;e];
    r_matrix_pre=[r_matrix_pre;r];
end
set(gca,'yscale','log');
set(gca,'XLim',[0 num_steps]);
set(gca,'YLim',[0 10]);
set(gca,'FontName','Times New Roman');
xlabel('Steps','interpreter','latex'); ylabel('$|u-u_{ex}|$','interpreter','latex');
legend(strcat('\alpha=',split(num2str(alpha_0))),'Location','southwest');
legend('boxoff');

% output figure
img=gcf;
filename='../images/3_convergence_curve_for_every_alpha_after_preconditioning';
print(img,[filename,'.png'], '-dpng', '-r600');  

%% output errors & residuals to table
filename='errors&residuals.xlsx';
result=[[1:num_steps]',e_matrix_pre'];
xlswrite(filename,result,'e_matrix_pre');
result=[[1:num_steps]',r_matrix_pre'];
xlswrite(filename,result,'r_matrix_pre');

%% compare preconditioned and optimal
figure
[v,e,r]=Richardson_iteration(A,b,alpha_opt,num_steps);
plot(e); hold on;
[v,e,r]=Richardson_iteration(A_pre,b_pre,1,num_steps);
plot(e); 
set(gca,'yscale','log');
set(gca,'XLim',[0 num_steps]);
set(gca,'YLim',[0 10]);
set(gca,'FontName','Times New Roman');
xlabel('Steps','interpreter','latex'); ylabel('$|u-u_{ex}|$','interpreter','latex');
legend('$\alpha_{opt}=0.200000$','preconditioned','Location','southwest','interpreter','latex');
legend('boxoff');

% output figure
img=gcf;
filename='../images/4_convergence_curve_for_every_alpha_after_preconditioning';
print(img,[filename,'.png'], '-dpng', '-r600');  

%% calculate the spectral radius of preconditioned iteration
B_=I-1*A_pre;
lambda_=eig(B_);                          % get eigenvalues
rho_=max(abs(lambda_));                    % get the max eigenvaluerho = max(abs(lambda));                 
%% output results on screen
fprintf('lambda_max is %.6f \n',max_eig);
fprintf('lambda_min is %.6f \n',min_eig);
fprintf('range of alpha for convergence is (%.6f,%.6f) \n',0,2/max_eig);
fprintf('alpha_opt is %.6f \n',alpha_opt);
fprintf('rho_opt is %.6f \n',rho_opt);
fprintf('rate of convergence is %.6f \n',convergence_rate(end));
fprintf('spectral radius of (I-aA) is %.6f \n',rho_);
