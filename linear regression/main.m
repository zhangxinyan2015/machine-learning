clc
clear all
close all
x=load('ex2x.dat');
y=load('ex2y.dat');
plot(x,y,'o')
ylabel('Height in meters')
xlabel('Age in years')
m = length(y); % store the number of training examples
x = [ones(m, 1), x]; % Add a column of ones to x
%linear regression
alpha=0.07;
theta=[0 0]';
max_iter=1500;
for i=1:max_iter
    grad=1/m*x'*(x*theta-y);    %梯度下降
    theta=theta-alpha*grad;
end
hold on % Plot new data without clearing old plot
plot(x(:,2), x*theta, '-') % remember that x is now a matrix with 2 columns
                           % and the second column contains the time info
legend('Training data', 'Linear regression')
% understanding J(theta)
J_vals = zeros(100, 100);   % initialize Jvals to 100x100 matrix of 0's
theta0_vals = linspace(-3, 3, 100);
theta1_vals = linspace(-1, 1, 100);
for i = 1:length(theta0_vals)
	  for j = 1:length(theta1_vals)
	  t = [theta0_vals(i); theta1_vals(j)];
	  %J_vals(i,j) = 1/(2*m)*sum((x*t-y).^2);
      J_vals(i,j) = 1/(2*m)*(x*t-y)'*(x*t-y);%更好的利用矩阵的方法：vertorized implementation
    end
end
J_vals = J_vals';
figure;
surf(theta0_vals, theta1_vals, J_vals)%surf函数横是x,纵是y
xlabel('\theta_0'); ylabel('\theta_1')
figure;
contour(theta0_vals, theta1_vals, J_vals, logspace(-2, 2, 15))
xlabel('\theta_0'); ylabel('\theta_1')