clear all

max_order = 6;% 估计时保留的最大光子数

% =======================输入双方使用的光强===============================
mu = [0.5,0.2,0.1,0];% Alice光强，按从大到小顺序排列，最大的是信号态
nu = [0.5,0.2,0.1,0];% Bob光强，按从大到小顺序排列，最大的是信号态


% =======================输入实验测试数据=================================

% 注意：如果双方使用的光强相同增益与误码率应该是一个对称阵
% Z基的增益Q_z（符合计数/基匹配下总bit数），其中第i行第j列代表使用mu_i，nu_j时的增益
Q_z = [5.06e-5,2.02e-5,1.01e-5,2.51e-8;...
       1.99e-5,7.97e-6,4.02e-6,8.67e-9;...
       9.87e-6,4.01e-6,1.98e-6,4.49e-9;
       6.22e-8,1.38e-8,4.19e-9,0];
   
% X基的增益，同上
Q_x = [1.06e-4,5.18e-5,3.85e-5,2.73e-5;...
       5.20e-5,1.67e-5,9.51e-6,4.38e-6;...
       3.81e-5,9.30e-6,4.20e-6,1.12e-6;
       2.68e-5,4.16e-6,1.00e-6,3.0e-10];
   
   

% Z基的误码率E_z（错误bit数/符合计数），其中第i行第j列代表使用mu_i，nu_j时的误码率
E_z = [0.14,0.17,0.22,50;...
       0.22,0.19,0.21,41.38;...
       0.38,0.24,0.33,40;...
       53.37,50,35.71,50]*10^-2;
   
% X基误码率，同上
E_x = [27.70,32.13,37.84,50.13;...
       31.90,27.67,31.02,51.07;...
       37.93,30.41,27.08,50.63;...
       49.85,49.74,50.36,50]*10^-2;
   

% 注意：如果使用任意一方使用真空态，其对应的误码率应是0.5
% if mu(end)==0
%     E_x(end,:) = 0.5;
%     E_z(end,:) = 0.5;
% end
% 
% if nu(end) == 0
%     E_x(:,end)=0.5;
%     E_z(:,end)=0.5;
% end

% 估计Y11和e11
[Y11_x,e11_x] = evaluate_Y11_and_e11_simp(Q_x,E_x,mu,nu,max_order,1e-5,0.01);
[Y11_z,e11_z] = evaluate_Y11_and_e11_simp(Q_z,E_z,mu,nu,max_order,1e-5,0.01);

% 计算R，注意如果R过小则视为1e-12，表示没有有效的安全码率
R = max([mu(1)*nu(1)*exp(-mu(1)-nu(1))*Y11_z*(1-Binary_Shannon_Entropy(e11_x))-...
    Q_z(1,1)*Correction_Efficiency(E_z(1,1))*Binary_Shannon_Entropy(E_z(1,1)),1e-12]);

   