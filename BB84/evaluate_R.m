clear all

% =======================����˫��ʹ�õĹ�ǿ===============================
mu = 0.6;% �ź�̬��ǿ
nu1 = 0.2;% ��ƭ̬��ǿ1
nu2 = 0.1;% ��ƭ̬��ǿ2

% =======================����ʵ���������=================================
% ����
Q_mu = 0;
Q_nu1 = 0;
Q_nu2 = 0;

% ������
E_mu = 0;
E_nu1 = 0;
E_nu2 = 0;

% ����Y_11��e_11
Y0_Est = max([(nu1*Q_nu2*exp(nu2)-nu2*Q_nu1*exp(nu1))/(nu1-nu2),0]);
Y1_Est = max([mu/(mu*(nu1-nu2)-(nu1^2-nu2^2))*(Q_nu1*exp(nu1)-Q_nu2*exp(nu2)-(nu1^2-nu2^2)/mu^2*(Q_mu*exp(mu)-Y0_Est)),0]);
e1_Est = max([min([(E_nu1*Q_nu1*exp(nu1)-E_nu2*Q_nu2*exp(nu2))/(nu1-nu2)/Y1_Est,1]),0]);

% ����R
R_Est = max([mu*exp(-mu)*Y1_Est*(1-Binary_Shannon_Entropy(e1_Est))-Q_mu*Correction_Efficiency(E_mu)*Binary_Shannon_Entropy(E_mu),0]);
R_Est = 0.5*R_Est; % q = 0.5