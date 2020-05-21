%% ������ʼ��
clear all

core_num = 4;%ʹ�ö�CPU���У������ΪCPU������
parpool('local',core_num);

mu_ori = [0.5,0.1,0];% Alice��ǿ�����Ӵ�С˳�����У��������ź�̬
nu_ori = [0.5,0.1,0];% Bob��ǿ�����Ӵ�С˳�����У��������ź�̬

num_mu = length(mu_ori);
num_nu = length(nu_ori);

max_order = 6;% �Ż��б�������������


decoy_ratio = 0.8;% ��ƭ̬ռ��
alpha_0 = 0.2;% ���˺���ϵ��(��λdB/km)

ita_d = 0.1;% ̽����Ч��
Y0 = 0.00001;% ��������
p_bsm = 0.95;% BSM�ɹ���

L = 50000000;% ���ַ������г��ȣ����ù�����ܵ����ڴ治��
loop_num = 400*9;% �ܷ���bit��Ϊloop_num*L

% Attack = 0;% �Ƿ����������

% R_data = [];




%% ����bit����,������ͳ��
l = 72;% ���˳���(��λkm)
t = 10^(-alpha_0*l/10);% ����͸����
mu = mu_ori*t;
nu = nu_ori*t;


Q_z_num_temp_up = zeros(num_mu,num_nu,loop_num);
Q_z_num_temp_down = zeros(num_mu,num_nu,loop_num);
Q_z_num = zeros(num_mu,num_nu,2);
Q_x_num_temp_up = zeros(num_mu,num_nu,loop_num);
Q_x_num_temp_down = zeros(num_mu,num_nu,loop_num);
Q_x_num = zeros(num_mu,num_nu,2);
E_z_num_temp_up = zeros(num_mu,num_nu,loop_num);
E_z_num_temp_down = zeros(num_mu,num_nu,loop_num);
E_z_num = zeros(num_mu,num_nu,2);
E_x_num_temp_up = zeros(num_mu,num_nu,loop_num);
E_x_num_temp_down = zeros(num_mu,num_nu,loop_num);
E_x_num = zeros(num_mu,num_nu,2);



parfor loop = 1:loop_num
    [bit_a,base_a,r_a,s_a,decoy_a] = Decoy_MDI_Photon_Generator(mu(1),mu(2:end),L,'decoy_ratio',decoy_ratio);
    [bit_b,base_b,r_b,s_b,decoy_b] = Decoy_MDI_Photon_Generator(nu(1),nu(2:end),L,'decoy_ratio',decoy_ratio);
    % ================���Ҫ���Ե�����Դ���������ע��ȥ��==================
    % photon_num_a = r_a+s_a;
    % photon_num_b = r_b+s_b;
    % r_a = r_a(bitand((photon_num_a==1),(photon_num_b==1)));
    % s_a = s_a(bitand((photon_num_a==1),(photon_num_b==1)));
    % r_b = r_b(bitand((photon_num_a==1),(photon_num_b==1)));
    % s_b = s_b(bitand((photon_num_a==1),(photon_num_b==1)));
    % base_a = base_a(bitand((photon_num_a==1),(photon_num_b==1)));
    % base_b = base_b(bitand((photon_num_a==1),(photon_num_b==1)));
    % bit_a = bit_a(bitand((photon_num_a==1),(photon_num_b==1)));
    % bit_b = bit_b(bitand((photon_num_a==1),(photon_num_b==1)));
    % ====================================================================
    [result_seq,valid_index] = MDI_Detector(r_a,r_b,s_a,s_b,base_a,base_b,bit_a,bit_b,p_bsm,ita_d,Y0);

    % ͳ�Ƹ��ֽ��
    valid_measure_index = valid_index;
    valid_measure_index(valid_index==1) = result_seq~=-1;% ԭbit�����л���ͬ����Ч�����Ĳ���
    valid_result_seq = result_seq(result_seq~=-1);% ȡ�������������Ч�����Ĳ���
    bit_a_v = bit_a(valid_measure_index);% ȡ��alice����Ч�����һ���ͬ�Ĳ���
    bit_b_v = bit_b(valid_measure_index);% ȡ��bob����Ч�����һ���ͬ�Ĳ���
    base_same = base_a(valid_index);% ȡ��������ͬѡ��Ļ�����
    base = base_a(valid_measure_index);% ȡ������ѡ�õĻ�������ͬ������Ч�����Ĳ���
    % ȡ���������л���ͬ���ֵ���ƭ̬ʹ�����
    decoy_seq_a = decoy_a(valid_index);
    decoy_seq_b = decoy_b(valid_index);
    % ȡ����������Ч�������ֵ���ƭ̬ʹ�����
    decoy_seq_a_v = decoy_a(valid_measure_index);
    decoy_seq_b_v = decoy_b(valid_measure_index);

    % ͳ��Q
    id_z = zeros(num_mu,num_nu,length(decoy_seq_a_v));
    id_x = zeros(num_mu,num_nu,length(decoy_seq_a_v));
    Q_z_temp_up = zeros(num_mu,num_nu);
    Q_z_temp_down = zeros(num_mu,num_nu);
    Q_x_temp_up = zeros(num_mu,num_nu);
    Q_x_temp_down = zeros(num_mu,num_nu);
    for i = 1:num_mu
        for j = 1:num_nu
            % ��Ч����������ʹ��Z������ƭ̬���Ϊij�Ĳ���
            id_z(i,j,:)=bitand(bitand((decoy_seq_a_v==i),(decoy_seq_b_v==j)),(base==0));
            % ��Ч����������ʹ��X������ƭ̬���Ϊij�Ĳ���
            id_x(i,j,:)=bitand(bitand((decoy_seq_a_v==i),(decoy_seq_b_v==j)),(base==1));
            % ʹ��Z������ƭ̬���Ϊij����Ч���������ܺ�
%             Q_z_num(i,j,1) = Q_z_num(i,j,1)+sum(id_z(i,j,:));
%             % ʹ��Z������ƭ̬���Ϊij����bit��
%             Q_z_num(i,j,2) = Q_z_num(i,j,2)+sum(bitand(bitand((decoy_seq_a==i),(decoy_seq_b==j)),(base_same==0)));
%             Q_x_num(i,j,1) = Q_x_num(i,j,1)+sum(id_x(i,j,:));
%             Q_x_num(i,j,2) = Q_x_num(i,j,2)+sum(bitand(bitand((decoy_seq_a==i),(decoy_seq_b==j)),(base_same==1)));
            Q_z_temp_up(i,j) = sum(id_z(i,j,:));
            Q_z_temp_down(i,j) = sum(bitand(bitand((decoy_seq_a==i),(decoy_seq_b==j)),(base_same==0)));
            Q_x_temp_up(i,j) = sum(id_x(i,j,:));
            Q_x_temp_down(i,j) = sum(bitand(bitand((decoy_seq_a==i),(decoy_seq_b==j)),(base_same==1)));
        end
    end
    Q_z_num_temp_up(:,:,loop) = Q_z_temp_up;
    Q_z_num_temp_down(:,:,loop) = Q_z_temp_down;
    Q_x_num_temp_up(:,:,loop) = Q_x_temp_up;
    Q_x_num_temp_down(:,:,loop) = Q_x_temp_down;
    id_z = logical(id_z);
    id_x = logical(id_x);
    
    % ����bit��ת
    % ����X�������ݲ�����������Ƿ�תbit
    bit_a_v(base==1) = mod(bit_a_v(base==1)+valid_result_seq(base==1),2);% x-base
    % ����Z����ֻҪΪ��Ч����������bit��ת
    bit_a_v(base==0) = mod(bit_a_v(base==0)+1,2);% z-base
    
    E_z_temp_up = zeros(num_mu,num_nu);
    E_z_temp_down = zeros(num_mu,num_nu);
    E_x_temp_up = zeros(num_mu,num_nu);
    E_x_temp_down = zeros(num_mu,num_nu);
    % ͳ��E
    for i =1:num_mu
        for j =1:num_nu
%             % ʹ��Z������ƭ̬���Ϊij�Ĳ����д��������
%             E_z_num(i,j,1) = E_z_num(i,j,1)+sum(bit_a_v(id_z(i,j,:))~=bit_b_v(id_z(i,j,:)));
%             % ʹ��Z������ƭ̬���Ϊij�Ĳ��ֵ�������
%             E_z_num(i,j,2) = E_z_num(i,j,2)+sum(id_z(i,j,:));
%             % ʹ��X������ƭ̬���Ϊij�Ĳ����д��������
%             E_x_num(i,j,1) = E_x_num(i,j,1)+sum(bit_a_v(id_x(i,j,:))~=bit_b_v(id_x(i,j,:)));
%             % ʹ��X������ƭ̬���Ϊij�Ĳ��ֵ�������
%             E_x_num(i,j,2) = E_x_num(i,j,2)+sum(id_x(i,j,:));
            E_z_temp_up(i,j) = sum(bit_a_v(id_z(i,j,:))~=bit_b_v(id_z(i,j,:)));
            E_z_temp_down(i,j) = sum(id_z(i,j,:));
            E_x_temp_up(i,j) = sum(bit_a_v(id_x(i,j,:))~=bit_b_v(id_x(i,j,:)));
            E_x_temp_down(i,j) = sum(id_x(i,j,:));
        end
    end
    E_z_num_temp_up(:,:,loop) = E_z_temp_up;
    E_z_num_temp_down(:,:,loop) = E_z_temp_down;
    E_x_num_temp_up(:,:,loop) = E_x_temp_up;
    E_x_num_temp_down(:,:,loop) = E_x_temp_down;
    
end


%% ��������
Q_z_num(:,:,1) = sum(Q_z_num_temp_up,3);
Q_z_num(:,:,2) = sum(Q_z_num_temp_down,3);
Q_x_num(:,:,1) = sum(Q_x_num_temp_up,3);
Q_x_num(:,:,2) = sum(Q_x_num_temp_down,3);
E_z_num(:,:,1) = sum(E_z_num_temp_up,3);
E_z_num(:,:,2) = sum(E_z_num_temp_down,3);
E_x_num(:,:,1) = sum(E_x_num_temp_up,3);
E_x_num(:,:,2) = sum(E_x_num_temp_down,3);


% ����ǿ��ͬ���ֵ��Ӿ���Ϊ�Գ���
[~,ia,ib] = intersect(mu,nu);
Q_x_num(ia,ib,1) = Q_x_num(ia,ib,1)+Q_x_num(ia,ib,1)';
Q_x_num(ia,ib,2) = Q_x_num(ia,ib,2)+Q_x_num(ia,ib,2)';
Q_z_num(ia,ib,1) = Q_z_num(ia,ib,1)+Q_z_num(ia,ib,1)';
Q_z_num(ia,ib,2) = Q_z_num(ia,ib,2)+Q_z_num(ia,ib,2)';
E_x_num(ia,ib,1) = E_x_num(ia,ib,1)+E_x_num(ia,ib,1)';
E_x_num(ia,ib,2) = E_x_num(ia,ib,2)+E_x_num(ia,ib,2)';
E_z_num(ia,ib,1) = E_z_num(ia,ib,1)+E_z_num(ia,ib,1)';
E_z_num(ia,ib,2) = E_z_num(ia,ib,2)+E_z_num(ia,ib,2)';

Q_z = Q_z_num(:,:,1)./Q_z_num(:,:,2);
Q_x = Q_x_num(:,:,1)./Q_x_num(:,:,2);
E_z = E_z_num(:,:,1)./E_z_num(:,:,2);
E_x = E_x_num(:,:,1)./E_x_num(:,:,2);

% �����̬��������ǿ����Ϊ0.5
if mu(end)==0
    E_x(end,:) = 0.5;
    E_z(end,:) = 0.5;
end

if nu(end) == 0
    E_x(:,end)=0.5;
    E_z(:,end)=0.5;
end

%% ����Y11��e11,R
% [Y11_x,e11_x] = evaluate_Y11_and_e11(Q_x,E_x,mu_ori,nu_ori,max_order,100000,1e-8,1e-5,0.01,0.01);
% [Y11_z,e11_z] = evaluate_Y11_and_e11(Q_z,E_z,mu_ori,nu_ori,max_order,100000,1e-8,1e-5,0.01,0.01);
[Y11_x,e11_x] = evaluate_Y11_and_e11_simp(Q_x,E_x,mu_ori,nu_ori,max_order,1e-5,0.01);
[Y11_z,e11_z] = evaluate_Y11_and_e11_simp(Q_z,E_z,mu_ori,nu_ori,max_order,1e-5,0.01);

% ����R
R = max([mu_ori(1)*nu_ori(1)*exp(-mu_ori(1)-nu_ori(1))*Y11_z*(1-Binary_Shannon_Entropy(e11_x))-...
    Q_z(1,1)*Correction_Efficiency(E_z(1,1))*Binary_Shannon_Entropy(E_z(1,1)),1e-10]);
% R_data = [R_data,R];

delete(gcp('nocreate'));
save(['mdisim_',num2str(l),'km.mat']);