%% ������ʼ��
clear all
mu = 0.6;% �ź�̬��ǿ
nu1 = 0.2;% ��ƭ̬��ǿ1
nu2 = 0.1;% ��ƭ̬��ǿ2
decoy_ratio = 0.7;% ��ƭ̬ռ��
alpha_0 = 0.2;% ���˺���ϵ��(��λdB/km)
% l = 100;% ���˳���(��λkm)
% t_ab = 10^(-alpha_0*l/10);% ����͸����
ita_b = 0.2;% bob��̽��Ч��
% ita = t_ab*ita_b;% ��̽��Ч�ʣ�̽�⵽�������ӵĸ��ʣ�
Y0 = 0.00001;% ��������
e_detector = 0.001;% ��������������
L = 10000000; %ÿ��������ķ������г���
Attack = 0;% �Ƿ����������

R_data = [];
%% ���з��ͣ����������
for l=10:10:140
t_ab = 10^(-alpha_0*l/10);% ����͸����
ita = t_ab*ita_b;% ��̽��Ч�ʣ�̽�⵽�������ӵĸ��ʣ�
[bit_seq,base_choices,photon_nums,decoy_seq] = Decoy_Photon_Generator(mu,[nu1,nu2],L,'decoy_ratio',decoy_ratio);
if Attack == 1
    % ִ�й���
    photon_nums = BB84_PNS_Attacker(photon_nums);
end
result_seq = BB84_Detector([bit_seq;base_choices;photon_nums],ita,Y0,e_detector);
match_seq = result_seq==bit_seq;
rate = length(find(match_seq==1))/length(bit_seq);
valid_index = result_seq~=-2;
valid_seq = result_seq(valid_index);% ȡ��˫��ʹ����ͬ���Ĳ��ֵĽ��ս��
valid_src_seq = bit_seq(valid_index);% ȡ��˫��ʹ����ͬ���Ĳ��ֵķ�������
valid_decoy_seq = decoy_seq(valid_index);% ȡ��˫��ʹ����ͬ���Ĳ��ֵ���ƭ̬����
signal_seq = valid_seq(valid_decoy_seq==1);% ȡ���ź�̬�Ľ��ս��
decoy_seq1 = valid_seq(valid_decoy_seq==2);% ȡ��ʹ�õ�һ����ƭ̬�Ľ��ս��
decoy_seq2 = valid_seq(valid_decoy_seq==3);% ȡ��ʹ�õڶ�����ƭ̬�Ľ��ս��
signal_src_seq = valid_src_seq(valid_decoy_seq==1);% ȡ���ź�̬�ķ�������
decoy_src_seq1 = valid_src_seq(valid_decoy_seq==2);% ȡ��ʹ�õ�һ����ƭ̬�ķ�������
decoy_src_seq2 = valid_src_seq(valid_decoy_seq==3);% ȡ��ʹ�õڶ�����ƭ̬�ķ�������
Q_mu = length(find(signal_seq~=-1))/length(signal_seq);
Q_nu1 = length(find(decoy_seq1~=-1))/length(decoy_seq1);
Q_nu2 = length(find(decoy_seq2~=-1))/length(decoy_seq2);
detected_signal_index = signal_seq~=-1;% ���շ�̽�⵽���ź�̬����
detected_decoy1_index = decoy_seq1~=-1;% ���շ�̽�⵽�ĵ�һ����ƭ̬����
detected_decoy2_index = decoy_seq2~=-1;% ���շ�̽�⵽�ĵڶ�����ƭ̬����
detected_signal_seq = signal_seq(detected_signal_index);
detected_decoy_seq1 = decoy_seq1(detected_decoy1_index);
detected_decoy_seq2 = decoy_seq2(detected_decoy2_index);
detected_signal_src_seq = signal_src_seq(detected_signal_index);% ���ͷ���Ӧ���ź�̬����
detected_decoy_src_seq1 = decoy_src_seq1(detected_decoy1_index);% ���ͷ���Ӧ�ĵ�һ����ƭ̬����
detected_decoy_src_seq2 = decoy_src_seq2(detected_decoy2_index);% ���ͷ���Ӧ�ĵڶ�����ƭ̬����
match_signal_index = detected_signal_seq==detected_signal_src_seq;% ˫����ͬ���ź�̬����
match_decoy1_index = detected_decoy_seq1==detected_decoy_src_seq1;% ˫����ͬ�ĵ�һ��ƭ̬����
match_decoy2_index = detected_decoy_seq2==detected_decoy_src_seq2;% ˫����ͬ�ĵڶ���ƭ̬����
E_mu = length(find(match_signal_index==0))/length(detected_signal_src_seq);
E_nu1 = length(find(match_decoy1_index==0))/length(detected_decoy_src_seq1);
E_nu2 = length(find(match_decoy2_index==0))/length(detected_decoy_src_seq2);
Y0_Est = max([(nu1*Q_nu2*exp(nu2)-nu2*Q_nu1*exp(nu1))/(nu1-nu2),0]);
Y1_Est = max([mu/(mu*(nu1-nu2)-(nu1^2-nu2^2))*(Q_nu1*exp(nu1)-Q_nu2*exp(nu2)-(nu1^2-nu2^2)/mu^2*(Q_mu*exp(mu)-Y0_Est)),0]);
e1_Est = max([min([(E_nu1*Q_nu1*exp(nu1)-E_nu2*Q_nu2*exp(nu2))/(nu1-nu2)/Y1_Est,1]),0]);

R_Est = max([mu*exp(-mu)*Y1_Est*(1-Binary_Shannon_Entropy(e1_Est))-Q_mu*Correction_Efficiency(E_mu)*Binary_Shannon_Entropy(E_mu),0]);
R_Est = 0.5*R_Est; % q = 0.5
R_data = [R_data,R_Est];
end

%% ��ͼ
semilogy(10:10:140,R_data);
xlabel("Transmission Distance(km)");
ylabel("Key Generation Rate");
set(gca,'Xlim',[10 150]);
