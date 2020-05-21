%% 参数初始化
clear all
mu = 0.6;% 信号态光强
nu1 = 0.2;% 诱骗态光强1
nu2 = 0.1;% 诱骗态光强2
decoy_ratio = 0.7;% 诱骗态占比
alpha_0 = 0.2;% 光纤耗损系数(单位dB/km)
% l = 100;% 光纤长度(单位km)
% t_ab = 10^(-alpha_0*l/10);% 光纤透过率
ita_b = 0.2;% bob端探测效率
% ita = t_ab*ita_b;% 总探测效率（探测到单个光子的概率）
Y0 = 0.00001;% 暗计数率
e_detector = 0.001;% 正常计数误码率
L = 10000000; %每个采样点的仿真序列长度
Attack = 0;% 是否存在窃听者

R_data = [];
%% 进行发送，接收与估算
for l=10:10:140
t_ab = 10^(-alpha_0*l/10);% 光纤透过率
ita = t_ab*ita_b;% 总探测效率（探测到单个光子的概率）
[bit_seq,base_choices,photon_nums,decoy_seq] = Decoy_Photon_Generator(mu,[nu1,nu2],L,'decoy_ratio',decoy_ratio);
if Attack == 1
    % 执行攻击
    photon_nums = BB84_PNS_Attacker(photon_nums);
end
result_seq = BB84_Detector([bit_seq;base_choices;photon_nums],ita,Y0,e_detector);
match_seq = result_seq==bit_seq;
rate = length(find(match_seq==1))/length(bit_seq);
valid_index = result_seq~=-2;
valid_seq = result_seq(valid_index);% 取出双方使用相同基的部分的接收结果
valid_src_seq = bit_seq(valid_index);% 取出双方使用相同基的部分的发送序列
valid_decoy_seq = decoy_seq(valid_index);% 取出双方使用相同基的部分的诱骗态序列
signal_seq = valid_seq(valid_decoy_seq==1);% 取出信号态的接收结果
decoy_seq1 = valid_seq(valid_decoy_seq==2);% 取出使用第一种诱骗态的接收结果
decoy_seq2 = valid_seq(valid_decoy_seq==3);% 取出使用第二种诱骗态的接收结果
signal_src_seq = valid_src_seq(valid_decoy_seq==1);% 取出信号态的发送序列
decoy_src_seq1 = valid_src_seq(valid_decoy_seq==2);% 取出使用第一种诱骗态的发送序列
decoy_src_seq2 = valid_src_seq(valid_decoy_seq==3);% 取出使用第二种诱骗态的发送序列
Q_mu = length(find(signal_seq~=-1))/length(signal_seq);
Q_nu1 = length(find(decoy_seq1~=-1))/length(decoy_seq1);
Q_nu2 = length(find(decoy_seq2~=-1))/length(decoy_seq2);
detected_signal_index = signal_seq~=-1;% 接收方探测到的信号态部分
detected_decoy1_index = decoy_seq1~=-1;% 接收方探测到的第一种诱骗态部分
detected_decoy2_index = decoy_seq2~=-1;% 接收方探测到的第二种诱骗态部分
detected_signal_seq = signal_seq(detected_signal_index);
detected_decoy_seq1 = decoy_seq1(detected_decoy1_index);
detected_decoy_seq2 = decoy_seq2(detected_decoy2_index);
detected_signal_src_seq = signal_src_seq(detected_signal_index);% 发送方相应的信号态部分
detected_decoy_src_seq1 = decoy_src_seq1(detected_decoy1_index);% 发送方相应的第一种诱骗态部分
detected_decoy_src_seq2 = decoy_src_seq2(detected_decoy2_index);% 发送方相应的第二种诱骗态部分
match_signal_index = detected_signal_seq==detected_signal_src_seq;% 双方相同的信号态部分
match_decoy1_index = detected_decoy_seq1==detected_decoy_src_seq1;% 双方相同的第一诱骗态部分
match_decoy2_index = detected_decoy_seq2==detected_decoy_src_seq2;% 双方相同的第二诱骗态部分
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

%% 绘图
semilogy(10:10:140,R_data);
xlabel("Transmission Distance(km)");
ylabel("Key Generation Rate");
set(gca,'Xlim',[10 150]);
