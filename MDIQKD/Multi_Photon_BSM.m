function [r0,r1,s0,s1] = Multi_Photon_BSM(r_a,r_b,s_a,s_b,base,phase,p_bsm,varargin)
    p = inputParser;
    p.addParameter('HOM_Table',containers.Map(),@ismap);
    p.addParameter('Int_Table',containers.Map(),@ismap);
    p.addParameter('t',0.5,@isnumeric);
    p.addParameter('r',0.5,@isnumeric);
    p.parse(varargin{:});
    HOM_Table = p.Results.HOM_Table;
    Int_Table = p.Results.Int_Table;
    t = p.Results.t;% transmissivity
    r = p.Results.r;% reflectivity
    L = length(r_a);
    r0 = int8(zeros(1,L));
    r1 = int8(zeros(1,L));
    s0 = int8(zeros(1,L));
    s1 = int8(zeros(1,L));
    
    M = r_a+s_a;
    N = r_b+s_b;
    R = r_a+r_b;
    S = s_a+s_b;
    
    zeros_input_index = bitand(M==0,N==0);
    r0(zeros_input_index)=0;
    r1(zeros_input_index)=0;
    s0(zeros_input_index)=0;
    s1(zeros_input_index)=0;
    
    successful_BSM = binornd(1,p_bsm,1,L);
    r0(successful_BSM==0) = binornd(double(r_a(successful_BSM==0)),r)+binornd(double(r_b(successful_BSM==0)),t);
    s0(successful_BSM==0) = binornd(double(s_a(successful_BSM==0)),r)+binornd(double(s_b(successful_BSM==0)),t);
    r1(successful_BSM==0) = R(successful_BSM==0)-r0(successful_BSM==0);
    s1(successful_BSM==0) = S(successful_BSM==0)-s0(successful_BSM==0);

    
    comma = ',';
    MultiPhoton = [getchar(M(:)),comma(ones(L,1),:),getchar(N(:)),comma(ones(L,1),:),getchar(phase(:))];
    MultiPhotonR = [getchar(r_a(:)),comma(ones(L,1),:),getchar(r_b(:))];
    MultiPhotonS = [getchar(s_a(:)),comma(ones(L,1),:),getchar(s_b(:))];
    comb_zR = containers.Map();
    comb_zS = containers.Map();
    comb_x = containers.Map();
    
    linear_sp = 1:L;
    res_bit_index_x = bitand(bitand(successful_BSM==1,base==1),~zeros_input_index);
    res_bit_index_zR = bitand(bitand(successful_BSM==1,base==0),~zeros_input_index);
    res_bit_index_zS = bitand(bitand(successful_BSM==1,base==0),~zeros_input_index);
    
    i = find(res_bit_index_x,1);
    while ~isempty(i)        
        comb_x_index = bitand(bitand(bitand(bitand(M==M(i),N==N(i)),phase==phase(i)),base==1),successful_BSM==1);
        comb_x(MultiPhoton(i,:)) = linear_sp(comb_x_index);
        res_bit_index_x = bitand(res_bit_index_x,~comb_x_index);
                
        if ~isKey(Int_Table,MultiPhoton(i,:))
            Int_Table(MultiPhoton(i,:)) = Multi_Photon_Interference(double(M(i)),double(N(i)),phase(i));
        end
        i = find(res_bit_index_x,1);
    end
    
    i = find(res_bit_index_zR,1);
    while ~isempty(i)% z-base-R mode
        comb_zR_index = bitand(bitand(bitand(bitand(r_a==r_a(i),r_b==r_b(i)),base==0),...
            successful_BSM==1),~zeros_input_index);
        comb_zR(MultiPhotonR(i,:)) = linear_sp(comb_zR_index);
        res_bit_index_zR = bitand(res_bit_index_zR,~comb_zR_index);
        if ~isKey(HOM_Table,MultiPhotonR(i,:))
            HOM_Table(MultiPhotonR(i,:))=Multi_Photon_HOM(double(r_a(i)),double(r_b(i)),t,r);
        end
        i = find(res_bit_index_zR,1);
    end
    
    i = find(res_bit_index_zS,1);
    while ~isempty(i)
        comb_zS_index = bitand(bitand(bitand(bitand(s_a==s_a(i),s_b==s_b(i)),base==0),...
            successful_BSM==1),~zeros_input_index);
        comb_zS(MultiPhotonS(i,:)) = linear_sp(comb_zS_index);
        res_bit_index_zS = bitand(res_bit_index_zS,~comb_zS_index);
        if ~isKey(HOM_Table,MultiPhotonS(i,:))
            HOM_Table(MultiPhotonS(i,:))=Multi_Photon_HOM(double(s_a(i)),double(s_b(i)),t,r);
        end
        i = find(res_bit_index_zS,1);
    end
    
    X_Key = keys(comb_x);
    ZR_Key = keys(comb_zR);
    ZS_Key = keys(comb_zS);
    for i = 1:length(X_Key)
        State = Int_Table(X_Key{i});
        index = randsrc(1,length(comb_x(X_Key{i})),State(1:2,:));
        r0(comb_x(X_Key{i})) = State(3,index);
        r1(comb_x(X_Key{i})) = State(4,index);
        s0(comb_x(X_Key{i})) = State(5,index);
        s1(comb_x(X_Key{i})) = State(6,index);
    end
    for i = 1:length(ZR_Key)
        R_State = HOM_Table(ZR_Key{i});
        r0(comb_zR(ZR_Key{i})) = randsrc(1,length(comb_zR(ZR_Key{i})),R_State);
        r1(comb_zR(ZR_Key{i})) = R(comb_zR(ZR_Key{i}))-r0(comb_zR(ZR_Key{i}));
    end
    for i = 1:length(ZS_Key)
        S_State = HOM_Table(ZS_Key{i});
        s0(comb_zS(ZS_Key{i})) = randsrc(1,length(comb_zS(ZS_Key{i})),S_State);
        s1(comb_zS(ZS_Key{i})) = S(comb_zS(ZS_Key{i}))-s0(comb_zS(ZS_Key{i}));
    end
end

function res = ismap(obj)
    res = isa(obj,'containers.Map');
end