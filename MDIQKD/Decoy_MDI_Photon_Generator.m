function [bit_seq,base_choices,r_nums,s_nums,decoy_seq] = Decoy_MDI_Photon_Generator(mu,nu,L,varargin)
% note: r,s nums are only valid for z-base, not valid for x-base
    p = inputParser;
    p.addParameter('input_bit',int8([]),@isnumeric);
    p.addParameter('base',int8([]),@isnumeric);
    p.addParameter('decoy_bit',int8([]),@isnumeric);
    p.addParameter('photon_num',int8([]),@isnumeric);
    p.addParameter('r_nums',int8([]),@isnumeric);
    p.addParameter('s_nums',int8([]),@isnumeric);
    p.addParameter('decoy_ratio',0.1,@isnumeric);
    p.parse(varargin{:});
    decoy_state_num = length(nu);
    decoy_nu = nu;
    if decoy_state_num<2
        decoy_nu = [nu,0];
        decoy_state_num = decoy_state_num+1;
    end
    decoy_nu = [mu,decoy_nu];
    decoy_seq = int8(p.Results.decoy_bit);
    if length(p.Results.decoy_bit)<L-length(p.Results.photon_num)
        decoy_seq = [decoy_seq,randsrc(1,L-length(p.Results.decoy_bit)-length(p.Results.photon_num),...
        [1:decoy_state_num+1;1-p.Results.decoy_ratio,p.Results.decoy_ratio/decoy_state_num*ones(1,decoy_state_num)])];
    end
    decoy_mu_seq = decoy_nu(decoy_seq);
    [bit_seq,base_choices,r_nums,s_nums] = MDI_Photon_Generator(decoy_mu_seq,L,varargin{:});
end