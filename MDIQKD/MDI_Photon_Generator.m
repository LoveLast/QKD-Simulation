function [bit_seq,base_choices,r_nums,s_nums] = MDI_Photon_Generator(mu,L,varargin)
    p = inputParser;
    p.addParameter('input_bit',int8([]),@isnumeric);
    p.addParameter('base',int8([]),@isnumeric);
    p.addParameter('decoy_bit',int8([]),@isnumeric);
    p.addParameter('r_nums',int8([]),@isnumeric);
    p.addParameter('s_nums',int8([]),@isnumeric);
    p.addParameter('photon_num',int8([]),@isnumeric);
    p.addParameter('decoy_ratio',0.1,@isnumeric);
    p.parse(varargin{:});
    bit_seq = int8(p.Results.input_bit);
    if length(p.Results.input_bit)<L
        %make the bit seq the same length
        bit_seq = [bit_seq,randi(2,1,L-length(p.Results.input_bit))-1];
    end
    base_choices = int8(p.Results.base);
    if length(p.Results.base)<L
        %make the base choice seq the same length
        base_choices = [base_choices,randi(2,1,L-length(p.Results.base))-1];
    end
    photon_nums = int8(p.Results.photon_num);
    if length(p.Results.photon_num)<L
        %make the numbers of photon seq the same length
        photon_nums = [photon_nums,random('Poisson',mu,1,L-length(p.Results.photon_num))];
    end
    r_nums = int8(zeros(1,L));
    
    r_nums(base_choices==1)=binornd(double(photon_nums(base_choices==1)),0.5);%x-base
    r_nums(bitand(base_choices==0,bit_seq==0))=photon_nums(bitand(base_choices==0,bit_seq==0));%z-base
    r_nums(1:length(p.Results.r_nums)) = p.Results.r_nums;
    s_nums = photon_nums-r_nums;
    
    
end