function [bit_seq,base_choices,photon_nums] = Photon_Generator(mu,L,varargin)
    p = inputParser;
    p.addParameter('input_bit',[],@isnumeric);
    p.addParameter('base',[],@isnumeric);
    p.addParameter('decoy_bit',[],@isnumeric);
    p.addParameter('photon_num',[],@isnumeric);
    p.addParameter('decoy_ratio',0.1,@isnumeric);
    p.parse(varargin{:});
    bit_seq = p.Results.input_bit;
    if length(p.Results.input_bit)<L
        %make the bit seq the same length
        bit_seq = [p.Results.input_bit,randi(2,1,L-length(p.Results.input_bit))-1];
    end
    base_choices = p.Results.base;
    if length(p.Results.base)<L
        %make the base choice seq the same length
        base_choices = [p.Results.base,randi(2,1,L-length(p.Results.base))-1];
    end
    photon_nums = p.Results.photon_num;
    if length(p.Results.photon_num)<L
        %make the numbers of photon seq the same length
        photon_nums = [p.Results.photon_num,random('Poisson',mu,1,L-length(p.Results.photon_num))];
    end
end