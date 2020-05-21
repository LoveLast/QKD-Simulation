function [result_seq,valid_index] = MDI_Detector(r_a,r_b,s_a,s_b,base_a,base_b,bit_a,bit_b,p_bsm,ita_d,y0,varargin)
% return result_seq:
% -1 : not valid BSM; 
% 0  : parallel respond; 
% 1  : crossing respond
% valid_index: base_a==base_b
    p = inputParser;
    p.addParameter('HOM_Table',containers.Map(),@ismap);
    p.addParameter('Int_Table',containers.Map(),@ismap);
    p.addParameter('t',0.5,@isnumeric);
    p.addParameter('r',0.5,@isnumeric);
    p.parse(varargin{:});
    valid_index = (base_a==base_b);%get the part with same base selection
    base = base_a(valid_index);
    ra = r_a(valid_index);
    rb = r_b(valid_index);
    sa = s_a(valid_index);
    sb = s_b(valid_index);
    bita = bit_a(valid_index);
    bitb = bit_b(valid_index);
    phase = (bita~=bitb);
    [r_0,r_1,s_0,s_1] = Multi_Photon_BSM(ra,rb,sa,sb,base,phase,p_bsm,varargin{:});
%     r_0_det = zeros(1,length(r_0));
%     r_1_det = zeros(1,length(r_1));
%     s_0_det = zeros(1,length(s_0));
%     s_1_det = zeros(1,length(s_1));
    r_0_det = int8(binornd(double(r_0),ita_d)+binornd(1,y0,1,length(r_0)));
    r_1_det = int8(binornd(double(r_1),ita_d)+binornd(1,y0,1,length(r_1)));
    s_0_det = int8(binornd(double(s_0),ita_d)+binornd(1,y0,1,length(s_0)));
    s_1_det = int8(binornd(double(s_1),ita_d)+binornd(1,y0,1,length(s_1)));
    result_seq = int8(zeros(1,length(r_0)));
    for l = 1:length(r_0)
%         r_0_det(l) = sum(randsrc(1,r_0(l),[1,0;ita_d,1-ita_d]))+randsrc(1,1,[1,0;y0,1-y0]);
%         r_1_det(l) = sum(randsrc(1,r_1(l),[1,0;ita_d,1-ita_d]))+randsrc(1,1,[1,0;y0,1-y0]);
%         s_0_det(l) = sum(randsrc(1,s_0(l),[1,0;ita_d,1-ita_d]))+randsrc(1,1,[1,0;y0,1-y0]);
%         s_1_det(l) = sum(randsrc(1,s_1(l),[1,0;ita_d,1-ita_d]))+randsrc(1,1,[1,0;y0,1-y0]);
        if r_0_det(l)>0&&r_1_det(l)==0
            if s_0_det(l)>0&&s_1_det(l)==0% parallel
                result_seq(l) = 0;
            elseif s_0_det(l)==0&&s_1_det(l)>0% cross
                result_seq(l) = 1;
            else% invalid BSM
                result_seq(l) = -1;
            end
        elseif r_0_det(l)==0&&r_1_det(l)>0
            if s_0_det(l)>0&&s_1_det(l)==0% cross
                result_seq(l) = 1;
            elseif s_0_det(l)==0&&s_1_det(l)>0% parallel
                result_seq(l) = 0;
            else% invalid BSM
                result_seq(l) = -1;
            end
        else % invalid BSM
            result_seq(l) = -1;
        end
    end
end