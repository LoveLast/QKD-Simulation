function result_seq = BB84_Detector(photon_seq,ita,y0,e_detector,varargin)
% photon_seq = [bit_seq;base_seq;number_seq]
% ita = p(detected)
% e_detector = p(respond error|detected)
% result_seq: -1: no detected;-2: wrong base
    p = inputParser;
    p.addParameter('measure_base',[],@isnumeric);
    p.parse(varargin{:});
    L = size(photon_seq,2);
    bit_seq = photon_seq(1,:);
    base_seq = photon_seq(2,:);
    number_seq = photon_seq(3,:);
    result_seq = zeros(1,L);
    measure_base_seq = p.Results.measure_base;
    if length(p.Results.measure_base)<L
        measure_base_seq = [p.Results.measure_base,randi(2,1,L-length(p.Results.measure_base))-1];
    end
    for i = 1:L
        if measure_base_seq(i)~=base_seq(i)
            % wrong base
            result_seq(i) = -2;
        elseif number_seq(i)>0
            detected = randsrc(1,number_seq(i),[0,1,-1;1-ita,ita*(1-e_detector),ita*e_detector]);
            if ~isempty(find(detected==-1, 1))
                % wrong detected 
                result_seq(i) = mod(bit_seq(i)+1,2);
            elseif isempty(find(detected==1, 1))
                % no detected
                % consider dark count
                result_seq(i) = randsrc(1,1,[bit_seq(i),mod(bit_seq(i)+1,2),-1;y0*0.5,y0*0.5,1-y0]);
            else
                % right detected
                % consider dark count
                result_seq(i) = randsrc(1,1,[mod(bit_seq(i)+1,2),bit_seq(i);y0*0.5,1-y0*0.5]);
            end
        else 
            % no photon
            % consider dark count
            result_seq(i) = randsrc(1,1,[bit_seq(i),mod(bit_seq(i)+1,2),-1;y0*0.5,y0*0.5,1-y0]);
        end
    end
end