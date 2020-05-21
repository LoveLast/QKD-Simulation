function States = Multi_Photon_Interference(M,N,phase_diff)
% A simulation of multi-photon Interference with X base
% phase_diff = 1 means an different bit choice 
% return States = [Prob,Dict']
% Prob = [State_Key;Possibility], 
% Dict = {State_Key:State}, State_Key is an int and State is
% [r0,r1,s0,s1]
    state_num = nchoosek(M+N+3,3);
    Prob = zeros(2,state_num);
    Prob(1,:) = 1:state_num;
    Dict = zeros(state_num,4);
    Res = containers.Map();
    A_Res = zeros(M+1,M+1,M+1);
    B_Res = zeros(N+1,N+1,N+1);
    for i = 0:M
        for u = 0:i
            for v = 0:M-i
                A_Res(i+1,u+1,v+1) = nchoosek(M,i)*nchoosek(i,u)*nchoosek(M-i,v)*(-1)^(i);
            end
        end
    end
    for j = 0:N
        for s = 0:j
            for t = 0:N-j
                B_Res(j+1,s+1,t+1) = nchoosek(N,j)*nchoosek(j,s)*nchoosek(N-j,t);
                if phase_diff==1
                    B_Res(j+1,s+1,t+1) = B_Res(j+1,s+1,t+1)*(-1)^(s+t);
                end
            end
        end
    end
    key = 1;
    for i=0:M
        for j=0:N
            for u = 0:i
                for v = 0:M-i
                    for s = 0:j
                        for t = 0:N-j
                            r0 = M-i-v+N-j-t;
                            r1 = i-u+j-s;
                            s0 = v+t;
                            s1 = u+s;
                            str_key = [getchar(r0),',',getchar(r1),',',getchar(s0),',',getchar(s1)];
                            if ~isKey(Res,str_key)
                                Res(str_key)= key;
                                Dict(key,:) = [r0,r1,s0,s1];
                                key = key+1;
                            end
                            Prob(2,Res(str_key)) = Prob(2,Res(str_key))+A_Res(i+1,u+1,v+1)*B_Res(j+1,s+1,t+1);
                        end
                    end
                end
            end
        end
    end
    Prob(2,:) = Prob(2,:).^2;
    Prob(2,:) = Prob(2,:)./sum(Prob(2,:));
    States = [Prob;Dict'];
end