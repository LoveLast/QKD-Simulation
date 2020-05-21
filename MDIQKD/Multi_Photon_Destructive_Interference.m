function States = Multi_Photon_Destructive_Interference(M,N)
% A simulation of multi-photon Destructive Interference
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
                A_Res(i+1,u+1,v+1) = nchoosek(M,i)*nchoosek(i,u)*nchoosek(M-i,v);
            end
        end
    end
    for j = 0:N
        for s = 0:j
            for t = 0:N-j
                B_Res(j+1,s+1,t+1) = nchoosek(N,j)*nchoosek(j,s)*nchoosek(N-j,t)*(-1)^(j+s+t);
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
                            str_key = [num2str(r0),',',num2str(r1),',',num2str(s0),',',num2str(s1)];
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