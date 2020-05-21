function States = Multi_Photon_HOM(M,N,t,r)
% A simulation of multi-photon HOM effect
% return States = [Left;Possibility], Left respresents |Left,M+N-Left>

    Left = 0:M+N;
    Possibility = zeros(1,M+N+1);
    
    for m=0:M
        for n = 0:N
            index = M-m+n+1;
            Possibility(index) = Possibility(index)+(-1)^m*sqrt(factorial(M)*factorial(N)...
                *factorial(M-m+n)*factorial(N-n+m))/(factorial(M-m)*factorial(m)*factorial(N-n)*factorial(n))...
                *t^(M+N-m-n)*r^(m+n);
        end
    end

    Possibility = Possibility.^2;
    Possibility = Possibility/sum(Possibility);
    States = [Left;Possibility];
end