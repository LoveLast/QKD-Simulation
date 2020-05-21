function [Y11,e11] = evaluate_Y11_and_e11_simp(Q,E,mu,nu,max_order,TolCon,deltaQ_perc)   
    options = optimoptions('linprog','Display','iter','Preprocess','basic');
    times1 = 10^(-floor(log10(Q(1,1))));
    if times1>1e-3/TolCon
        times1 = 1e-3/TolCon;
    end
    Q = Q*times1;
    options.ConstraintTolerance = TolCon*times1;

    % min Y11 with linear programming
    num_mu = length(mu);
    num_nu = length(nu);
    num_var = (max_order+1)^2;
    c = zeros(num_var*2,1);
    c(max_order+3) = 1;% the parameter of Y11
    lb = zeros(num_var*2,1);
    ub = ones(num_var*2,1);
    Aeq = zeros(num_mu*num_nu*2,num_var*2);
    beq = zeros(num_mu*num_nu*2,1);
    A = zeros(num_var,num_var*2);
    b = zeros(num_var,1);
    for m = 1:num_mu
        for n = 1:num_nu
            beq((m-1)*num_mu+n) = Q(m,n);
            beq(num_mu*num_nu+(m-1)*num_mu+n) = Q(m,n)*E(m,n);
            for i = 0:max_order
                for j = 0:max_order
                    Aeq((m-1)*num_mu+n,i*(max_order+1)+j+1) = mu(m)^i*nu(n)^j/...
                        (factorial(i)*factorial(j))*exp(-mu(m)-nu(n));
                end
            end
        end
    end
    Aeq(num_mu*num_nu+1:end,num_var+1:end) = Aeq(1:num_mu*num_nu,1:num_var);
    for i = 1:num_var
        %q - Y<=0
        A(i,i) = -1;
        A(i,i+num_var) = 1;
    end
    
    [Y,Y11] = linprog(c,A,b,Aeq,beq,lb,ub*times1,options);
    if isempty(Y11)
        A = [A;Aeq;-Aeq];
        b = [b;beq.*(1+deltaQ_perc);-beq.*(1-deltaQ_perc)];
        [Y,Y11] = linprog(c,A,b,[],[],lb,ub*times1,options);
        Y11 = Y11/times1;
        % max e_11 with linear programming
        c(max_order+3) = 0;
        c(max_order+3+num_var) = -1;

        [q,q11] = linprog(c,A,b,[],[],lb,ub*times1,options);
        e11 = -q11/times1/Y11;
    else
        Y11 = Y11/times1;
        % max e_11 with linear programming
        c(max_order+3) = 0;
        c(max_order+3+num_var) = -1;
        
        [q,q11] = linprog(c,A,b,Aeq,beq,lb,ub*times1,options);
        e11 = -q11/times1/Y11;
    end
end



