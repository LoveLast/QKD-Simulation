function [Y11,e11] = evaluate_Y11_and_e11(Q,E,mu,nu,max_order,max_iter_step,TolX,TolCon,deltaQ_perc,deltaE_perc)   
    options = optimoptions('linprog','Display','iter','Preprocess','basic');
    disp("solving for Y11...");
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
    c = zeros(num_var,1);
    c(max_order+3) = 1;% the parameter of Y11
    lb = zeros(num_var,1);
    ub = ones(num_var,1);
    Aeq = zeros(num_mu*num_nu,num_var);
    beq = zeros(num_mu*num_nu,1);
    for m = 1:num_mu
        for n = 1:num_nu
            beq((m-1)*num_mu+n) = Q(m,n);
            for i = 0:max_order
                for j = 0:max_order
                    Aeq((m-1)*num_mu+n,i*(max_order+1)+j+1) = mu(m)^i*nu(n)^j/...
                        (factorial(i)*factorial(j))*exp(-mu(m)-nu(n));
                end
            end
        end
    end

    [Y,Y11] = linprog(c,[],[],Aeq,beq,lb,ub*times1,options);
    if isempty(Y11)
        A = [Aeq;-Aeq];
        b = [beq.*(1+deltaQ_perc);-beq.*(1-deltaQ_perc)];
        [Y,Y11] = linprog(c,A,b,[],[],lb,ub*times1,options);
    end
    Y11 = Y11/times1;

    % max e11 with linear programming
%     options = optimoptions('linprog','Display','iter','Algorithm','dual-simplex',...
%                        'Preprocess','none');
%     options.ConstraintTolerance = TolCon;
%     delta=options.ConstraintTolerance/2;
%     options.ConstraintTolerance=delta;
    disp("solving for e11...");
    times2 = 10^(-floor(log10(E(1,1))));
    if times2>1e-3/TolCon/times1
        times2 = 1e-3/TolCon/times1;
    end
    E = E*times2;
    times = times1*times2;
    options.ConstraintTolerance = TolCon*times;
    
    
    c = -c;
    ub(max_order+3) = Y11;% e11<=1 => e11*y11<=y11
    beq1 = zeros(num_mu*num_nu,1);
    for m = 1:num_mu
        for n = 1:num_nu
            beq1((m-1)*num_mu+n) = Q(m,n)*E(m,n);
        end
    end

    [e,e11] = linprog(c,[],[],Aeq,beq1,lb,ub*times,options);
    if isempty(e11)
        A = [Aeq;-Aeq];
        b = [beq1.*(1+deltaE_perc);-beq1.*(1-deltaE_perc)];
        [e,e11] = linprog(c,A,b,[],[],lb,ub*times,options);
    end
    e11 = -e11/times/Y11;
    
    
    
%     if isempty(e11)
%         % 当E_z(1,1)过小时需要忽视平凡约束
%         if mu(end)==0
%             num_mu = num_mu-1;
%         end
%         if nu(end)==0
%             num_nu = num_nu-1;
%         end
%         Aeq1 = zeros(num_mu*num_nu,num_var);
%         beq1 = zeros(num_mu*num_nu,1);
%         for m = 1:num_mu
%             for n = 1:num_nu
%                 beq1((m-1)*num_mu+n) = Q(m,n)*E(m,n);
%                 for i = 0:max_order
%                     for j = 0:max_order
%                         Aeq1((m-1)*num_mu+n,i*(max_order+1)+j+1) = mu(m)^i*nu(n)^j/...
%                             (factorial(i)*factorial(j))*exp(-mu(m)-nu(n));
%                     end
%                 end
%             end
%         end
%         [e,e11] = linprog(c,[],[],Aeq1,beq1,lb,ub*times,options);
%         if isempty(e11)
%             A = [Aeq1;-Aeq1];
%             b = [beq1*(1+deltaE_perc);-beq1*(1-deltaE_perc)];
%             [e,e11] = linprog(c,A,b,[],[],lb,ub*times,options);
%         end
%         e11 = -e11/times/Y11;
%     end
    
    if isempty(e11)
        % 使用全局优化
        Q = Q/times1;
        E = E/times2;
        option = optimset('MaxFunEvals',max_iter_step,'MaxIter',max_iter_step,'TolX',TolX,'TolCon',TolCon);

        [e,e11,exitflag] = fmincon(@(t) target_function(t,max_order),rand(2*num_var,1),...
            [],[],[Aeq,0*Aeq],beq,zeros(2*num_var,1),ones(2*num_var,1),@(v) constraint_function_eq(v,Aeq,Q,E),option);
        if exitflag ~= 1
            A = [Aeq;-Aeq];
            b = [beq.*(1+deltaQ_perc);-beq.*(1-deltaQ_perc)];
            [e,e11,exitflag] = fmincon(@(t) target_function(t,max_order),rand(2*num_var,1),...
                [A,0*A],b,[],[],zeros(2*num_var,1),ones(2*num_var,1),@(v) constraint_function_neq(v,Aeq,Q,E,deltaE_perc),option);
        end
        e11 = -e11;
    end
end

function opt = target_function(vars,max_order)
    opt = -vars((max_order+1)^2+max_order+3);
end

function [res_neq,res_eq] = constraint_function_eq(vars,Aeq,Q,E)
    num_var = size(Aeq,2);
    num_mu = size(Q,1);
    num_nu = size(Q,2);
    res_neq = 0;
    res_eq = zeros(size(Aeq,1),1);

    for m = 1:num_mu
        for n = 1:num_nu
            res_eq((m-1)*num_mu+n) = -E(m,n)*Q(m,n);
        end
    end
    res_eq=res_eq+Aeq*(vars(1:num_var).*vars(num_var+1:end));
    
end

function [res_neq,res_eq] = constraint_function_neq(vars,Aeq,Q,E,delta)
    num_var = size(Aeq,2);
    num_mu = size(Q,1);
    num_nu = size(Q,2);
    res_eq = 0;
    res_neq = zeros(size(Aeq,1),1);

    for m = 1:num_mu
        for n = 1:num_nu
            res_neq((m-1)*num_mu+n) = -E(m,n)*Q(m,n);
        end
    end
    res_neq = [res_neq.*(1+delta);-res_neq.*(1-delta)];
    res_neq=res_neq+[Aeq*(vars(1:num_var).*vars(num_var+1:end));-Aeq*(vars(1:num_var).*vars(num_var+1:end))];
    
end