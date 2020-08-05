function [theta1_list,theta2_list,theta1_est,theta2_est] = SAPG(y,theta1,theta2,sig,c,Nbiter,N0,T0,t)
    theta1_list = [theta1];
    theta2_list = [theta2];
    [H,W] = size(y);
    d = H*W;  
    L = (0.95/sig)^2;
    lam = min(5/L,2);
    gamma = 0.98/(L+1/lam);
    gamma_prime = 0.98*lam;
    X = y ;
    X_bar = y;
    epsilon = 10^(-10);
    theta1_est = [];
    theta2_est = [];
    start = cputime;
    % warmup
    for i = 1:T0
        g_bar = 2*epsilon*X_bar ;
        [xbar,gbar1,gbar2] = TGVprox(X_bar,theta1_list(end),theta2_list(end),0.003,50,lam);
        X_bar = X_bar-gamma_prime*g_bar-gamma_prime*(X_bar-xbar)/lam + sqrt(2*gamma_prime)*randn(H,W);
        g = (X-y)/(sig^2) + 2*epsilon*X ;
        [x,g1,g2] = TGVprox(X,theta1_list(end),theta2_list(end),0.003,50,lam);
        X = X-gamma*g-gamma*(X-x)/lam + sqrt(2*gamma)*randn(H,W);
    end
    fprintf("Fin de warm-up initialisations, temps : %f \n",cputime-start);
    n = 1;
    t1 = false;
    t2  = false;
    while  n<=Nbiter & t1 == false & t2 == false
        %thining
        for i=1:t
            g_bar = 2*epsilon*X_bar ;
            [xbar,gbar1,gbar2] = TGVprox(X_bar,theta1_list(end),theta2_list(end),0.003,50,lam);
            X_bar = X_bar-gamma_prime*g_bar-gamma_prime*(X_bar-xbar)/lam + sqrt(2*gamma_prime)*randn(H,W);
        end

        g = (X-y)/(sig^2) + 2*epsilon*X ;
        [x,g1,g2] = TGVprox(X,theta1_list(end),theta2_list(end),0.003,50,lam);
        X = X-gamma*g-gamma*(X-x)/lam + sqrt(2*gamma)*randn(H,W);
        delta2 = 20*((n)^(-0.8))/d;
        delta1 = 20*((n)^(-0.8))/d;
        theta1_list = [theta1_list,max(1e-3,theta1_list(end)+delta1*(gbar1-g1))];
        theta2_list = [theta2_list,max(1e-3,theta2_list(end)+delta2*(gbar2-g2))];
        if mod(n,30)==0
            fprintf('nb iter:%4d theta1 : %f  theta2 : %f g1bar-g1 : %f, g2bar-g2 : %f, temps : %f  \n',n,theta1_list(end),theta2_list(end),gbar1-g1,gbar2-g2,cputime-start);
        end
        if n>=N0
            theta1_est = [theta1_est,mean(theta1_list(20:end))];
            theta2_est = [theta2_est,mean(theta2_list(20:end))];
            if size(theta1_est,2) >= 2
                if c == 2
                    if max(abs(theta1_est(end)-theta1_est(end-1)),abs(theta2_est(end)-theta2_est(end-1)))/max(abs(theta1_est(end-1)),abs(theta2_est(end-1))) < 10^(-4)
                        t1 = true;
                        fprintf("Critère 2 à  l'itération %f : \n",n);
                    end
                end
                if c == 3
                    if max(abs(theta1_est(end)-theta1_est(end-1)),abs(theta2_est(end)-theta2_est(end-1))) < 10^(-3)
                        t2 = true;
                        fprintf("Critère 2 à  l'itération %f : \n",n);
                    end
                end
            end
        end

        n = n+1;
    end
end

