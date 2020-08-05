function [x,g1,g2] = TGVprox(y,theta1,theta2,tau,Nbiter,lam)

	rho = 1.99;		
	sigma = 1/tau/72; 
	[H,W]=size(y);

	
	opD = @(x) cat(3,[diff(x,1,1);zeros(1,W)],[diff(x,1,2) zeros(H,1)]);
	opDadj = @(u) [-u(1,:,1);-diff(u(1:end-1,:,1),1,1);u(end-1,:,1)]+...
		[-u(:,1,2) -diff(u(:,1:end-1,2),1,2) u(:,end-1,2)];	
	opJ = @(r) cat(3,...
		[r(1,:,1);diff(r(:,:,1),1,1)],[diff(r(:,:,1),1,2) zeros(H,1)],...
		[r(:,1,2) diff(r(:,:,2),1,2)],[diff(r(:,:,2),1,1);zeros(1,W)]);
	opJadj = @(u) cat(3,...
		[-diff(u(:,:,1),1,1);u(end,:,1)]-[u(:,1,2) diff(u(:,:,2),1,2)],...
		[-diff(u(:,:,3),1,2) u(:,end,3)]-[u(1,:,4);diff(u(:,:,4),1,1)]);		
	prox_tau_fx = @(x) (x+tau*y/lam)/(1+tau/lam);
	prox_tau_fr = @(r) r-bsxfun(@rdivide,r,max(sqrt(sum(r.^2,3))/...
		(tau*theta1),1));
	prox_sigma_g_conj = @(u) bsxfun(@rdivide,u,max(sqrt(sum(u.^2,3))/...
		theta2,1));
	
	x2 = y; 		
	r2 = zeros([H,W,2]);	
	u2 = zeros([H,W,4]);	
% Décommenter les lignes commentées si test de convergence
%     cost = [];
	for iter = 0:Nbiter-1
		tmp = tau*opJadj(u2);
		x = prox_tau_fx(x2-opDadj(tmp));
		r = prox_tau_fr(r2+tmp);
		u = prox_sigma_g_conj(u2+sigma*opJ(opD(2*x-x2)-(2*r-r2)));
		x2 = x2+rho*(x-x2);
		r2 = r2+rho*(r-r2);
		u2 = u2+rho*(u-u2);
%         g1 = sum(sum(sqrt(sum(r.^2,3))));
%         g2 = sum(sum(sqrt(sum(opJ(opD(x)-r).^2,3))));
%         cost = [cost,norm(x-y,'fro')^2/(2*lam)+theta1*sum(sum(sqrt(sum(r.^2,3))))+theta2*sum(sum(sqrt(sum(opJ(opD(x)-r).^2,3))))];
%         if mod(iter,20)==0
%             fprintf('iter : %4d, g1 : %f, g2 : %f \n',iter,g1,g2);
%         end
    end
	g1 = sum(sum(sqrt(sum(r.^2,3))));
	g2 = sum(sum(sqrt(sum(opJ(opD(x)-r).^2,3))));
end