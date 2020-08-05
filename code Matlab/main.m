function [y,x,theta1,theta2,theta1_est,theta2_est] = main(im_name,SNR,patch)

	N= 600;	

	y  = double(imread(strcat('../im/',im_name,'.png')))/255;
    a = size(size(y));
    if a(2)== 3 
        y = rgb2gray(y);
    end
    figure(1);
	imshow(y);
	rng(0);
    
    sig = mean(y,'all') * (10^(-SNR/10));
	y = y+randn(size(y))*sig; 
	figure(2);
	imshow(y);
    
    y_prime = y(patch(1):patch(2),patch(3):patch(4));
    N0 = 20;
    T0 = 25;
    t = 6;
    [theta1,theta2,theta1_est,theta2_est] = SAPG(y_prime,10,10,sig,1,N,N0,T0,t);
    figure(3);
    plot(theta1,'r');
    hold on;
    plot(theta2,'b');
    title(strcat(im_name,' ',int2str(SNR),'db ',int2str(N),' itérations'))
    legend('theta1','theta2');
    grid on;
    
    figure(4);
    a1 = abs(diff(theta1)./theta1(1:end-1));
    a2 = abs(diff(theta2)./theta2(1:end-1));
    plot(a1,'r');
    hold on;
    plot(a2,'b');
    title(strcat(im_name,' ',int2str(SNR),'db ',int2str(N),' itérations'))
    legend('vr theta1','vr theta2');
    grid on;
    
    figure(5);
    b1 = abs(diff(theta1_est)./theta1_est(1:end-1));
    b2 = abs(diff(theta2_est)./theta2_est(1:end-1));
    plot(b1,'r');
    hold on;
    plot(b2,'b');
    title(strcat(im_name,' ',int2str(SNR),'db ',int2str(N),' itérations'))
    legend('vr theta1 à partir de N0','vr theta2 à partir de N0');
    grid on;
    
   
    
    
    
	x = TGVdenoising(y,theta1_est(end),theta2_est(end),sig,100);
	figure(6);
	imshow(x);
    % à décommenter si on veut sauvegarder les images, changer le nom du
    % dossier si nécessaire
	%imwrite(y,strcat('../im_noise/',im_name,int2str(SNR),'2.png'));
	%imwrite(x,strcat('../im_denoised/',im_name,int2str(SNR),'_2d.png'));
end