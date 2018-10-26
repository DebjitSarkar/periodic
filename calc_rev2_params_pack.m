clear all
%% Parameters
dmv=linspace(.5,.9,100);
Ncasmin=4;
Ncasmax=12;
tol=1e-2;
max_iter=1e5;
Z0=50;

%% Derived Parameters
Ncasv=Ncasmin:Ncasmax;
Nc_str=string(Ncasv);
for ind=1:length(Ncasv)
    Nc_str(ind)=sprintf('N = %u',Ncasv(ind));
end

%% Loop through netwton raphson algorithm
x1m=zeros(length(Ncasv),length(dmv));
x2m=zeros(length(Ncasv),length(dmv));
tic
fprintf('Computing Ncas = ')
for nc=1:length(Ncasv)
    fprintf('%u, ',Ncasv(nc));
    [x1,x2]=calc_rev2_b0Z0(tol,max_iter,dmv,Ncasv(nc));
    x1m(nc,:)=x1;
    x2m(nc,:)=x2;
end
fprintf('\n');
toc
del1m=atan(1./x1m)./2./pi;
del2m=atan(1./x2m)./2./pi;

%% Loop Through Ncas to check min power
rn=linspace(-1,1,150);
minP_1=zeros(length(Ncasv),length(dmv));
minP_2=zeros(length(Ncasv),length(dmv));
tic
for nc=1:length(Ncasv)
    [NDM,RN]=meshgrid(1:length(dmv),rn);
    [DM,~]=meshgrid(dmv,rn);
    R=DM.*RN;
    % small b0Z0
    b0Z0_1=x1m(nc,:);
    beZ0_1=b0Z0_1(NDM).*R./(1+(1+R).*b0Z0_1(NDM).^2);
    dPHI_1=asin(beZ0_1);
    PHI_1=pi/2+dPHI_1;
    YBn_1=sqrt(1-(beZ0_1).^2);
    S21_man_cas=zeros(size(NDM));
    T_1=Ncasv(nc).*PHI_1;
    ct1=cos(T_1); st1=sin(T_1);
    S21_1=2.*YBn_1./(2.*YBn_1.*ct1+1j.*(1+YBn_1.^2).*st1);
    S21sqr_1=abs(S21_1).^2;
    minP_1(nc,:)=min(S21sqr_1,[],1);
    % Large b0Z0
    b0Z0_2=x2m(nc,:);
    beZ0_2=b0Z0_2(NDM).*R./(1+(1+R).*b0Z0_2(NDM).^2);
    dPHI_2=asin(beZ0_2);
    PHI_2=pi/2+dPHI_2;
    YBn_2=sqrt(1-(beZ0_2).^2);
    T_2=Ncasv(nc).*PHI_2;
    ct2=cos(T_2); st2=sin(T_2);
    S21_2=2.*YBn_2./(2.*YBn_2.*ct2+1j.*(1+YBn_2.^2).*st2);
    S21sqr_2=abs(S21_2).^2;  
    minP_2(nc,:)=min(S21sqr_2,[],1);  
end
toc
minPdB_1=10.*log10(minP_1);
minPdB_2=10.*log10(minP_2);

%% Plotting
figure(1)
subplot(2,2,1) % small bZ
hold off
plot(dmv,x1m)
grid on
ylabel('b0Z0')
xlabel('\delta_m')
ylim([0,2])
xlim([dmv(1),dmv(end)])
legend(Nc_str,'Location','northwest')
title('Smaller b0Z0 Parameters')
subplot(2,2,2) % large bZ
hold off
plot(dmv,x2m)
grid on
ylabel('b0Z0')
xlabel('\delta_m')
ylim([1.5,15])
xlim([dmv(1),dmv(end)])
legend(Nc_str,'Location','northwest')
title('Larger b0Z0 Parameters')
subplot(2,2,3) % del for small bZ 
hold off
plot(dmv,del1m)
grid on
ylabel('\Delta [fractional wavelengths]')
xlabel('\delta_m')
xlim([dmv(1),dmv(end)])
legend(Nc_str,'Location','northwest')
ylim([.08,.24])
subplot(2,2,4) % del for large bZ 
hold off
plot(dmv,del2m)
grid on
ylabel('\Delta [fractional wavelengths]')
xlabel('\delta_m')
xlim([dmv(1),dmv(end)])
legend(Nc_str,'Location','northeast')
ylim([0,.12])

figure(2)
subplot(2,2,1) % small bZ
hold off
plot(dmv,x1m)
grid on
ylabel('b0Z0')
xlabel('\delta_m')
ylim([0,1.5])
xlim([dmv(1),dmv(end)])
legend(Nc_str,'Location','northeast')
title('Smaller b0Z0 Parameters')
subplot(2,2,2) % large bZ
hold off
plot(dmv,x2m)
grid on
ylabel('b0Z0')
xlabel('\delta_m')
ylim([1.5,15])
xlim([dmv(1),dmv(end)])
legend(Nc_str,'Location','northwest')
title('Larger b0Z0 Parameters')
subplot(2,2,3)
hold off
plot(dmv,minPdB_1)
ylim([-.5,0])
grid on
legend(Nc_str,'Location','southwest')
title('Smaller b0Z0 Parameters')
ylabel('Minimum Power [dB]')
xlabel('\delta_m')
xlim([dmv(1),dmv(end)])
subplot(2,2,4)
hold off
plot(dmv,minPdB_2)
ylim([-.5,0])
legend(Nc_str,'Location','southwest')
title('Larger b0Z0 Parameters')
ylabel('Minimum Power [dB]')
xlabel('\delta_m')
xlim([dmv(1),dmv(end)])
grid on

function [x1,x2]=calc_rev2_b0Z0(tol,max_iter,dm,Ncas)
%% Inputs
% tol=1e-3;
% dm=linspace(.2,.9,5);
% Ncas=4;
% max_iter=1e4;
x0_1=.2;
x0_2=10;
alpha=.125; % relaxation parameter

%% Derived Parameters
opd=(1+dm);
omd=(1-dm);

%% Newton-Raphson
x=x0_1;
% er1=zeros(1,max_iter);
for n=1:max_iter
    f=asin(dm.*x./(1+opd.*x.^2))+asin(dm.*x./(1+omd.*x.^2))-2*pi/Ncas;
    er=max(abs(f(~(isnan(f)|isinf(f)))));
    if isempty(er)
        break
    end
%     er1(n)=er;
    if er<tol
        break
    end
    fap=    dm.*(1-opd.*x.^2) ...
            .*(1+opd.*x.^2-dm.*x).^-0.5 ...
            .*(1+opd.*x.^2).^-1.5;
    fbp=    dm.*(1-omd.*x.^2) ...
            .*(1+omd.*x.^2-dm.*x).^-0.5 ...
            .*(1+omd.*x.^2).^-1.5;
    fp=fap+fbp;
    x=x-alpha.*f./fp;
end
x1=real(x);
x1(x1<2*tol)=NaN;

x=x0_2;
% er2=zeros(1,max_iter);
for n=1:max_iter
    f=asin(dm.*x./(1+opd.*x.^2))+asin(dm.*x./(1+omd.*x.^2))-2*pi/Ncas;
    er=max(abs(f(~(isnan(f)|isinf(f)))));
    if isempty(er)
        break
    end
%     er2(n)=er;
    if er<tol
        break
    end
    fap=    dm.*(1-opd.*x.^2) ...
            .*(1+opd.*x.^2-dm.*x).^-0.5 ...
            .*(1+opd.*x.^2).^-1.5;
    fbp=    dm.*(1-omd.*x.^2) ...
            .*(1+omd.*x.^2-dm.*x).^-0.5 ...
            .*(1+omd.*x.^2).^-1.5;
    fp=fap+fbp;
    x=x-alpha.*f./fp;
end
x2=real(x);
x2(x2<2*tol)=NaN;
%% Check Converged Value
% Z0b1_1=x1.*dm./(1+opd.*x1.^2);
% Z0b2_1=-x1.*dm./(1+omd.*x1.^2);
% dphi_1=Ncas.*(asin(Z0b1_1)-asin(Z0b2_1)).*180/pi
% Z0b1_2=x2.*dm./(1+opd.*x2.^2);
% Z0b2_2=-x2.*dm./(1+omd.*x2.^2);
% dphi_2=Ncas.*(asin(Z0b1_2)-asin(Z0b2_2)).*180/pi
% 
% %% Calculate fractional length Required
% del1=atan(1./x1)./2./pi;
% del2=atan(1./x2)./2./pi;

% %% Plotting
% st=find(~(isnan(x2)|isinf(x2)),1);
% figure(1)
% hold off
% plot(dm,x1)
% hold on
% plot(dm,x2)
% xlim([dm(st),dm(end)])
% 
% figure(2)
% hold off
% semilogy(er1)
% hold on
% semilogy(er2)
% ylim([tol,10])
% 
% figure(3)
% hold off
% plot(dm,del1)
% hold on
% plot(dm,del2)
% xlim([dm(st),dm(end)])
end
