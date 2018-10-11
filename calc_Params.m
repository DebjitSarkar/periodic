% Cody Scarborough
clear all
%% Parameters
Ptol_dB=-.5; % Lowest acceptable power in the passband
f0=10e9; Z0=50; % Code does not depend on these until printing
dm=[.1,.2,.4,.8]; % Values of modulation depth to checkout
plot_me=true;

%% Derived Parameters
w0=2*pi*f0;
n=2:10;
ms=linspace(-1,1,501);
[N,DM,Ms]=meshgrid(n,dm,ms);
R=Ms.*DM;
dm_str=string(dm);
for ind=1:length(dm)
    dm_str(ind)=sprintf('\\delta_m = %0.1f',dm(ind));
end

%% Calc Center Capacitance
b0Z0=sin(pi./N)./DM;
bZ0=b0Z0.*R;

%% Determine Minimum Power
phi=acos(-bZ0); % Bloch phase progression for all modulation points
ZBn=1./sqrt(1-(bZ0).^2);
theta=N.*phi; % Cascaded phase progression
ct=cos(theta); st=sin(theta);
S21sqr=abs(2./(2.*ct+1j.*(ZBn+1./ZBn).*st)).^2;
% Handle N=2 Case (Zb blows up since inf series impedances)
S21sqr(:,(n==2),[1,length(ms)])=0.5;
S21dB=10.*log10(S21sqr);
S21dB_min=min(S21dB,[],3);
% Check all modulation reaches 2*pi
check_sum=sum(sum((abs(2*pi-(phi(:,:,end)-phi(:,:,1)).*N(:,:,1)))>1e-3));
if check_sum
    fprintf('Uh-oh, phi(dm)-phi(-dm)~=2*pi\n')
end

%% Find first N to reach power tolerance
S21dB_min_vec=S21dB_min(1,:);
lowest_N_ind=find(S21dB_min_vec>Ptol_dB,1);

%% Printing
fprintf('Lowest N to achieve power tolerance: %u\n',n(lowest_N_ind))
fprintf('Required C0 vs dm...\n')
fprintf('|')
fprintf('  dm =%10.3f  |',dm)
fprintf('\n')
fprintf('|')
fprintf('  C0 =%10.3e  |',b0Z0(:,lowest_N_ind,1)/w0/Z0)
fprintf('\n')

%% Plotting
if plot_me
figure(1)
subplot(3,1,1)
plot(n,S21dB_min_vec,'-*')
grid on
title('Minimum Power Over Modulating Range')
xlabel('N')
ylabel('S_{21} [dB]')
subplot(3,1,2:3)
plot(n,b0Z0(:,:,1),'-*')
title('Required Susceptance')
xlabel('N')
ylabel('\omega C Z_0')
legend(dm_str)
grid on

figure(2)
imagesc(ms,n,abs(reshape(S21dB(1,:,:)+3,length(n),length(ms))))
caxis([2.995,3])
xlabel('Modulation Signal')
ylabel('Number of Stages')
title('Transmitted Power of Modulation Range')
end

