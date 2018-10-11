% syms BL ZC b;
% eq1 = 0 == cos(BL) - ZC*b*sin(BL);
% eq2 = 1j*50 == 1j*ZC*sin(BL);
% eq3 = 1j/50 == 1j*b*cos(BL) + b*(cos(BL) - ZC*b*sin(BL))*1j + (sin(BL)*1j)/ZC;
% eq = [eq1; eq2; eq3];

num = 10;
Crange = linspace(1e-13, 1e-11, num);
ZCrange = linspace(10, 150, num);
BL = pi / 10;

figure;
hold on;

err = zeros(num);
for m = 1:num
    for n = 1:num
        ZC = ZCrange(n);
        b = 2 * pi * 1e10 * Crange(m);
        err1 = cos(BL) - ZC*b*sin(BL);
        err2 = 1j*ZC*sin(BL) - 1j*50;
        err3 = 1j*b*cos(BL) + b*(cos(BL) - ZC*b*sin(BL))*1j + (sin(BL)*1j)/ZC; -1j*50;
        
        err(m,n) = (err1^2 + err2^2 + err3^2)^0.5;
        scatter3(ZC, b / (2 * pi * 1e10), err(m,n));
    end
end

%%
f = 1e10;
w = 2 * pi * f;

ZC = 1;
b = 1;

num = 100;
Crange = linspace(1e-13, 1e-11, num);
ZCrange = linspace(0.1, 150, num);
BL = pi/10;
%BL = linspace(0, pi/4, num);

mag = zeros(num);
ang = mag;

for m = 1:num
    for n = 1:num
        ZC = ZCrange(m);
        b = w * Crange(n);
        
        A = cos(BL)-sin(BL)*b*ZC;
        B = 1j*sin(BL)*ZC;
        C = (1j.*b*ZC*2*cos(BL)-1j*(ZC*b).^2*sin(BL)+1j*sin(BL))/ZC;
        D = A;
        %assert(abs(A*D - B*C - 1) < 0.1);
        S21 = 2 / (A + B/ZC + C*ZC + D);
        
        mag(m,n) = abs(S21)^2;
        ang(m,n) = phase(S21);
    end
end

% X = ZCrange
% Y = Crange
% Z = mag, abs

X = repmat(ZCrange', 1, num);
Y = repmat(Crange, num, 1);

figure;
hold on;
title('Magnitude');
xlabel('ZC');
ylabel('C [F]');
surfc(X, Y, mag);

figure;
hold on;
title('Angle');
xlabel('ZC');
ylabel('C [F]');
surfc(X, Y, rad2deg(ang));



%ABCD_0(:,:,i) = [cos(BL)-sin(BL)*b(i)*Z0, 1j*sin(BL)*Z0;
%        (1j.*b(i)*Z0*2*cos(BL)-1j*(Z0*b(i)).^2*sin(BL)+1j*sin(BL))/Z0,cos(BL)-sin(BL)*b(i)*Z0];
%S21_0(i) = 2/(ABCD_0(1,1,i)+ABCD_0(1,2,i)/Z0+ABCD_0(2,1,i)*Z0+ABCD_0(2,2,i));


% scatter3(Crange, ZCrange, err);

%%
num = 100;
Y1 = linspace(-1,1, num); % Change ranges
Y2 = linspace(-1,1, num);
betad = pi / 10;
ZC = 50;
abcd_tl = [cos(betad), 1j*ZC*sin(betad); 1j/ZC*sin(betad), cos(betad)];

mag = zeros(num);
ang = mag;

for m = 1:num
    Y1val = Y1(m);
    for n = 1:num
        Y2val = Y2(n);
        
        abcd = [1 0; Y1val 1];
        abcd = abcd * abcd_tl;
        abcd = abcd * [1 0; Y2val 1];
        abcd = abcd * [1 0; Y1val 1];
        
        A = abcd(1,1);
        B = abcd(1,2);
        C = abcd(2,1);
        D = abcd(2,2);
        
        S21 = 2 / (A + B/ZC + C*ZC + D);
        mag(m,n) = abs(S21);
        ang(m,n) = angle(S21);
    end
end

X = repmat(Y1', 1, num);
Y = repmat(Y2, num, 1);

figure;
hold on;
title('Magnitude');
xlabel('Y1');
ylabel('Y2');
mesh(X, Y, mag);

figure;
hold on;
title('Angle');
xlabel('Y1');
ylabel('Y2');
mesh(X, Y, rad2deg(ang));
