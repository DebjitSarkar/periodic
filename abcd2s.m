function s = abcd2s(abcd, Z0)

A = abcd(1,1);
B = abcd(1,2);
C = abcd(2,1);
D = abcd(2,2);

denom = A + B / Z0 + C * Z0 + D;
S11 = (A + B / Z0 - C * Z0 - D) / denom;
S12 = 2 * (A * D - B * C) / denom;
S21 = 2 / denom;
S22 = (-A + B / Z0 - C * Z0 + D) / denom;

s = [S11, S12; S21, S22];
