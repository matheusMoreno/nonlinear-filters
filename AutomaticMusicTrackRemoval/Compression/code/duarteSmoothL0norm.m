function stillong = duartSmoothL0norm(xlong,N,P,sigmaL0)
   

[a b] = max(abs(xlong));
x = xlong(b-N/2:b+N/2-1);
[y_SL0_fmincon, w_SL0_fmincon] = fmincon_NL_comp_SL0_DCT(x,P,sigmaL0);
    

Nlong = length(xlong);
Xlong = zeros(Nlong,P);
for k=1:P
    Xlong(:,k) = [xlong.^(2*k - 1)];
end;

stillong = Xlong*w_SL0_fmincon;