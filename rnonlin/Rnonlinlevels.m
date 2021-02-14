function [yl2l, yr2l, maxl, maxr] = Rnonlinlevels (yl2, yr2, fs)

framesize = fs / 1000 * 30;
c = 1;
for b = 1 : framesize : length(yl2)-framesize
    for a = 1 : 1 : 40
        yl2l(a,c) = 10*log10((1/framesize)*sum(yl2(a,b:b+framesize).*yl2(a,b:b+framesize)));
        
    end
    c = c + 1;
end
c = 1;
for b = 1 : framesize : length(yr2)-framesize
    for a = 1 : 1 : 40
        yr2l(a,c) = 10*log10((1/framesize)*sum(yr2(a,b:b+framesize).*yr2(a,b:b+framesize)));
       
    end
    c = c + 1;
end
for a = 1 : 1 : length(yl2l)
maxl(a,1) = max(yl2l(:,a));
maxr(a,1) = max(yr2l(:,a));
end
end