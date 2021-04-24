function xcorrw = rnonlinweighting(levels, maxlevel, xcorrtobeweighted)
% 
maxcoeffsave = 0;
coeffacum = 0;
c = 1;
offsetaccume = 0;
numacc = 0;
numacc2 = 0;
[ds1, ds2] = size(xcorrtobeweighted);
tsize = ds1 * ds2;
coeffsum = 0;
for a = 1 : 1 : length(xcorrtobeweighted)
    for b = 1 : 1 : 40
        if levels(b,a) < (maxlevel(c,1) - 80)
            xcorrw(b,a) = 0;
            numacc2 = numacc2 + 1;
            coeffacum(b,a) = 0;
        end
        if (levels(b,a) <= ((maxlevel(c,1) - 40))) && (levels(b,a) >= ((maxlevel(c,1) - 80)))
            xcorrw(b,a) = (xcorrtobeweighted(b,a)*((1/40)*(((levels(b,a) - (maxlevel(c,1)-80)) / (((maxlevel(c,1)-40) - (maxlevel(c,1)-80)))))));
            offsetaccume = offsetaccume + (1/40 *(((levels(b,a) - (maxlevel(c,1)-80)) / (((maxlevel(c,1)-40) - (maxlevel(c,1)-80))))));
            numacc = numacc + 1; 
            coeffacum(b,a) = (1/40 *(((levels(b,a) - (maxlevel(c,1)-80)) / (((maxlevel(c,1)-40) - (maxlevel(c,1)-80))))));
        end
    end
na = 40 - numacc - numacc2;
odiff = numacc * (1/40);
xc = (odiff) - offsetaccume;
    for b = 1 : 1 : 40
        if levels(b,a) > (maxlevel(c,1) - 40)
            xcorrw(b,a) = (xcorrtobeweighted(b,a)*((1/40)+(xc/na)));
%             maxcoeffsave(b,a) = ((1/40)+(xc/na));
            coeffacum(b,a) = ((1/40)+(xc/na));
        end
    end
coeffsum(a) = (((offsetaccume/numacc) + (na*(1/40)) + (xc/na)));
% coeffacume(a) = offsetaccume;
numacc = 0;
numacc2 = 0;
offsetaccume = 0;
    c = c + 1;
end


end
