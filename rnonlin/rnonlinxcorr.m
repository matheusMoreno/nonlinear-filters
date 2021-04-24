function [xcorrval lagstore] = rnonlinxcorr (clean, dirty, fs)
steplength = fs/1000*30;
corrlength = fs/1000*10;
c = 1;
lagstore = zeros(length(1 : steplength : length(clean) - steplength), 1);
for a = 1 : steplength : length(clean) - steplength
    for b = 1 : 1 : 40
        [xcorrseq, lags] = xcorr(dirty(b,a : a + steplength), clean(b,a : a + steplength), corrlength, 'coeff');
        xcorrval(b,c) = max(xcorrseq);
        for d = 1 : 1 : length(xcorrseq)
            if xcorrseq(d) == xcorrval(b,c)
                lagstore(c) = lags(d);
            end
        end
    end
    c = c + 1;
end
end
