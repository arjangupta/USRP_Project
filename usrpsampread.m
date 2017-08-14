close all; clc;

% Author: Arjan Gupta
% Description: USRP rx packet samples analysis

% This code reads and plots the USRP samples
fh = fopen('usrp_samples2412.dat', 'rb');
data = fread(fh,'int16');
cdata = complex(data(1:2:end),data(2:2:end));
fclose(fh);
re = real(cdata);

%figure
%plot(re)

%figure
%plot(xcorr(re,re));

% we need to generate an 'ideal' tx signal

% expected short preamble
spreamb = [zeros(1,8), (1+1i), zeros(1,3), (-1-1i), zeros(1,3)];
spreamb = [spreamb, (1+1i), zeros(1,3), (-1-1i), zeros(1,3), (-1-1i), zeros(1,3), (1+1i), zeros(1,3)];
spreamb = [spreamb, zeros(1,4), (-1-1i), zeros(1,3), (-1-1i), zeros(1,3), (1+1i), zeros(1,3)];
spreamb = [spreamb, (1+1i), zeros(1,3), (1+1i), zeros(1,3), (1+1i), zeros(1,7)];
spreamb = ((13/6)^(1/2)).*spreamb;

% inverse fast fourier transform of spreamb
ifft_spreamb = ifft(spreamb);

% extend ifft_spreamb
ifft_spreamb = [ifft_spreamb, ifft_spreamb, ifft_spreamb(1:33)];

% multiply ifft_spreamb by window function
ifft_spreamb(1) = 0.5*ifft_spreamb(1);
ifft_spreamb(161) = 0.5*ifft_spreamb(161);

% all even values of ifft_spreamb seem to be negative, so let us rectify that
for k = 1:length(ifft_spreamb)
    if (mod(k,2) == 0)
       ifft_spreamb(k) = ifft_spreamb(k)*-1; 
    end
end

% expected long preamble
lpreamb = [zeros(1,6), ones(1,2), -1.*ones(1,2), ones(1,2), -1, 1, -1, 1];
lpreamb = [lpreamb, ones(1,5), -1.*ones(1,2), ones(1,2), -1, 1, -1, ones(1,4)];
lpreamb = [lpreamb, 0, 1, -1.*ones(1,2), ones(1,2), -1, 1, -1, 1, -1.*ones(1,5), 1];
lpreamb = [lpreamb, 1, -1.*ones(1,2), 1, -1, 1, -1, ones(1,4), zeros(1,5)];

% inverse fast fourier transform of lpreamb
ifft_lpreamb = ifft(lpreamb);

% put the second half of the lpreamb in front
ifft_lpreamb = [ifft_lpreamb(33:end), ifft_lpreamb(1:32)];
% corred = xcorr(cdata,ifft_lpreamb);
% abs_corred = abs(corred);

% extend ifft_spreamb
ifft_lpreamb = [ifft_lpreamb, ifft_lpreamb, ifft_lpreamb(1:33)];

% multiply ifft_lpreamb by window function
ifft_lpreamb(1) = 0.5*ifft_lpreamb(1);
ifft_lpreamb(161) = 0.5*ifft_lpreamb(161);

% all even values of ifft_lpreamb seem to be negative, so let us rectify
% that
for k = 1:length(ifft_lpreamb)
    if (mod(k,2) == 0)
       ifft_lpreamb(k) = -1*ifft_lpreamb(k); 
    end
end

preamble = [ifft_spreamb(1:160), (ifft_spreamb(161)+ifft_lpreamb(1)), ifft_lpreamb(2:160)];

% corred = xcorr(cdata,ifft_lpreamb);
% abs_corred = abs(corred);
% 
% figure
% plot(abs_corred);
% abs_cdata = abs(cdata);
% cdata500 = cdata(2901200:2901700); 

% figure
% subplot(2,1,1)
% plot(real(cdata500));
% subplot(2,1,2)
% plot(imag(cdata500));

% ideal preamble in 802.11!
% figure
% abs_preamb = abs(preamble);
% plot(1:320,abs_preamb);
% legend('Ideal preamble in 802.11');

figure
plot(xcorr(abs(cdata),abs(preamble)));