# Mobile-Communications-
Gunjan Mane, Just for you baby

Lost Call
________________________________________________________________________________________________________
clc;
clear all;
close all;
kk=1;
N=input("Enter N(Number of Trunks: ");
A=input("Enter A: ");
for n=1:N
    num=power(A,n)/factorial(n);
    deno=0;
    for k=0:n
        deno=deno+power(A,k)/factorial(k);
    end
    final(kk)=num/deno;
    kk=kk+1;
end
disp(final);
n=1:N
    stem(n,final());
    xlabel("N");
    ylabel("P(x)");
    hold on;
plot(n,final);
hold on;
 
 
BER AWGN
________________________________________________________________________________________________________
clc;
close all;
clear all;
bitcount=100000;
SNR=0:1:10;
for k=1:1:length(SNR)
    tote=0;
    totb=0;
    while tote<100
        rbits=round(rand(1,bitcount));
        N0=1/10^(SNR(k)/10);
        tx=-2*(rbits-0.5);
        rx=tx+sqrt(N0/2)*(randn(1,length(tx)))+i*(randn(1,length(tx)));
        rx2=rx<0;
        diff=rbits-rx2;
        tote=tote+sum(abs(diff));
        totb=totb+length(rbits);
    end
    BER(k)=tote/totb;
end
semilogy(SNR,BER,'*r');
hold on;
xlabel('Eb/N0 (dB)');
ylabel('BER');
title('SNR Vs BER when modulation with BPSK in AWGN channel')
thber=0.5*erfc(sqrt(10.^(SNR/10)));
semilogy(SNR,thber);
grid on;
legend('Simulated curve','Theoritical Curve');
 
 
PCM: Speech Encoding-Decoding
________________________________________________________________________________________________________
clc;
close all;
clear all;
n=input('Enter n for n-bit PCM');
n1=input('Enter number of samples per period');
L=2^n;
x=0:2*pi/n1:4*pi;
s=8*sin(x);
subplot(5,1,1);
plot(s);
title('Analog Signal');
grid on;
xlabel('Time');
ylabel('Amplitude');
subplot(5,1,2);
stem(s);
grid on;
title('Sampled Signal');
xlabel('Time');
ylabel('Amplitude');
vmax=8;
vmin=-vmax;
del=(vmax-vmin)/L;
part=vmin:del:vmax;
code=vmin-(del/2):del:vmax+(del/2);
[ind,q]=quantiz(s,part,code);
l1=length(ind);
l2=length(q);
for i=1:l1
    if(ind(i)~=0)
        ind(i)=ind(i)-1;
    end
    i=i+1;
end
for i=1:l2
    if(q(i)==vmin-(del/2))
        q(i)=vmin+(del/2);
    end
end
subplot(5,1,3)
stem(q);
grid on;
title('Quantized Signal');
xlabel('Time');
ylabel('Amplitude');
code=de2bi(ind,'left-msb');
k=1;
for i=1:l1
    for j=1:n
        coded(k)=code(i,j);
        j=j+1;
        k=k+1;
    end
    i=i+1;
end
subplot(5,1,4);
stairs(coded);
axis([0 100 -2 3]);
grid on;
title('Encoded Signal')
xlabel('Time');
ylabel('Amplitude');
qunt=reshape(coded,n,length(coded)/n);
index=bi2de(qunt','left-msb');
q=del*index+vmin+(del/2);
subplot(5,1,5);
plot(q);
grid on;
title('Decoded Signal');
xlabel('Time');
ylabel('Amplitude');
 
 
GMSK
________________________________________________________________________________________________________
clear all;
close all;
clc;
nrz_data=[0 1 0.5 1 0 0.5 0.5]; %sample code
Tb=1; %bit duration
BT=0.3; %BT product of filter
sps=32;
Ts=Tb/sps; %sample period
t=0.5;
Tb=(-2*Tb:Ts:2*Tb);
alpha=2*pi*BT/(sqrt(log(2)))
gauss=(alpha*(2*t-0.5))-(alpha*(2*t+0.5));
K=pi/2;
gauss=K*gauss;
nrz=upsample(nrz_data,sps);
nrz_gauss=conv(gauss,nrz);
subplot(3,1,1);
stem(nrz_data);
title('sample input data');
xlabel('time');
ylabel('amp');
subplot(3,1,2);
plot(nrz_gauss);
title('gmsk output');
xlabel('time');
ylabel('amp');
MSK = [];
x = [];
for i=1:1000
    f = i/100;  % f is frequency normalized to to 1/(bit duration)
    x = [x, f ];
    ymsk = 16/(pi^2) * (cos(6.2832 * f))^2/ (1- 16 * f^2)^2;
    MSK = [MSK, 10 * log10(ymsk)];
end
subplot(3,1,3);
plot(x,MSK, 'r-');
axis([0 10 -60 10]);
title('GMSK PSD');
ylabel('Spectral Power Level in dB');
xlabel('Frequency Offset / Bit Rate');
 
 
HATA
________________________________________________________________________________________________________
clc;
close all;
clear all;
d=1:0.01:20;
hm=5;
hb1=30;
hb2=100;
hb3=200;
fc=1000;
abh=3.2*(log10(11.75*hm)).^2-4.97;

L50urban1=69.55+26.16*(log10(fc))+(44.9-6.55*log10(hb1))*log10(d)-13.82*log10(hb1)-abh;
L50urban2=69.55+26.16*(log10(fc))+(44.9-6.55*log10(hb2))*log10(d)-13.82*log10(hb2)-abh;
L50urban3=69.55+26.16*(log10(fc))+(44.9-6.55*log10(hb3))*log10(d)-13.82*log10(hb3)-abh;

L50suburban1=L50urban1-2*(log10(fc/28)).^2-5.4;
L50suburban2=L50urban2-2*(log10(fc/28)).^2-5.4;
L50suburban3=L50urban3-2*(log10(fc/28)).^2-5.4;

L50rural1=L50urban1-4.78*(log10(fc)).^2+18.33*(log10(fc))-40.94;
L50rural2=L50urban2-4.78*(log10(fc)).^2+18.33*(log10(fc))-40.94;
L50rural3=L50urban3-4.78*(log10(fc)).^2+18.33*(log10(fc))-40.94;

figure(1);
plot(d,L50urban1,'-r',d,L50urban2,'--r',d,L50urban3,':r');
hold on;
legend('hb1=30','hb2=100','hb3=200');
title('Urban');
grid on;
xlabel('d(km)');
ylabel('L(dB)');

figure(2);
plot(d,L50suburban1,'-g',d,L50suburban2,'--g',d,L50suburban3,':g');
hold on;
legend('hb1=30','hb2=100','hb3=200');
title('Suburban');
grid on;
xlabel('d(km)');
ylabel('L(dB)');

figure(3);
plot(d,L50urban1,'-r',d,L50urban2,'--r',d,L50urban3,':r');
hold on;
legend('hb1=30','hb2=100','hb3=200');
title('Rural');
grid on;
xlabel('d(km)');
ylabel('L(dB)');

