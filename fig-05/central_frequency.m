%Function to find the biggest energy between 0.18 and 3.3 Hz

function [central] = central_frequency(ECG,fs)

%Resample to 256 HZ to get a modulation frequency range between 0 and 16 Hz. 
%This will help to visualize the modulation spectral components at higher
%heart rates
if fs~=256
    ECG=resample(ECG,256,fs);
    fs=256;
end

% BASEWINDOW=32;
% BASEOVERLAP=BASEWINDOW*0.75;
% BASENFFT=512;
% MODWINDOW=128;
% MODNFFT=512;
% 
% [B,F,T] = modspecgram(ECG,fs,BASEWINDOW,BASEOVERLAP,BASENFFT,MODWINDOW,MODNFFT);
% 
% %Modulation spectrogram figure
% modspecgram(ECG,fs,BASEWINDOW,BASEOVERLAP,BASENFFT,MODWINDOW,MODNFFT);

modspect = strfft_modulation_spectrogram(ECG(1: 1048), fs, 32, 8, 16, [], 4, []);
B = sqrt(modspect.power_modulation_spectrogram);
F = modspect.freq_axis';
T = modspect.freq_mod_axis'; 
figure()
plot_modulation_spectrogram_data(modspect);

%Scale the matrix between 0 and 1
B1=B(:,11:end);
max_value=max(B1);
max_value=max(max_value);
B1=B1/max_value;
A=[F B1];
T1=T(11:end);
T1=[0 T1'];

%Modulation spectrogram matrix with the frequencies axis 
A=[T1;A];

maxima_col=0;
[m,n]=size(B1);      %number of columns

sum1=0;  
sum_max=0;

% Find where 3 Hz is:
for i=2:n
    if A(1,i)==3
       limit3=i;
       break
    end
end
 
%Find central frequency

for i=5:limit3+8
    for j=2:1:82                %number of rows
       sum1=A(j,i)+sum1;
    end
    
    if sum1>sum_max
            sum_max=sum1;
            maxima_col=i;
    end
    
    central_freq=A(1,maxima_col);
    sum1=0;      
end
MaxEnergy=sum_max;

%Biggest frequency detected between 0.18 and 3.3 Hz
central=central_freq;

%Checking if there is a bigger energy after the central energy:
centralAfter=0;
final=limit3+8;        %column with 3.3 Hz
sum2=0;
sum_max1=0;
b=maxima_col+5;
if b < final
    for i=b:final
        for j=2:1:82                %number of rows
           sum2=A(j,i)+sum2;
        end

        if sum2>sum_max1
                sum_max1=sum2;
                maxima_col1=i;
        end

        centralAfter=A(1,maxima_col1);
        sum2=0;      
    end
    AftEnergy=sum_max1; 
end

%Checking if there is a biggest energy before the central energy: 
centralBefore=0;
sum3=0;
sum_max2=0;
BefEnergy=0;
final=maxima_col-5;

    for i=5:final
        for j=2:1:82                %number of rows
           sum3=A(j,i)+sum3;
        end

        if sum3>sum_max2
                sum_max2=sum3;
                maxima_col2=i;
        end

        centralBefore=A(1,maxima_col2);
        sum3=0;      
    end 
    BefEnergy=sum_max2; 

%Check if the other detected frequencies have a bigger (> threshold) 
%energy that the central
thres=MaxEnergy*0.9;

if centralBefore ~= 0
    if BefEnergy>=thres 
        central=centralBefore;
    end
end

if centralAfter ~= 0
    if AftEnergy>=thres 
        central=centralAfter;
    end
end
 
end