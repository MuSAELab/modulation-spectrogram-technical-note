function enhanced = Enhancement_ECG(signal,fs)

%Resample the signal to 256 Hz. 
if fs~=256
    signal=resample(signal,256,fs);
    fs=256;
end

windowsize=32;          %Number of points of the window for the spectrogram
fhop=windowsize/4;      %75% overlap

%ECG signal normalization
signal = (signal-mean(signal))./std(signal);

%Heart rate detection (biggest energy between 0.18 and 3.3 Hz in the 
%modulation spectral frequency axis)
[fc] = central_frequency(signal,fs);

%Band pass filter configuration
%The width for each modulation frequency lobe is 0.625 Hz, so BW2=0.625/2
BW1=0.4;            %Stop frequency
BW2 = 0.3125;       %Pass frequency

%Number of filters according to the detected central frequency
NumFilters = fc*10+BW1;
NumFilters1 = fc*9+BW1;
NumFilters2 = fc*8+BW1;
NumFilters3 = fc*7+BW1;
NumFilters4 = fc*6+BW1;
NumFilters5 = fc*5+BW1;

if NumFilters <= 16
    NF=10;
elseif NumFilters1 <= 16
    NF=9;
elseif NumFilters2 <= 16
    NF=8;
elseif NumFilters3 <= 16
    NF=7;
elseif NumFilters4 <= 16
    NF=6;
elseif NumFilters5 <= 16
    NF=5;
end

switch NF
    
    case 10
        %Band pass filters for all the modulation spectral components
        h_band1 = bandpass_func(1*fc-BW1, 1*fc-BW2, 1*fc+BW2, 1*fc+BW1);
        h_band2 = bandpass_func(2*fc-BW1, 2*fc-BW2, 2*fc+BW2, 2*fc+BW1);
        h_band3 = bandpass_func(3*fc-BW1, 3*fc-BW2, 3*fc+BW2, 3*fc+BW1);
        h_band4 = bandpass_func(4*fc-BW1, 4*fc-BW2, 4*fc+BW2, 4*fc+BW1);
        h_band5 = bandpass_func(5*fc-BW1, 5*fc-BW2, 5*fc+BW2, 5*fc+BW1);
        h_band6 = bandpass_func(6*fc-BW1, 6*fc-BW2, 6*fc+BW2, 6*fc+BW1);
        h_band7 = bandpass_func(7*fc-BW1, 7*fc-BW2, 7*fc+BW2, 7*fc+BW1);
        h_band8 = bandpass_func(8*fc-BW1, 8*fc-BW2, 8*fc+BW2, 8*fc+BW1);
        h_band9 = bandpass_func(9*fc-BW1, 9*fc-BW2, 9*fc+BW2, 9*fc+BW1);
        h_band10 = bandpass_func(10*fc-BW1, 10*fc-BW2, 10*fc+BW2, 10*fc+BW1);

        %computing magnitude + phase (75% overlap)
        [mag,phase,fpad] = basetran_mod(signal,fhop);

        [nfreqbins, nframes] = size(mag);
        
        stepfreq=fs/nfreqbins;       %frequency step
        f40=40/stepfreq;            

        for k=1:f40       %Filtering between 0 and 40 Hz in the frequency vertical axis

           m_F1(k,:) = (filtfilt(h_band1.Numerator, 1, mag(k,:)));
           m_F2(k,:) = (filtfilt(h_band2.Numerator, 1, mag(k,:)));
           m_F3(k,:) = (filtfilt(h_band3.Numerator, 1, mag(k,:)));
           m_F4(k,:) = (filtfilt(h_band4.Numerator, 1, mag(k,:)));
           m_F5(k,:) = (filtfilt(h_band5.Numerator, 1, mag(k,:)));
           m_F6(k,:) = (filtfilt(h_band6.Numerator, 1, mag(k,:)));
           m_F7(k,:) = (filtfilt(h_band7.Numerator, 1, mag(k,:)));
           m_F8(k,:) = (filtfilt(h_band8.Numerator, 1, mag(k,:)));
           m_F9(k,:) = (filtfilt(h_band9.Numerator, 1, mag(k,:)));
           m_F10(k,:) = (filtfilt(h_band10.Numerator, 1, mag(k,:)));

           mag_filt(k,:) = m_F1(k,:) + m_F2(k,:) + m_F3(k,:) + m_F4(k,:) + m_F5(k,:) + m_F6(k,:) + m_F7(k,:) + m_F8(k,:) + m_F9(k,:) + m_F10(k,:);

        end

        %Half wave rectification to eliminate possible negative values in
        %the filtered magnitude
        for i=1:f40
            for j=1:nframes
                if mag_filt(i,j) < 0
                    mag_filt(i,j)=0;
                end
            end
        end

        %Fill with zeros the the frequencies bigger than 40 Hz
        a=windowsize; %8*fhop
        b=windowsize-f40;
        c=nframes;

        mag_filt(f40+1:a,:) = zeros(b,c);
        
        %Double size spectrum for the inverse FFT
        mag_filt_tot(1:windowsize,:) = [mag_filt(1:(windowsize/2+1),:); (mag_filt([(windowsize/2):-1:2],:))];

        %Inverse FFT to get the time domain signal
        enhanced = invbasetran_mod(mag_filt_tot,phase,fpad);
        
        %Normalization
        enhanced = (enhanced-mean(enhanced))./std(enhanced); 
   
    case 9
    
        h_band1 = bandpass_func(1*fc-BW1, 1*fc-BW2, 1*fc+BW2, 1*fc+BW1);
        h_band2 = bandpass_func(2*fc-BW1, 2*fc-BW2, 2*fc+BW2, 2*fc+BW1);
        h_band3 = bandpass_func(3*fc-BW1, 3*fc-BW2, 3*fc+BW2, 3*fc+BW1);
        h_band4 = bandpass_func(4*fc-BW1, 4*fc-BW2, 4*fc+BW2, 4*fc+BW1);
        h_band5 = bandpass_func(5*fc-BW1, 5*fc-BW2, 5*fc+BW2, 5*fc+BW1);
        h_band6 = bandpass_func(6*fc-BW1, 6*fc-BW2, 6*fc+BW2, 6*fc+BW1);
        h_band7 = bandpass_func(7*fc-BW1, 7*fc-BW2, 7*fc+BW2, 7*fc+BW1);
        h_band8 = bandpass_func(8*fc-BW1, 8*fc-BW2, 8*fc+BW2, 8*fc+BW1);
        h_band9 = bandpass_func(9*fc-BW1, 9*fc-BW2, 9*fc+BW2, 9*fc+BW1);

        [mag,phase,fpad] = basetran_mod(signal,fhop);

        [nfreqbins, nframes] = size(mag);
        stepfreq=fs/nfreqbins;
        f40=40/stepfreq;

        for k=1:f40

           m_F1(k,:) = (filtfilt(h_band1.Numerator, 1, mag(k,:)));
           m_F2(k,:) = (filtfilt(h_band2.Numerator, 1, mag(k,:)));
           m_F3(k,:) = (filtfilt(h_band3.Numerator, 1, mag(k,:)));
           m_F4(k,:) = (filtfilt(h_band4.Numerator, 1, mag(k,:)));
           m_F5(k,:) = (filtfilt(h_band5.Numerator, 1, mag(k,:)));
           m_F6(k,:) = (filtfilt(h_band6.Numerator, 1, mag(k,:)));
           m_F7(k,:) = (filtfilt(h_band7.Numerator, 1, mag(k,:)));
           m_F8(k,:) = (filtfilt(h_band8.Numerator, 1, mag(k,:)));
           m_F9(k,:) = (filtfilt(h_band9.Numerator, 1, mag(k,:)));

           mag_filt(k,:) = m_F1(k,:) + m_F2(k,:) + m_F3(k,:) + m_F4(k,:) + m_F5(k,:) + m_F6(k,:) + m_F7(k,:) + m_F8(k,:) + m_F9(k,:);

        end

        for i=1:f40
            for j=1:nframes
                if mag_filt(i,j) < 0
                    mag_filt(i,j)=0;
                end
            end
        end


        a=windowsize; 
        b=windowsize-f40;
        c=nframes;

        mag_filt(f40+1:a,:) = zeros(b,c);

        mag_filt_tot(1:windowsize,:) = [mag_filt(1:(windowsize/2+1),:); (mag_filt([(windowsize/2):-1:2],:))];

        enhanced = invbasetran_mod(mag_filt_tot,phase,fpad);
        enhanced = (enhanced-mean(enhanced))./std(enhanced);
               

   case 8
    
        h_band1 = bandpass_func(1*fc-BW1, 1*fc-BW2, 1*fc+BW2, 1*fc+BW1);
        h_band2 = bandpass_func(2*fc-BW1, 2*fc-BW2, 2*fc+BW2, 2*fc+BW1);
        h_band3 = bandpass_func(3*fc-BW1, 3*fc-BW2, 3*fc+BW2, 3*fc+BW1);
        h_band4 = bandpass_func(4*fc-BW1, 4*fc-BW2, 4*fc+BW2, 4*fc+BW1);
        h_band5 = bandpass_func(5*fc-BW1, 5*fc-BW2, 5*fc+BW2, 5*fc+BW1);
        h_band6 = bandpass_func(6*fc-BW1, 6*fc-BW2, 6*fc+BW2, 6*fc+BW1);
        h_band7 = bandpass_func(7*fc-BW1, 7*fc-BW2, 7*fc+BW2, 7*fc+BW1);
        h_band8 = bandpass_func(8*fc-BW1, 8*fc-BW2, 8*fc+BW2, 8*fc+BW1);

        [mag,phase,fpad] = basetran_mod(signal,fhop);

        [nfreqbins, nframes] = size(mag);
        stepfreq=fs/nfreqbins;
        f40=40/stepfreq;

        for k=1:f40

           m_F1(k,:) = (filtfilt(h_band1.Numerator, 1, mag(k,:)));
           m_F2(k,:) = (filtfilt(h_band2.Numerator, 1, mag(k,:)));
           m_F3(k,:) = (filtfilt(h_band3.Numerator, 1, mag(k,:)));
           m_F4(k,:) = (filtfilt(h_band4.Numerator, 1, mag(k,:)));
           m_F5(k,:) = (filtfilt(h_band5.Numerator, 1, mag(k,:)));
           m_F6(k,:) = (filtfilt(h_band6.Numerator, 1, mag(k,:)));
           m_F7(k,:) = (filtfilt(h_band7.Numerator, 1, mag(k,:)));
           m_F8(k,:) = (filtfilt(h_band8.Numerator, 1, mag(k,:)));

           mag_filt(k,:) = m_F1(k,:) + m_F2(k,:) + m_F3(k,:) + m_F4(k,:) + m_F5(k,:) + m_F6(k,:) + m_F7(k,:) + m_F8(k,:);

        end

        for i=1:f40
            for j=1:nframes
                if mag_filt(i,j) < 0
                    mag_filt(i,j)=0;
                end
            end
        end

        a=windowsize; 
        b=windowsize-f40;
        c=nframes;

        mag_filt(f40+1:a,:) = zeros(b,c);

        mag_filt_tot(1:windowsize,:) = [mag_filt(1:(windowsize/2+1),:); (mag_filt([(windowsize/2):-1:2],:))];

        enhanced = invbasetran_mod(mag_filt_tot,phase,fpad);
        enhanced = (enhanced-mean(enhanced))./std(enhanced);
  
   case 7
    
        h_band1 = bandpass_func(1*fc-BW1, 1*fc-BW2, 1*fc+BW2, 1*fc+BW1);
        h_band2 = bandpass_func(2*fc-BW1, 2*fc-BW2, 2*fc+BW2, 2*fc+BW1);
        h_band3 = bandpass_func(3*fc-BW1, 3*fc-BW2, 3*fc+BW2, 3*fc+BW1);
        h_band4 = bandpass_func(4*fc-BW1, 4*fc-BW2, 4*fc+BW2, 4*fc+BW1);
        h_band5 = bandpass_func(5*fc-BW1, 5*fc-BW2, 5*fc+BW2, 5*fc+BW1);
        h_band6 = bandpass_func(6*fc-BW1, 6*fc-BW2, 6*fc+BW2, 6*fc+BW1);
        h_band7 = bandpass_func(7*fc-BW1, 7*fc-BW2, 7*fc+BW2, 7*fc+BW1);

        [mag,phase,fpad] = basetran_mod(signal,fhop);

        [nfreqbins, nframes] = size(mag);
        stepfreq=fs/nfreqbins;
        f40=40/stepfreq;

        for k=1:f40

           m_F1(k,:) = (filtfilt(h_band1.Numerator, 1, mag(k,:)));
           m_F2(k,:) = (filtfilt(h_band2.Numerator, 1, mag(k,:)));
           m_F3(k,:) = (filtfilt(h_band3.Numerator, 1, mag(k,:)));
           m_F4(k,:) = (filtfilt(h_band4.Numerator, 1, mag(k,:)));
           m_F5(k,:) = (filtfilt(h_band5.Numerator, 1, mag(k,:)));
           m_F6(k,:) = (filtfilt(h_band6.Numerator, 1, mag(k,:)));
           m_F7(k,:) = (filtfilt(h_band7.Numerator, 1, mag(k,:)));

           mag_filt(k,:) = m_F1(k,:) + m_F2(k,:) + m_F3(k,:) + m_F4(k,:) + m_F5(k,:) + m_F6(k,:) + m_F7(k,:);

        end

        for i=1:f40
            for j=1:nframes
                if mag_filt(i,j) < 0
                    mag_filt(i,j)=0;
                end
            end
        end

        a=windowsize; 
        b=windowsize-f40;
        c=nframes;

        mag_filt(f40+1:a,:) = zeros(b,c);

        mag_filt_tot(1:windowsize,:) = [mag_filt(1:(windowsize/2+1),:); (mag_filt([(windowsize/2):-1:2],:))];

        enhanced = invbasetran_mod(mag_filt_tot,phase,fpad);
        enhanced = (enhanced-mean(enhanced))./std(enhanced);
   
   case 6
    
        h_band1 = bandpass_func(1*fc-BW1, 1*fc-BW2, 1*fc+BW2, 1*fc+BW1);
        h_band2 = bandpass_func(2*fc-BW1, 2*fc-BW2, 2*fc+BW2, 2*fc+BW1);
        h_band3 = bandpass_func(3*fc-BW1, 3*fc-BW2, 3*fc+BW2, 3*fc+BW1);
        h_band4 = bandpass_func(4*fc-BW1, 4*fc-BW2, 4*fc+BW2, 4*fc+BW1);
        h_band5 = bandpass_func(5*fc-BW1, 5*fc-BW2, 5*fc+BW2, 5*fc+BW1);
        h_band6 = bandpass_func(6*fc-BW1, 6*fc-BW2, 6*fc+BW2, 6*fc+BW1);

        [mag,phase,fpad] = basetran_mod(signal,fhop);

        [nfreqbins, nframes] = size(mag);
        stepfreq=fs/nfreqbins;
        f40=40/stepfreq;

        for k=1:f40

           m_F1(k,:) = (filtfilt(h_band1.Numerator, 1, mag(k,:)));
           m_F2(k,:) = (filtfilt(h_band2.Numerator, 1, mag(k,:)));
           m_F3(k,:) = (filtfilt(h_band3.Numerator, 1, mag(k,:)));
           m_F4(k,:) = (filtfilt(h_band4.Numerator, 1, mag(k,:)));
           m_F5(k,:) = (filtfilt(h_band5.Numerator, 1, mag(k,:)));
           m_F6(k,:) = (filtfilt(h_band6.Numerator, 1, mag(k,:)));

           mag_filt(k,:) = m_F1(k,:) + m_F2(k,:) + m_F3(k,:) + m_F4(k,:) + m_F5(k,:) + m_F6(k,:);

        end

        for i=1:f40
            for j=1:nframes
                if mag_filt(i,j) < 0
                    mag_filt(i,j)=0;
                end
            end
        end


        a=windowsize; 
        b=windowsize-f40;
        c=nframes;

        mag_filt(f40+1:a,:) = zeros(b,c);

        mag_filt_tot(1:windowsize,:) = [mag_filt(1:(windowsize/2+1),:); (mag_filt([(windowsize/2):-1:2],:))];

        enhanced = invbasetran_mod(mag_filt_tot,phase,fpad);
        enhanced = (enhanced-mean(enhanced))./std(enhanced);

   case 5
    
        h_band1 = bandpass_func(1*fc-BW1, 1*fc-BW2, 1*fc+BW2, 1*fc+BW1);
        h_band2 = bandpass_func(2*fc-BW1, 2*fc-BW2, 2*fc+BW2, 2*fc+BW1);
        h_band3 = bandpass_func(3*fc-BW1, 3*fc-BW2, 3*fc+BW2, 3*fc+BW1);
        h_band4 = bandpass_func(4*fc-BW1, 4*fc-BW2, 4*fc+BW2, 4*fc+BW1);
        h_band5 = bandpass_func(5*fc-BW1, 5*fc-BW2, 5*fc+BW2, 5*fc+BW1);

        [mag,phase,fpad] = basetran_mod(signal,fhop);

        [nfreqbins, nframes] = size(mag);
        stepfreq=fs/nfreqbins;
        f40=40/stepfreq;

        for k=1:f40

           m_F1(k,:) = (filtfilt(h_band1.Numerator, 1, mag(k,:)));
           m_F2(k,:) = (filtfilt(h_band2.Numerator, 1, mag(k,:)));
           m_F3(k,:) = (filtfilt(h_band3.Numerator, 1, mag(k,:)));
           m_F4(k,:) = (filtfilt(h_band4.Numerator, 1, mag(k,:)));
           m_F5(k,:) = (filtfilt(h_band5.Numerator, 1, mag(k,:)));

           mag_filt(k,:) = m_F1(k,:) + m_F2(k,:) + m_F3(k,:) + m_F4(k,:) + m_F5(k,:);

        end

        for i=1:f40
            for j=1:nframes
                if mag_filt(i,j) < 0
                    mag_filt(i,j)=0;
                end
            end
        end


        a=windowsize; 
        b=windowsize-f40;
        c=nframes;

        mag_filt(f40+1:a,:) = zeros(b,c);

        mag_filt_tot(1:windowsize,:) = [mag_filt(1:(windowsize/2+1),:); (mag_filt([(windowsize/2):-1:2],:))];

        enhanced = invbasetran_mod(mag_filt_tot,phase,fpad);
        enhanced = (enhanced-mean(enhanced))./std(enhanced);

end

end

