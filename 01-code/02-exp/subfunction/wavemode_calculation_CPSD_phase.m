%i_method = 1: 时间谱
%i_method = 2: 阶比谱
function [GAMMA,freq]=wavemode_calculation_CPSD_phase(signal,fs,nk,mode,phase,rotor_speed,save_directory,i_file,fname,i_method)
    %[Pulse,rotor_speed]=keyRotation(data(:,end),fs);
    filename_label = ['method ', num2str(i_method)];
    mode=[-50:50];
    if i_method == 1
        data_fft = 2^floor(log2(length(signal)));
        data = signal(1:data_fft,1:nk);
        freq = [0:data_fft/2.56 - 1]*fs/data_fft;  %数据频域离散刻度
        data_freq = fft(data)*2/data_fft;
        data_freq = data_freq(1:data_fft/2.56,:); 
    for k=1:length(mode)
        GAMMA(:,k)=1/nk*data_freq(:,1:nk)*exp(mode(k)*i*phase).'; 
    end
    
    else   
    L_signal = length(signal);
    L_seg = round(L_signal/10);
    Wind = hamming(L_seg);
    Noverlap = round(L_seg/2);
    Nfft = 2^(ceil(log2(L_seg))+1);  
%     GAMMA = zeros(Nfft/2+1,nk+1);
    for k=1:nk
        for l = 1:nk
            [C{k,l},freq] = cpsd(signal(:,k),signal(:,l),Wind,Noverlap,Nfft,fs);          
        end
    end
    for m =1:length(mode)
        temp_f = zeros(Nfft/2+1,1);
        for k = 1:nk
            for l = 1:nk%2*pi*k/nk
                temp_f = temp_f + 0.5*C{k,l}*exp(i*mode(m)*phase(k))*exp(-i*mode(m)*phase(l));
            end
        end
        GAMMA(:,m) = temp_f/(nk*nk);
    end
    
    
%       GAMMA = abs(GAMMA);
    end
    GAMMA = 10*log10(abs(GAMMA)/4e-10);
    h=figure('Visible', 'on');
    set(gcf,'outerposition',get(0,'screensize'));%最大化
    colormap(jet);
    imagesc(mode,freq',GAMMA); 
     axis xy; %ylim([1,30]);
     testTime='试验-2019-12-13';
    title({[testTime,'-截止模态分析',' -CPSD method '];[char(fname(i_file)),'-转速: ',num2str(round(rotor_speed))]},'FontSize',14)
    saveas(h,[save_directory,'\','Image','-',strrep(char(fname(i_file)),'.wav','-'),'.fig'])
    saveas(h,[save_directory,'\','Image','-',strrep(char(fname(i_file)),'.wav','-'),'.png'])

 end
   
