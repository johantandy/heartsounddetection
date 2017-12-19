function gtfeatures = gtfeatures(ctl_list)
% Generate gammatone features (GF) and gammatone frequency cepstral
% coefficients (GFCC).

% ctl_list: list of files to be processed
% sampFreq: sampling frequency, default is 8000
% bWav: boolean variable indicating whether speech data is stored in the
%       WAV format
% numChannel: number of channels

% Written by Yang Shao, and adapted by Xiaojia Zhao in Oct'11

% if ~exist('ctl_list', 'var')
%     fprintf(1, 'Input file list is missing!\n');
%     exit;
% end

if ~exist('sampFreq', 'var')
    sampFreq = 16000;
end

if ~exist('bWav', 'var')
    bWav = 1;
end

if ~exist('numChannel', 'var')
    numChannel = 64;
end


gt = gen_gammaton(sampFreq, numChannel);  % get gammatone filterbank


    


    if bWav == 1
        rawfd = fopen(ctl_list);     % Training files are stored in the WAV format
        [sig, count] = fread(rawfd, inf, 'int16');
        sig = sig(23:end);
        fclose(rawfd);
    else
        sig = dlmread(entry);  % Test files are stored in the ASCII format
    end
    
    
    sig = reshape(sig, 1, length(sig));  % gammatone filtering function requires input to be row vector
    
    %fprintf('file read'); fprintf('\n');
    
    g=fgammaton(sig, gt, sampFreq, numChannel);     % gammatone filter pass and decimation to get GF features
    %fprintf('gammatone filtered'); fprintf('\n');
    
    gtfeatures = gtf2gtfcc(g, 2, 23);  % apply dct to get GFCC features with 0th coefficient removed

%     writeHTK(strcat(entry,'.gtf'), g);     % write GF features
% 
%     writeHTK(strcat(entry,'.gtfcc'), gfcc);     % write GFCC features
%     
%     fprintf('output written'); fprintf('\n');
   % gtcc = loadData(entry,'gtfcc');
    %gtcc=reshape(gtcc,1,[]);
    
        dlmwrite('feature.csv',gtfeatures,'delimiter',',');
   
    clear g sig
    %fprintf('time consumed: %f (sec)',toc); fprintf('\n');
      
end

%files = dir('C:\Users\Lenovo\Documents\Johan\data\1\*.wav');
%for file = files'
%    disp(file.name)
    % Do some stuff
%end
