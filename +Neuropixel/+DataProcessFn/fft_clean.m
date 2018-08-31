function data = fft_clean(imec, data, varargin)
% from Will Allen
% data is passed in as nChannels x nSamples, and is processed internally as
% this transposed and then transposed on the way out

p = inputParser();
p.addParameter('thresh', 10, @isscalar);
p.addParameter('nbins', 20, @isscalar);
p.addParameter('showFigures', false, @islogical);
p.parse(varargin{:});

thresh = p.Results.thresh;
nbins = p.Results.nbins;
fDebug = p.Results.showFigures;

% data must be single
% data is passed in as nChannels x nSamples, and is processed internally as
% this transposed and then transposed on the way out
data = single(data)';

in = data;
if fDebug
    in = data;
end

nSkip_med = 4;
nw = 3; %frequency neighbors to set to zero

if thresh==0, thresh = []; end
if thresh == 0
    return
end

vrMu = mean(data, 1);
data = bsxfun(@minus, data, vrMu);

% if max(data(:)) > 1000
%     a = 1;
% end

n = size(data,1);
n_pow2 = 2^nextpow2(n);
if n < n_pow2
    data = fft(data, n_pow2);
else
    data = fft(data);
end
n1 = floor(n_pow2/2);
viFreq = (1:n1)';
% vrFft1 = abs(mean(bsxfun(@times, data(1+viFreq,:), viFreq), 2));
vrFft1 = (mean(bsxfun(@times, abs(data(1+viFreq,:)), viFreq), 2));
if fDebug
    vrFft0 = mean(abs(data(1+viFreq,:)), 2); vrFreq = (1:numel(vrFft1))/numel(vrFft1)*imec.fs/2;  
    figure; subplot(221); plot(vrFreq, 2*pow2db(vrFft0),'k.','MarkerSize',8); xlabel('Freq (Hz)'); ylabel('Power (dB)'); grid on; xlim([0 500]);     
end

n2 = round(n1/nbins); 
for ibin=1:nbins
    vi1 = (n2*(ibin-1) : n2*ibin) + 1;
    if ibin==nbins, vi1(vi1>n1)=[]; end
    vrFft2 = vrFft1(vi1);
%     vrFft2 = detrend(vrFft2);
    vrFft2 = vrFft2 - median(vrFft2(1:nSkip_med:end)); %mad transform
    vrFft1(vi1) = vrFft2 / median(abs(vrFft2(1:nSkip_med:end)));
%     vrFft1(vi1) = zscore((vrFft1(vi1)));    
end

% broaden spectrum
vl_noise = vrFft1>thresh;
vi_noise = find(vl_noise);
for i_nw=1:nw
    viA = vi_noise-i_nw;    viA(viA<1)=[];
    viB = vi_noise+i_nw;    viB(viB>n1)=[];
    vl_noise(viA)=1;
    vl_noise(viB)=1;
end
vi_noise = find(vl_noise);
if fDebug
    subplot(222); plot(vrFreq, vrFft1,'k.','MarkerSize',8); xlabel('Freq (Hz)'); ylabel('z-score (detrended)'); 
    grid on; xlim([0 500]); ylim([0 50]); set(gca,'YScale', 'linear');
    hold on; plot(get(gca,'XLim'),thresh*[1 1], 'r-','MarkerSize',8);
    hold on; plot(vrFreq(vi_noise), vrFft1(vi_noise), 'r.','MarkerSize',8);
    % hold on; plot(vrFreq(vi_noise), vrPowDb(vi_noise), 'r.');
end

data(1+vi_noise,:) = 0;
data(end-vi_noise+1,:) = 0;
data = real(ifft(data, n_pow2, 'symmetric')); %~30% faster than below
% data = real(ifft(data));
if n < n_pow2, data = data(1:n,:); end
data = bsxfun(@plus, data, vrMu); %add mean back

if fDebug
   % figure; set(gcf,'Color','w');
    subplot(223); plot(vrFreq, vrFft1,'k.','MarkerSize',8); xlabel('Freq (Hz)');
    axis([0 10000 -10 40]);
    n = min(10000, size(data,1));
    subplot(224); plot(in(1:n,1)); hold on; plot(data(1:n,1));
%     plot(data(1:n,1)-data(1:n,1));    
    
%     return;
%     sRateHz = 25000;
%     vrFreq = (1:size(data,1))/size(data,1)*sRateHz;
%     figure; plot(viFreq, vrFft0, '.');
    figure
    plot(viFreq, vrFft1, '.');
    hold on; plot(viFreq(vi_noise), vrFft1(vi_noise), 'o');
    grid on; ylim([-5 20]);
    
    plot(get(gca,'XLim'), thresh*[1,1]);
    axis([0 10000 -10 40]);
end

% transpose back to nChannels x nSamples
data = data';
