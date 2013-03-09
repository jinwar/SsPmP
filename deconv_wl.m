function decon = deconv_wl(data,wavelet,water_level,timeshift)
% function to calculate the decon using water level method

if ~exist('water_level','var')
	water_level = 0.05;
end
if ~exist('timeshift','var')
	timeshift = 0;
end

data = data(:).*tukeywin(length(data),0.25);
wavelet = wavelet(:).*tukeywin(length(data),0.25);

if length(data)~=length(wavelet)
	wavelet = wavelet(:).*tukeywin(length(wavelet),0.25);
	wavelet(end:length(data)) = 0;
end

fftdata = fft(data);
fftwavelet = fft(wavelet);

N = length(data);
faxis_data = [0:floor(N/2), -floor(N/2):-1].*1/N;
if length(faxis_data) > N
	faxis_data = faxis_data(1:N);
end

maxamp = max(abs(fftwavelet));

ind = find(abs(fftwavelet) < maxamp*water_level);
fftwavelet(ind) = fftwavelet(ind)./abs(fftwavelet(ind)).*maxamp.*water_level;

deconfft = fftdata(:)./fftwavelet(:);
i = sqrt(-1);
deconfft = deconfft.*exp(-i.*2*pi*faxis_data(:).*timeshift);

decon = real(ifft(deconfft));

end
