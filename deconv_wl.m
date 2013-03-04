function decon = deconv_wl(data,wavelet)
% function to calculate the decon using water level method

water_level = 0.05;
data = data(:);
wavelet = wavelet(:);

fftdata = fft(data);
fftwavelet = fft(wavelet);

N = length(data);
faxis_data = [0:floor(N/2), -floor(N/2):-1].*1/N;
if length(faxis_data) > N
	faxis_data = faxis_data(1:N);
end

N = length(wavelet);
faxis_wavelet = [0:floor(N/2), -floor(N/2):-1].*1/N;
if length(faxis_wavelet) > N
	faxis_wavelet = faxis_wavelet(1:N);
end

res_fftwavelet = interp1(faxis_wavelet,fftwavelet,faxis_data,'spline','extrap');
maxamp = max(abs(res_fftwavelet));

ind = find(abs(res_fftwavelet) < maxamp*water_level);
res_fftwavelet(ind) = res_fftwavelet(ind)./abs(res_fftwavelet(ind)).*maxamp.*water_level;

deconfft = fftdata(:)./res_fftwavelet(:);

decon = real(ifft(deconfft));

end
