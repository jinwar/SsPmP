function pick_S(event)
% script to pick sspmp
% written by Ge Jin


beforetime = 10;
aftertime = 20;
avgtime = [-3 3];
filt = [0.02 0.5];
Vp = 6.5;
amp = 1.2;

load( ['eventmat/',event])

if ~isfield(sac,'isgood')
	for i=1:length(sac)
		sac(i).isgood = 1;
	end
end

stlas = [sac.STLA];
stlos = [sac.STLO];
[avgdist avgazi] = distance(mean(stlas),mean(stlos),sac(1).EVLA,sac(1).EVLO);
stnmR = {sacR.KSTNM};
taup_com = ['taup_time -ph S -rayp -deg ',num2str(avgdist),' -h ',num2str(sac(1).EVDP),'> taup_temp'];
system(taup_com);
rayp = load('taup_temp');
rayp = rayp/deg2km(1);
mt = [8.00, 4.6, 3.38 ]; % mt is the half space with transmitted waves
mi = [6.50, 3.75, 2.667]; % mi is the half space with incident wave
R = PSVRTmatrix(rayp,mi,mt);
Rpp = R(:,1); % P-to-P reflection
phaseshift = angle(Rpp);



sort_key = [sac.STLA];
ids = 1:length(sac);

mat = [ids(:),sort_key(:)];

mat = sortrows(mat,2);

figure(46)
clf
hold on
offset = 0;
for i=1:size(mat,1)
	id = mat(i,1);
	if isnan(sac(id).T2)
		sac(id).T2 = sac(id).T1;
	end
	if ~isnan(sac(id).T2)
		taxis = sac(id).B:sac(id).DELTA:sac(id).B+sac(id).DELTA*(sac(id).NPTS-1);
		ind = find(taxis > sac(id).T2 - beforetime & taxis < sac(id).T2 + aftertime);
		data = sac(id).DATA1(ind);
        if isempty(data)
            continue;
        end
		fN = 1/sac(id).DELTA/2;
		[b,a] = butter(2,[filt(1)/fN filt(2)/fN]);
		data = filtfilt(b,a,data);
		data = data./max(abs(data))*amp;
%		data = data.*(max(mat(:,2))-min(mat(:,2)))/20;
		t = -beforetime:sac(id).DELTA:aftertime;
		syndt = sac(id).T1 - sac(id).T2;
		if length(t) > length(data)
			t = t(1:length(data));
		end
%		h = t./2./(Vp^(-2) - rayp^2);
%		offset = mat(i,2);
		offset = i;
		plot(t,data + offset);
		plot(syndt,offset,'rx','markersize',15);
		data(find(data<0)) = 0;
		area(t,data + offset,offset);
		text(t(1)-3,offset,sac(id).KSTNM);
	end
end
%title(['Dist: ',num2str(avgdist),' Azi: ',num2str(avgazi)],'fontsize',20);
title(event,'fontsize',20);

isplotorigin = 1;
while 1
	[x,y,bot] = ginput(1);
	pick_id = round(y/2);
	if bot == 's'
		sac(mat(pick_id,1)).T2 = sac(mat(pick_id,1)).T2 + x;
	end
	if bot == 'q'
		break;
	end
	if bot == 'd'
		sac(mat(pick_id,1)).isgood = ~sac(mat(pick_id,1)).isgood;
	end
	if bot == 'r'
		sac(mat(pick_id,1)).T3 = sac(mat(pick_id,1)).T2 + x;
	end
	if bot == 'x'  % plot original xcor
		isplotorigin = 0;
	end
	if bot == 'c'  %plot original waveform
		isplotorigin = 1;
	end
	if bot == 'z'  % plot phase shift xcor
		isplotorigin = -1;
    end
    if bot == 'm' 
        save(['eventmat/',event],'sac','sacR')
		plot_S(event);
	end
	figure(46)
	clf
	hold on
	if isplotorigin == 1
		avgSS = 0; goodsta = 0; avgSSR = 0; goodstaR = 0;
		for i=1:size(mat,1)
			id = mat(i,1);
			taxis = sac(id).B:sac(id).DELTA:sac(id).B+sac(id).DELTA*(sac(id).NPTS-1);
			ind = find(taxis > sac(id).T2 - beforetime & taxis < sac(id).T2 + aftertime);
			data = sac(id).DATA1(ind);
			if isempty(data)
				continue;
			end
			fN = 1/sac(id).DELTA/2;
			[b,a] = butter(2,[filt(1)/fN filt(2)/fN]);
			data = filtfilt(b,a,data);
			data = data./max(abs(data))*amp;
			idR = find(ismember(stnmR,sac(id).KSTNM));
			t = -beforetime:sac(id).DELTA:aftertime;
			syndt = sac(id).T1 - sac(id).T2;
			t3dt = sac(id).T3 - sac(id).T2;
			if length(t) > length(data)
				t = t(1:length(data));
			end
			offset = 2*i;
			plot(t,data + offset);
			plot(syndt,offset,'rx','markersize',15);
			data(find(data<0)) = 0;
			if sac(id).isgood
				area(t,data + offset,offset);
				plot(t3dt,offset,'rv','markersize',15);
			end
			if ~isempty(idR)
				dataR = sacR(idR).DATA1(ind);
				dataR = filtfilt(b,a,dataR);
				dataR = dataR./max(abs(dataR))*amp;
				plot(t,dataR + offset,'r');
				dataR(find(dataR>0)) = 0;
				if sac(id).isgood
					area(t,dataR + offset,offset,'facecolor','r');
				end
			end
			text(t(1)-3,offset,sac(id).KSTNM);
			if sac(id).isgood
				ind = find(t>avgtime(1) & t < avgtime(2));
				avgSS = avgSS + data(ind);
                if ~isempty(idR)
                    avgSSR = avgSSR + dataR(ind);
                    goodstaR = goodstaR+1;
                end
				goodsta = goodsta+1;
			end
		end
	%	title(['Dist: ',num2str(avgdist),' Azi: ',num2str(avgazi)],'fontsize',20);
		title(event,'fontsize',20);
		avgSS = avgSS(:)./goodsta;
		avgSS = detrend(avgSS);
		f_avgSS = fft(avgSS);
		N = length(f_avgSS);
		i = sqrt(-1);
%		phaseshift = pi/2;
		ps_avgSS = real(ifft([f_avgSS(1:round(N/2))*exp(-i*phaseshift); f_avgSS(round(N/2)+1:end)*exp(+i*phaseshift)]));

		avgSSR = avgSSR(:)./goodstaR;
		avgSSR = detrend(avgSSR);
		f_avgSSR = fft(avgSSR);
		N = length(f_avgSSR);
		i = sqrt(-1);
%		phaseshift = pi/2;
		ps_avgSSR = real(ifft([f_avgSSR(1:round(N/2))*exp(+i*phaseshift); f_avgSSR(round(N/2)+1:end)*exp(-i*phaseshift)]));
		
	elseif isplotorigin == -1
		for i=1:size(mat,1)
			id = mat(i,1);
			taxis = sac(id).B:sac(id).DELTA:sac(id).B+sac(id).DELTA*(sac(id).NPTS-1);
			ind = find(taxis > sac(id).T2 - beforetime & taxis < sac(id).T2 + aftertime);
			data = sac(id).DATA1(ind);
			if isempty(data)
				continue;
			end
			fN = 1/sac(id).DELTA/2;
			[b,a] = butter(2,[filt(1)/fN filt(2)/fN]);
			data = filtfilt(b,a,data);
			[data lag] = xcorr(data,ps_avgSS);
            data = data./max(abs(data))*amp;
			lag = lag*sac(id).DELTA;
			lag = lag - beforetime -avgtime(1);
			t = -beforetime:sac(id).DELTA:aftertime;
			syndt = sac(id).T1 - sac(id).T2;
			t3dt = sac(id).T3 - sac(id).T2;
			if length(t) > length(data)
				t = t(1:length(data));
			end
			offset = 2*i;
			t = lag;
			plot(t,data + offset);
            idR = find(ismember(stnmR,sac(id).KSTNM));
			if ~isempty(idR)
                dataR = sacR(idR).DATA1(ind);
                dataR = filtfilt(b,a,dataR);
                dataR = xcorr(dataR,ps_avgSSR);
                dataR = dataR./max(abs(dataR))*amp;
                plot(t,dataR + offset,'r');
            end
			plot(syndt,offset,'rx','markersize',15);
			data(find(data<0)) = 0;
			if sac(id).isgood
				area(t,data + offset,offset);
				plot(t3dt,offset,'rv','markersize',15);
			end
			text(t(1)-3,offset,sac(id).KSTNM);
		end
		title(event,'fontsize',20);
	elseif isplotorigin == 0;
		for i=1:size(mat,1)
			id = mat(i,1);
			taxis = sac(id).B:sac(id).DELTA:sac(id).B+sac(id).DELTA*(sac(id).NPTS-1);
			ind = find(taxis > sac(id).T2 - beforetime & taxis < sac(id).T2 + aftertime);
			data = sac(id).DATA1(ind);
			if isempty(data)
				continue;
			end
			fN = 1/sac(id).DELTA/2;
			[b,a] = butter(2,[filt(1)/fN filt(2)/fN]);
			data = filtfilt(b,a,data);
			[data lag] = xcorr(data,avgSS);
            data = data./max(abs(data))*amp;
			lag = lag*sac(id).DELTA;
			lag = lag - beforetime -avgtime(1);
			t = -beforetime:sac(id).DELTA:aftertime;
			syndt = sac(id).T1 - sac(id).T2;
			t3dt = sac(id).T3 - sac(id).T2;
			if length(t) > length(data)
				t = t(1:length(data));
			end
			offset = 2*i;
			t = lag;
			plot(t,data + offset);
            idR = find(ismember(stnmR,sac(id).KSTNM));
			if ~isempty(idR)
                dataR = sacR(idR).DATA1(ind);
                dataR = filtfilt(b,a,dataR);
                dataR = xcorr(dataR,avgSSR);
                dataR = dataR./max(abs(dataR))*amp;
                plot(t,dataR + offset,'r');
            end
			plot(syndt,offset,'rx','markersize',15);
			data(find(data<0)) = 0;
			if sac(id).isgood
				area(t,data + offset,offset);
				plot(t3dt,offset,'rv','markersize',15);
			end
			text(t(1)-3,offset,sac(id).KSTNM);
		end
		title(event,'fontsize',20);
    end
    plot([0 0],[0 offset],'r--');
	xlim([-beforetime aftertime])
    ylim([0 offset+2])
	drawnow
	

end

save(['eventmat/',event],'sac','sacR')

