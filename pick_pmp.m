function pick_S(event)
% script to pick sspmp
% written by Ge Jin


beforetime = 20;
aftertime = 30;
avgtime = [-3 3];
filt = [0.02 0.5];
Vp = 6.5;

load( ['eventmat/',event,'.mat'])

if ~isfield(sac,'isgood')
	for i=1:length(sac)
		sac(i).isgood = 1;
	end
end

stlas = [sac.STLA];
stlos = [sac.STLO];
[avgdist avgazi] = distance(mean(stlas),mean(stlos),sac(1).EVLA,sac(1).EVLO);
stnmR = {sacR.KSTNM};
%taup_com = ['taup_time -ph S -rayp -deg ',num2str(avgdist),' -h ',num2str(sac(1).EVDP),'> taup_temp'];
%system(taup_com);
%rayp = load('taup_temp');
%rayp = rayp/deg2km(1);

sort_key = [sac.STLA];
ids = 1:length(sac);

mat = [ids(:),sort_key(:)];

mat = sortrows(mat,2);

figure(46)
clf
hold on
offset = 0;
avgSS = 0;
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
		data = data./max(abs(data));
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
		ind = find(t>avgtime(1) & t < avgtime(2));
		avgSS = avgSS + data(ind);
	end
end
%title(['Dist: ',num2str(avgdist),' Azi: ',num2str(avgazi)],'fontsize',20);
title(event,'fontsize',20);
avgSS = avgSS(:)./offset;
avgSS = detrend(avgSS);
f_avgSS = fft(avgSS);
N = length(f_avgSS);
i = sqrt(-1);
avgSS = real(ifft([f_avgSS(1:round(N/2))*exp(-i*pi/2); f_avgSS(round(N/2)+1:end)*exp(+i*pi/2)]));

while 1
	[x,y,bot] = ginput(1);
	pick_id = round(y);
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
	figure(46)
	clf
	hold on
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
        lag = lag*sac(id).DELTA;
        lag = lag - beforetime -avgtime(1);
		data = data./max(abs(data));
		t = -beforetime:sac(id).DELTA:aftertime;
		syndt = sac(id).T1 - sac(id).T2;
		t3dt = sac(id).T3 - sac(id).T2;
		if length(t) > length(data)
			t = t(1:length(data));
		end
		offset = 1*i;
		t = lag;
		plot(t,data + offset);
		plot(syndt,offset,'rx','markersize',15);
		plot(t3dt,offset,'rv','markersize',15);
		data(find(data<0)) = 0;
		if sac(id).isgood
			area(t,data + offset,offset);
		end
		text(t(1)-3,offset,sac(id).KSTNM);
	end
%	title(['Dist: ',num2str(avgdist),' Azi: ',num2str(avgazi)],'fontsize',20);
	title(event,'fontsize',20);
	xlim([-beforetime aftertime])
	drawnow

end

save(['eventmat/',event,'.mat'],'sac','sacR')

