% script to plot sspmp
% written by Ge Jin

clear;

beforetime = 10;
aftertime = 20;
filt = [0.02 0.5];
Vp = 6.5;

event = '20102050535'

sacfiles = dir('data/',event,'/*.sac.z');
for i = 1:length(sacfiles)
	sac(i) = readsac(sacfiles(i).name);
end

stlas = [sac.STLA];
stlos = [sac.STLO];
avgdist = distance(mean(stlas),mean(stlos),sac(1).EVLA,sac(1).EVLO);
taup_com = ['taup_time -ph S -rayp -deg ',num2str(avgdist),' -h ',num2str(sac(1).EVDP),'> taup_temp'];
system(taup_com);
rayp = load('taup_temp');
rayp = rayp/deg2km(1);

sort_key = [sac.STLA];
ids = 1:length(sac);

mat = [ids(:),sort_key(:)];

mat = sortrows(mat,2);

figure(43)
clf
hold on
offset = 0;
for i=1:size(mat,1)
	id = mat(i,1);
	if ~isnan(sac(id).T2)
		taxis = sac(id).B:sac(id).DELTA:sac(id).B+sac(id).DELTA*(sac(id).NPTS-1);
		ind = find(taxis > sac(id).T2 - beforetime & taxis < sac(id).T2 + aftertime);
		data = sac(id).DATA1(ind);
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
		h = t./2./(Vp^(-2) - rayp^2);
%		offset = mat(i,2);
		offset = offset+1;
		plot(t,data + offset);
		plot(syndt,offset,'rx','markersize',15);
		data(find(data<0)) = 0;
		area(t,data + offset,offset);
		text(t(1)-3,offset,sac(id).KSTNM);
	end
end

figure(44)
clf
hold on
offset = 0;
for i=1:size(mat,1)
	id = mat(i,1);
	if ~isnan(sac(id).T2)
		taxis = sac(id).B:sac(id).DELTA:sac(id).B+sac(id).DELTA*(sac(id).NPTS-1);
		ind = find(taxis > sac(id).T2 - beforetime & taxis < sac(id).T2 + aftertime);
		data = sac(id).DATA1(ind);
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
		h = t./2./((Vp^(-2) - rayp^2).^.5);
%		offset = mat(i,2);
		offset = offset+1;
		plot(h,data + offset);
		data(find(data<0)) = 0;
		area(h,data + offset,offset);
		text(-25,offset,sac(id).KSTNM);
	end
end
xlim([-10 100])
