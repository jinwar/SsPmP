function plot_S(event)
% script to plot sspmp
% written by Ge Jin

beforetime = 10;
aftertime = 20;
filt = [0.02 0.5];
Vp = 6.5;

load( ['eventmat/',event,'.mat'])
load seiscmap

stlas = [sac.STLA];
stlos = [sac.STLO];
evla = sac(1).EVLA;
evlo = sac(1).EVLO;
[avgdist avgazi] = distance(mean(stlas),mean(stlos),sac(1).EVLA,sac(1).EVLO);
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
	if (sac(id).isgood)
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
		t3dt = sac(id).T3 - sac(id).T2;
		if length(t) > length(data)
			t = t(1:length(data));
		end
		h = t./2./(Vp^(-2) - rayp^2);
%		offset = mat(i,2);
		offset = offset+1;
		plot(t,data + offset);
		plot(syndt,offset,'rx','markersize',15);
		plot(t3dt,offset,'rv','markersize',15);
		data(find(data<0)) = 0;
		area(t,data + offset,offset);
		text(t(1)-3,offset,sac(id).KSTNM);
	end
end
%title(['Dist: ',num2str(avgdist),' Azi: ',num2str(avgazi)],'fontsize',20);
	title(event,'fontsize',20);

figure(44)
clf
hold on
offset = 0;
for i=1:size(mat,1)
	id = mat(i,1);
	if sac(id).isgood
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
		t3dt = sac(id).T3 - sac(id).T2;
		if length(t) > length(data)
			t = t(1:length(data));
		end
		h = t./2./((Vp^(-2) - rayp^2).^.5);
%		offset = mat(i,2);
		offset = offset+1;
		plot(h,data + offset);
		t3h = t3dt./2./((Vp^(-2) - rayp^2).^.5);
		depth(offset) = t3h;
		plot(t3h,offset,'rv','markersize',15);
		data(find(data<0)) = 0;
		area(h,data + offset,offset);
		text(-25,offset,sac(id).KSTNM);
	end
end
xlim([-10 100])
%title(['Dist: ',num2str(avgdist),' Azi: ',num2str(avgazi)],'fontsize',20);
title(event,'fontsize',20);

% plot the map distribution of Moho Thickness
lalim=[-11.2 -7.8];
lolim=[148.8 152.5];
[dists azis] = distance(stlas,stlos,evla,evlo);
hrange = [20 36];
seiscmap = seiscmap(5:end,:);
hx = linspace(hrange(1),hrange(2),size(seiscmap,1));

figure(45)
clf
	ax = worldmap(lalim, lolim);
	set(ax, 'Visible', 'off')
	load pngcoastline
	geoshow([S.Lat], [S.Lon], 'Color', 'black','linewidth',1)
	goodind = find([sac.isgood]);
	plotm(stlas(goodind),stlos(goodind),'bv');
	for i =1:length(goodind)
		id = goodind(i);
		pdist = km2deg(tan(asin(rayp*Vp))*30);
		[plat plon] = reckon(stlas(goodind(i)),stlos(goodind(i)),pdist,azis(goodind(i)));
		t3dt = sac(id).T3 - sac(id).T2;
		if isnan(t3dt)
			continue;
		end
		t3h = t3dt./2./((Vp^(-2) - rayp^2).^.5);
		pointcolor = interp1(hx,seiscmap,t3h,'nearest','extrap');
		plotm(plat,plon,'ro','markerfacecolor',pointcolor,'markersize',20);
		textm(plat,plon+0.05,sac(goodind(i)).KSTNM);
	end
	colorbar
	colormap(seiscmap)
	caxis(hrange);
end
