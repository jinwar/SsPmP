% script to plot sspmp
% written by Ge Jin

load seiscmap
lalim=[-11.2 -7.8];
lolim=[148.8 152.5];

hrange = [20 36];
seiscmap = seiscmap(5:end,:);
hx = linspace(hrange(1),hrange(2),size(seiscmap,1));

figure(45)
clf
ax = worldmap(lalim, lolim);
set(ax, 'Visible', 'off')
load pngcoastline
geoshow([S.Lat], [S.Lon], 'Color', 'black','linewidth',1)

eventmatfiles = dir('eventmat/*.mat');
for ie = 1:length(eventmatfiles)
	clear sac
	load( ['eventmat/',eventmatfiles(ie).name])

	if ~isfield(sac,'isgood')
		continue;
	end
	disp(eventmatfiles(ie).name);
	beforetime = 10;
	aftertime = 20;
	filt = [0.02 0.5];
	Vp = 6.5;

	stlas = [sac.STLA];
	stlos = [sac.STLO];
	evla = sac(1).EVLA;
	evlo = sac(1).EVLO;
	[avgdist avgazi] = distance(mean(stlas),mean(stlos),sac(1).EVLA,sac(1).EVLO);
	taup_com = ['taup_time -ph S -rayp -deg ',num2str(avgdist),' -h ',num2str(sac(1).EVDP),'> taup_temp'];
	system(taup_com);
	rayp = load('taup_temp');
	rayp = rayp/deg2km(1);

	[dists azis] = distance(stlas,stlos,evla,evlo);


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
		plotm(plat,plon,'ro','markerfacecolor',pointcolor,'markersize',10);
		textm(plat,plon+0.05,sac(goodind(i)).KSTNM);
	end
	colorbar
	colormap(seiscmap)
	caxis(hrange);

end
