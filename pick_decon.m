function pick_decon(event)
% script to pick sspmp
% written by Ge Jin

isdebug = 1;

beforetime = 5;
aftertime = 15;
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

n = 0;
for i = 5:5:85
	for j = -5:-5:-85
		n = n+1;
		search_angle(n,1) = i;
		search_angle(n,2) = j;
	end
end

% do the rotate and decon
for ista = 1:length(sac)
	idR = find(ismember(stnmR,sac(ista).KSTNM));
	if isempty(idR)
		continue;
	end
	taxis = sac(ista).B:sac(ista).DELTA:sac(ista).B+sac(ista).DELTA*(sac(ista).NPTS-1);
	ind = find(taxis > sac(ista).T2 - beforetime & taxis < sac(ista).T2 + aftertime);
	data = sac(ista).DATA1(ind);
	fN = 1/sac(ista).DELTA/2;
	[b,a] = butter(2,[filt(1)/fN filt(2)/fN]);
	data = filtfilt(b,a,data);
	dataR = sacR(idR).DATA1(ind);
	dataR = filtfilt(b,a,dataR);
	clear err
	for iangle = 1:size(search_angle(:,1))
		dist1 = abs(dataR.*sind(search_angle(iangle,1)) - data.*cosd(search_angle(iangle,1)));
		dist2 = abs(dataR.*sind(search_angle(iangle,2)) - data.*cosd(search_angle(iangle,2)));
		dist = min([dist1(:),dist2(:)],[],2);
		err(iangle) = sum(dist.^2);
	end
	[temp best_angle_ind] = min(err);
	theta1 = search_angle(best_angle_ind,1);
	theta2 = search_angle(best_angle_ind,2);
	sac(ista).dataP = data.*sind(theta1) + dataR.*cosd(theta1);
	dataP = data.*sind(theta1) + dataR.*cosd(theta1);
	sac(ista).dataS = data.*sind(theta2) + dataR.*cosd(theta2);
	dataS = data.*sind(theta2) + dataR.*cosd(theta2);
	sac(ista).theta1 = theta1;
	sac(ista).theta2 = theta2;
%	sac(ista).decon = deconv_wl(sac(ista).dataS,sac(ista).dataP,0.01,beforetime/sac(ista).DELTA);
	sac(ista).decon = deconv_wl(dataP,dataS,0.001,beforetime/sac(ista).DELTA);

	if isdebug
		figure(88)
		clf
		hold on
		plot(dataR,data);
		x = [min(dataR) max(dataR)];
		plot(x,[x(1)*tand(search_angle(best_angle_ind,1)) x(2)*tand(search_angle(best_angle_ind,1))],'r--');
		plot(x,[x(1)*tand(search_angle(best_angle_ind,2)) x(2)*tand(search_angle(best_angle_ind,2))],'r--');
		figure(89)
		clf
		subplot(3,1,1)
		plot(sac(ista).dataP);
		title('P')
		subplot(3,1,2)
		plot(sac(ista).dataS);
		title('S')
		subplot(3,1,3)
		plot(sac(ista).decon);
		title('P/S')
		pause
	end
end



figure(55)
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
	end
end
%title(['Dist: ',num2str(avgdist),' Azi: ',num2str(avgazi)],'fontsize',20);
title(event,'fontsize',20);

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
	figure(55)
	clf
	hold on
	for i=1:size(mat,1)
		id = mat(i,1);
		if isempty(sac(id).decon)
			continue;
		end
		taxis = sac(id).B:sac(id).DELTA:sac(id).B+sac(id).DELTA*(sac(id).NPTS-1);
		ind = find(taxis > sac(id).T2 - beforetime & taxis < sac(id).T2 + aftertime);
		data = sac(id).DATA1(ind);
		fN = 1/sac(id).DELTA/2;
		[b,a] = butter(2,[filt(1)/fN filt(2)/fN]);
%		data = filtfilt(b,a,data);
		data = sac(id).decon;
		data = data./max(abs(data));
		t = -beforetime:sac(id).DELTA:aftertime;
		syndt = sac(id).T1 - sac(id).T2;
		t3dt = sac(id).T3 - sac(id).T2;
		if length(t) > length(data)
			t = t(1:length(data));
		end
		offset = 2*i;
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
	drawnow

end

save(['eventmat/',event,'.mat'],'sac','sacR')

