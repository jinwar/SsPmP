clear

n = 0;
deptharray = 0:100:600;
distarray = 30:5:50;

for idep = 1:length(deptharray)
	for idist = 1:length(distarray)
		taup_com = ['taup_time -ph S -rayp -o taup_temp -h ',num2str(deptharray(idep)),' -deg ',num2str(distarray(idist))]; 
		system(taup_com);
		n = n+1;
		disp(n);
		rayp(idep,idist) = load('taup_temp');
	end
end
rayp = rayp./deg2km(1);

figure(145)
clf
contourf(distarray,deptharray,rayp);
xlabel('Epicenter Distance (deg)','fontsize',20);
ylabel('Depth (km)','fontsize',20);
colorbar
set(gca,'fontsize',20)
