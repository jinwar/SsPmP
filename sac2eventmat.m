% script to convert sac files to event mat for ppsms project
% written by Ge Jin

clear;

beforetime = 10;
aftertime = 20;
filt = [0.02 0.5];
Vp = 6.5;

event = dir('./data/2*');

for ie = 1:length(event)
	clear sac
	disp(event(ie).name);
	sacfiles = dir(['data/',event(ie).name,'/*.sac.z']);
	for i = 1:length(sacfiles)
		sac(i) = readsac(fullfile('data',event(ie).name,sacfiles(i).name));
	end
	save(fullfile('eventmat',[event(ie).name,'.mat']),'sac');
end
