% script to convert sac files to event mat for ppsms project
% written by Ge Jin

clear;

beforetime = 10;
aftertime = 20;
filt = [0.02 0.5];
Vp = 6.5;

event = dir('./data/2*');

for ie = 1:length(event)
	clear sac sacR
	disp(event(ie).name);
	sacfiles = dir(['data/',event(ie).name,'/*.sac.BHZ']);
	for i = 1:length(sacfiles)
		sac(i) = readsac(fullfile('data',event(ie).name,sacfiles(i).name));
	end
	sacfiles = dir(['data/',event(ie).name,'/*.sac.BHR']);
	for i = 1:length(sacfiles)
		sacR(i) = readsac(fullfile('data',event(ie).name,sacfiles(i).name));
	end
	save(fullfile('eventmat',[event(ie).name,'.mat']),'sac','sacR');
end
