function reset_S(event)
% script to pick sspmp
% written by Ge Jin

beforetime = 20;
aftertime = 40;
filt = [0.02 0.5];
Vp = 6.5;

load( ['eventmat/',event,'.mat'])

for i=1:length(sac)
	sac(i).isgood = 1;
	sac(i).T2 = sac(i).T1;
    sac(i).T3 = NaN;
end

save(['eventmat/',event,'.mat'],'sac','sacR')

