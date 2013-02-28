#!/bin/csh

foreach event (`ls -d data/20*`)
cd $event
sac<<!
r *.sac.?
ch o 0
wh
q
!
taup_setsac -mod prem -ph S-1 -evdpkm *.sac.z
cd ../..
end
