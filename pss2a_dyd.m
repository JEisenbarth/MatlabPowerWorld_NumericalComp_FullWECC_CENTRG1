function [flag] = pss2a_dyd(fileID,n,name,kv,id,J1,K1,J2,K2,Tw1,Tw2,Tw3,...
    Tw4,T6,T7,Ks2,Ks3,Ks4,T8,T9,N,m,Ks1,T1,T2,T3,T4,Vstmax,Vstmin,a,Ta,Tb)
%pss2a_dyd This function will write the PSLF pss2a model to a dyd file.
%   fileID must be created before calling function.
%   All input arguments must be strings besides fileID.
%   flag=0 if complete, flag=1 if not complete 
%   Parameters:
%   EPCL Default Variable Data Description
%   J1 1.0 Input signal #1 code
%   K1 0.0 Input signal #1 remote bus number
%   J2 3.0 Input signal #2 code
%    K2 0.0 Input signal #2 remote bus number
%   Tw1 2.0 First washout on signal #1, sec.
%    Tw2 2.0 Second washout on signal #1, sec.
%    Tw3 2.0 First washout on signal #2, sec.
%    Tw4 0.0 Second washout on signal #2, sec.
%   T6 0.0 Time constant on signal #1, sec.
%    T7 2.0 Time constant on signal #2, sec.
%    Ks2 0.2 Gain on signal #2
%    Ks3 1.0 Gain on signal #2
%    Ks4 1.0 Gain on signal #2
%    T8 0.5 Lead of ramp tracking filter
%    T9 0.1 Lag of ramp tracking filter
%   n 1.0 Order of ramp tracking filter
%   m 5.0 Order of ramp tracking filter
%   Ks1 10.0 Stabilizer gain
%   T1 0.25 Lead/lag time constant, sec.
%   T2 0.04 Lead/lag time constant, sec.
%   T3 0.2 Lead/lag time constant, sec.
%   T4 0.03 Lead/lag time constant, sec.
%   Vstmax 0.1 Stabilizer output max limit, p.u.
%   Vstmin -0.1 Stabilizer output min limit, p.u.
%   a 1. Lead/lag num. Gain. (not in IEEE model)
%   Ta 0. Lead/lag time constant, sec. (not in IEEE model)
%   Tb 0. Lead/lag time constant, sec. (not in IEEE model)
flag=1;
fprintf(fileID,'pss2a     %s "%s" %s "%s" : #9 %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n',...
    n,name,kv,id,J1,K1,J2,K2,Tw1,Tw2,Tw3,Tw4,T6,T7,Ks2,Ks3,Ks4,T8,T9,...
    N,m,Ks1,T1,T2,T3,T4,Vstmax,Vstmin,a,Ta,Tb);
flag=0;
end

