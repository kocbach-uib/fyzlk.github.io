%  Kepler trajectory; if parameter a>0, deviation from Kepler
%  Particle starts at perihelium; Verlet algorithm is used
%  Ladislav Kocbach and Suhail Lubbad, 2001-2009
r0y=0; r0x=input('Enter distance from center   > ');
vfact=input('Enter kinetic energy factor  > ');
a=input('a=0 is Coulomb;   > ');
%velocity
v0x=0; v0y=sqrt(1/r0x*vfact);
ntime=6000;
tmax=270;     t=0:tmax/ntime:tmax;
m=1;  % strength
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v0=[ v0x;    v0y];
r0=[ r0x;    r0y];
r1=r0+v0*tmax/ntime;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=zeros(2,ntime);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  start
     r(:,1)= r0;
     r(:,2)= r1 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ccount=0;      hold on;
pplo=plot(r(1,1:2),r(2,1:2),'r-'); axis equal; hold on;
set(gca,'xlim',[-r0x-0.3 r0x+0.3],'ylim',[-r0x-0.3 r0x+0.3])
for nt = 3:ntime 
    Rt=r(:,nt-1);  force=- m*Rt/sqrt(a^2+sum(Rt.^2))^3;
    r(:,nt)= 2*r(:,nt-1)- r(:,nt-2) +force *( t(nt-1)-t(nt-2) )^2;  
	if(ccount==30)  
	  set(pplo,'xdata',(r(1,1:nt)),'ydata',(r(2,1:nt)));
	  drawnow;ccount=0; 
	end
	ccount=ccount+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pplb=plot(r(1,1),r(2,1),'o');        
for kt=1:50:ntime 
   set(pplb,'xdata',(r(1,kt)),'ydata',(r(2,kt))); drawnow; 
end
%  DATA FOR kepler: ntime=6000; tmax=270;   5    0.253    0.9
%  to get a nice rosett





