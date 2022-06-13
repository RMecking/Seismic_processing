clear all; close all; clc;
%----------------------------------------------
% author: Rebekka Mecking, in June, 2022
% Script to double the available traces in the dataset by inserting 
% traces in the centre between two adjacent traces
% Works solely in the time domain and uses cross-correlation
% to find the wavelet-shift between two adjacent traces in a moving time
% window
% Needs the SU_matlab toolbox or similar tools to load the seismic data
%-----------------------------------------------
disk='K';
Dir=[disk,':/SU_matlab'];
path(path,strcat(Dir,'/segy'));
path(path,strcat(Dir,'/seismic_plots'));
path(path,strcat(Dir,'/scaling_tapering'));
path(path,strcat(Dir,'/bp_filter'));
path(path,strcat(Dir,'/SegyMAT'));
addpath(Dir);

% input data can be segy or su-file
[D,H,S]=ReadSegy('S2R2_geom.sgy');
datacube=zeros(length(D(1,:))*2,1,length(D(:,1)));
for i=1:length(D(1,:))
    datacube(i*2,1,:)=D(:,i);
end
dt=H(1).dt/1e6;
ns=H(1).ns;
t=dt:dt:ns*dt;
dx=1; % Trace distance
ntr=72; % Traces/ Shot
% window length of mvong window calculated from minimum frequency 
% and phase shift between two traces
vmin=70; % minimum velocity in dataset
fmin=20; % minimum dtaa frequency
twinsamp=floor((1/fmin+dx/vmin)/dt); % window length for the correlation
gapsamp=ceil((0.05*twinsamp)); % window gap of moving time window
D2=zeros(length(t),1); % empty new dataset for interpolated data
outputname='S2R2_geom_interp.su';
S2=S; 
l=1;
ll=1;

for j=1:length(D(1,:))-1
    if mod(j,100)==0
        disp(['Interpolation is at trace ', num2str(j), ' of ', num2str(length(D(1,:))), '.']);
    end
%H2(l)=InitSegyTraceHeader(ns,dt*1e6);
H2(l)=make_empty_header;
% wenn tr1 die letzte spur ist, springe zur ersten Spur im nächsten schuss
if H(j).TraceNumber==ntr
    D2(:,l)=tr1;
    H2(l).fldr=H(j).FieldRecord;
    H2(l).tracf=ll;
    H2(l).tracl=l;
    H2(l).gelev=H(j).ReceiverGroupElevation;
    H2(l).selev=H(j).SourceSurfaceElevation;
    H2(l).gx=H(j).GroupX;
    H2(l).gy=H(j).GroupY;
    H2(l).sx=H(j).SourceX;
    H2(l).sy=H(j).SourceY;
    %H2(l).offset=H(j).offset;
    H2(l).dt=dt*1e9;
    H2(l).ns=ns;
    l=l+1;
    ll=1;
    continue;
end
H2(l+1)=make_empty_header;
% choose adjacent traces
tr1=D(:,j);
tr2=D(:,j+1);

% corrcorrelation in moving window
k=1;
A=zeros(ns,1);
for i=1:gapsamp:ns
    clearvars -except i gapsamp tr1 tr2 ns dt t twinsamp k A D2 l D H j H2 S2 ntr ll
   
   if i<=twinsamp/2
       dsub1=tr1(1:i+floor(twinsamp/2));
       dsub2=tr2(1:i+floor(twinsamp/2));
       t1=t(1:i+floor(twinsamp/2));
   elseif i>ns-twinsamp/2
       dsub1=tr1(i-floor(twinsamp/2):end);
       dsub2=tr2(i-floor(twinsamp/2):end);   
       t1=t(i-floor(twinsamp/2):end);
   else
       dsub1=tr1(i-floor(twinsamp/2)+1:i+floor(twinsamp/2));
       dsub2=tr2(i-floor(twinsamp/2)+1:i+floor(twinsamp/2));
       t1=t(i-floor(twinsamp/2)+1:i+floor(twinsamp/2));
   end
   
       dcorr=xcorr(dsub1,dsub2);
       [maxi,I]=max(dcorr);
       shiftsamp=I-length(dsub1); % sample shift at maximum correlation
       if shiftsamp<0 
       dsub1b=[zeros(abs(shiftsamp),1); dsub1];
       dsub2b=[dsub2; zeros(abs(shiftsamp),1)];
       else
       dsub2b=[zeros(abs(shiftsamp),1); dsub2];
       dsub1b=[dsub1; zeros(abs(shiftsamp),1)];           
       end
       dsub1b=dsub1b(ceil(abs(shiftsamp)/2)+1:length(dsub1b)-floor(abs(shiftsamp)/2));
       dsub2b=dsub2b(ceil(abs(shiftsamp)/2)+1:length(dsub2b)-floor(abs(shiftsamp)/2));
       t1(dsub2b==0 | dsub1b==0)=[];
       dsub1b(dsub2b==0)=[];
       dsub2b(dsub2b==0)=[];
       dsub2b(dsub1b==0)=[];
       dsub1b(dsub1b==0)=[];
       trinterp=0.5*(dsub1b+dsub2b);
       %trinterp(1:shiftsamp)=dsub2b(1:shiftsamp);
       %trinterp(length(dsub1)+1:end)=dsub1b(length(dsub1)+1:end);
%        plot(t1,dsub1b);
%        hold on;
%        plot(t1,dsub2b);
%        plot(t1,trinterp);
       % um halben shiftsamp-wert verschieben
       %figure(2)
       %plot(t1,dsub1);
       %hold on;
       %plot(t1,dsub2);
       %plot(t1,trinterp);
       i1=find(t1(1)==t);
       i2=find(t1(end)==t);
       if i<=twinsamp/2
           A(i1:i2,k)=trinterp;
           k=k+1;
       elseif i>=(ns-twinsamp/2)
           A(i1:i2,k)=trinterp;
           k=k+1;           
       else
           A(i1:i2,k)=trinterp;
           k=k+1;
       end
end
A(A==0)=NaN;
% Hier könnten noch Ausreißer über ein (1-2 Lambda)
% Standardabweichungskriterium
% entfernt werden...
C=mean(A,2,'omitnan'); % mittelwert über zeitsample aller berechneter fenster
D2(:,l)=tr1;


H2(l).fldr=H(j).FieldRecord;
H2(l).tracf=ll;
H2(l).tracl=l;
H2(l).gelev=H(j).ReceiverGroupElevation;
H2(l).selev=H(j).SourceSurfaceElevation;
H2(l).gx=H(j).GroupX;
H2(l).gy=H(j).GroupY;
H2(l).sx=H(j).SourceX;
H2(l).sy=H(j).SourceY;
%H2(l).offset=H(j).offset;
H2(l).dt=dt*1e9;
H2(l).ns=ns;
D2(:,l+1)=C;
H2(l).dt=dt;
H2(l+1).dt=dt*1e9;
H2(l).ns=ns;
H2(l+1).ns=ns;
H2(l+1)=H2(l);
H2(l+1).tracl=l+1;
H2(l+1).tracf=ll+1;
H2(l+1).fldr=H(j).FieldRecord;
H2(l+1).gelev=H(j).ReceiverGroupElevation+0.5*(H(j+1).ReceiverGroupElevation-H(j).ReceiverGroupElevation);
H2(l+1).selev=H(j).SourceSurfaceElevation+0.5*(H(j+1).SourceSurfaceElevation-H(j).SourceSurfaceElevation);
H2(l+1).gx=H(j).GroupX+0.5*(H(j+1).GroupX-H(j).GroupX);
H2(l+1).gy=H(j).GroupY+0.5*(H(j+1).GroupY-H(j).GroupY);
H2(l+1).sx=H(j).SourceX+0.5*(H(j+1).SourceX-H(j).SourceX);
H2(l+1).sy=H(j).SourceY+0.5*(H(j+1).SourceY-H(j).SourceY);
%H2(l+1).offset=sqrt((H(j).GroupX-H(j).SourceX).^2+(H(j).GroupY-H(j).SourceY).^2);
l=l+2;
ll=ll+2;
end

H2(l)=make_empty_header;
D2(:,l)=tr1;
    D2(:,l)=tr1;
    H2(l).fldr=H(j).FieldRecord;
    H2(l).tracf=ll;
    H2(l).tracl=l;
    H2(l).gelev=H(j).ReceiverGroupElevation;
    H2(l).selev=H(j).SourceSurfaceElevation;
    H2(l).gx=H(j).GroupX;
    H2(l).gy=H(j).GroupY;
    H2(l).sx=H(j).SourceX;
    H2(l).sy=H(j).SourceY;
    %H2(l).offset=H(j).offset;
    H2(l).dt=dt*1e9;
    H2(l).ns=ns;
D2(isnan(D2))=0;   
writesegy(outputname,D2,H2);