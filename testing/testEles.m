clear all global;

global BEAMLINE WF

load lcls_machine Initial
Initial.Momentum = 1.3;
Initial.sigz = 0.5e-3;
mybeam = MakeBeam6DGauss(Initial,1e3,4,1);
mybeam.Bunch.x(1,:) = mybeam.Bunch.x(1,:) + 50e-3;
% mybeam = CreateBlankBeam(1,1,1.3,1);
% mybeam.Bunch.x = ...
%    [0,-1e-4,1e-4,0,0,0,0,0,0,0,0;
%     0,0,0,-1e-4,1e-4,0,0,0,0,0,0;
%     0,0,0,0,0,-1e-4,1e-4,0,0,0,0;
%     0,0,0,0,0,0,0,-1e-4,1e-4,0,0;
%     0,0,0,0,0,0,0,0,0,-1e-4,1e-4;
%     1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3];
% mybeam.Bunch.Q = ones(1,size(mybeam.Bunch.x,2)) * mybeam.Bunch.Q;
% mybeam.Bunch.stop = zeros(1,size(mybeam.Bunch.x,2));

% BEAMLINE{1} = MultStruc(0.1,1e6,0,2,[0 0],1,'sext_normal');
% BEAMLINE{1} = MultStruc(0.1,1e6,pi/8,2,[0 0],1,'sext_skew');
% BEAMLINE{1} = QuadStruc(1,1,0,1,'quad');
% BEAMLINE{1} = QuadStruc(1,1,pi/4,1,'quad_skew');
% BEAMLINE{1} = CorrectorStruc(0.1,0.002,0,1,'xcor');
% BEAMLINE{1} = CorrectorStruc(0.1,0.002,0,2,'ycor');
% bfield = 0.2;
% bend = bfield * 299792458 / 1.3e9;
% BEAMLINE{1} = SBendStruc(1,bfield,bend,0,0,0,0,0,'sbend');
% BEAMLINE{1} = SBendStruc(1,bfield,bend,46.12e-3*ones(1,2),0,0,0,0,'sbend');
% BEAMLINE{1} = SBendStruc(1,bfield,bend,0,0.1*ones(1,2),0,0,0,'sbend');
BEAMLINE{1} = RFStruc(1.5, 10, 0, 1.3e3, 0,0,0,0,'acc_struc');
SetDesignMomentumProfile(1,1,1e9,1.3);
[stat,WF.TSR] = ParseSRWF( 'srwf_lcls_tran.dat', 0.01 );
[stat,WF.ZSR] = ParseSRWF( 'srwf_lcls_long.dat', 0.01 );

cavs = findcells(BEAMLINE,'Class','LCAV');
for cavnum=1:length(cavs)
  BEAMLINE{cavs(cavnum)}.Wakes = [1 1];
end

SetTrackFlags('SRWF_T',0,1,length(BEAMLINE));
SetTrackFlags('SRWF_Z',0,1,length(BEAMLINE));
[stat beamout] = TrackThru(1,length(BEAMLINE),mybeam,1,1);
if stat{1}~=1; error(stat{2}); end

SetTrackFlags('SRWF_T',1,1,length(BEAMLINE));
SetTrackFlags('SRWF_Z',1,1,length(BEAMLINE));
[stat wakebeamout] = TrackThru(1,length(BEAMLINE),mybeam,1,1);
if stat{1}~=1; error(stat{2}); end

figure(4)
subplot(121);hold on
plot(mybeam.Bunch.x(1,:)*1e3,mybeam.Bunch.x(2,:)*1e3,'rx')
plot(wakebeamout.Bunch.x(1,:)*1e3,wakebeamout.Bunch.x(2,:)*1e3,'go')
plot(beamout.Bunch.x(1,:)*1e3,beamout.Bunch.x(2,:)*1e3,'bo')
xlabel('x space / mm')
ylabel('x'' space / mrad')
% axis([-0.15 0.15 -0.15 0.15])
subplot(122);hold on
plot(mybeam.Bunch.x(3,:)*1e3,mybeam.Bunch.x(4,:)*1e3,'rx')
plot(wakebeamout.Bunch.x(3,:)*1e3,wakebeamout.Bunch.x(4,:)*1e3,'go')
plot(beamout.Bunch.x(3,:)*1e3,beamout.Bunch.x(4,:)*1e3,'bo')
xlabel('y space / mm')
ylabel('y'' space / mrad')
% axis([-0.15 0.15 -0.15 0.15])

figure(5);hold on
% plot(mybeam.Bunch.x(5,:),mybeam.Bunch.x(6,:),'rx')
plot(wakebeamout.Bunch.x(5,:),wakebeamout.Bunch.x(6,:),'go')
plot(beamout.Bunch.x(5,:),beamout.Bunch.x(6,:),'bo')
