close all
clear all
clc
%set(0,'defaulttextinterpreter','none','defaultaxesposition',[0.10    0.10    .89    .8]);
figure('Position',[1531         203         707         400]);

%   For 150mm 1Hz In plane, use lambda3.66, out of plane use 3.77
load RunData_000000.mat
%clear all
%load Output.lambda3.77.mat





% BarltropMx1Hz150mm
% scatter(data_srt_int(:,1),smooth(data_srt_int(:,2)),'.'); hold all;
% plot(Times-4.0254,-1.*(1.5*Mx + 0.9*SelfMoment)-0.,'-r','linewidth',2)

% BodySurface2
% BarltropMy1Hz150mm
% scatter(data_srt_int(:,1),smooth(data_srt_int(:,2)),'.'); hold all;
% plot(Times-3.95,-2*Mout_of_plane,'-r','linewidth',2)
% axis([0 4    0.2    1.8])


C1 = [BodyPointsX(:,1) BodyPointsY(:,1) BodyPointsZ(:,1)];
C2 = [BodyPointsX(:,2) BodyPointsY(:,2) BodyPointsZ(:,2)];
C3 = [BodyPointsX(:,3) BodyPointsY(:,3) BodyPointsZ(:,3)];
C4 = [BodyPointsX(:,4) BodyPointsY(:,4) BodyPointsZ(:,4)];
Norms = cross(C3-C1,C4-C2);
Norms = [Norms(:,1)./sqrt(dot(Norms,Norms,2)) Norms(:,2)./sqrt(dot(Norms,Norms,2)) Norms(:,3)./sqrt(dot(Norms,Norms,2))];
CollocPts = 0.25*(C1 + C2 + C3 + C4);


Rads = Rloc;
Chrd = Cloc;
R0 = Rads(BodySurface0);
C0 = Chrd(BodySurface0);


Xcp0 = CollocPts_x(BodySurface0);
Ycp0 = CollocPts_y(BodySurface0);
Zcp0 = CollocPts_z(BodySurface0);
PressChord = zeros(size(C0(1,:)));

clc
%%  NACA 63-215
Foil.N63215.LS.x = [0,0.00601000000000000,0.00863000000000000,0.0138000000000000,0.0265200000000000,0.0517100000000000,0.0767700000000000,0.101770000000000,0.151660000000000,0.201480000000000,0.251250000000000,0.301000000000000,0.350740000000000,0.400480000000000,0.450230000000000,0.500000000000000,0.549810000000000,0.599650000000000,0.649530000000000,0.699470000000000,0.749450000000000,0.799490000000000,0.849570000000000,0.899700000000000,0.949860000000000,1;];
Foil.N63215.LS.z = [0,-0.0115000000000000,-0.0138800000000000,-0.0176600000000000,-0.0242000000000000,-0.0332800000000000,-0.0399900000000000,-0.0453500000000000,-0.0533600000000000,-0.0589500000000000,-0.0625900000000000,-0.0644800000000000,-0.0647000000000000,-0.0631500000000000,-0.0600400000000000,-0.0556200000000000,-0.0501300000000000,-0.0438200000000000,-0.0369100000000000,-0.0296200000000000,-0.0222400000000000,-0.0151300000000000,-0.00867000000000000,-0.00334000000000000,0.000160000000000000,0;];
Foil.N63215.US.x = [1,0.950140000000000,0.900300000000000,0.850430000000000,0.800510000000000,0.750550000000000,0.700530000000000,0.650470000000000,0.600350000000000,0.550190000000000,0.500000000000000,0.449770000000000,0.399520000000000,0.349260000000000,0.299000000000000,0.248750000000000,0.198520000000000,0.148340000000000,0.0982300000000000,0.0732300000000000,0.0482900000000000,0.0234800000000000,0.0112000000000000,0.00637000000000000,0.00399000000000000,0;];
Foil.N63215.US.z = [0,0.00616000000000000,0.0136800000000000,0.0221300000000000,0.0310500000000000,0.0401400000000000,0.0490600000000000,0.0575100000000000,0.0652400000000000,0.0720300000000000,0.0776800000000000,0.0819400000000000,0.0845700000000000,0.0853000000000000,0.0839200000000000,0.0804900000000000,0.0748700000000000,0.0668200000000000,0.0556900000000000,0.0484700000000000,0.0396000000000000,0.0279200000000000,0.0198000000000000,0.0152800000000000,0.0125000000000000,0;];
%%  S814 Geometry
Foil.S814.US.x = [0,...
    0.00116,0.0083,0.02064,0.03771,0.05918,0.08475,0.11409,0.14685,...
    0.18266,0.22111,0.26177,0.30418,0.34829,0.39439,0.44237,0.49169,...
    0.54177,0.59199,0.64174,0.69037,0.73723,0.78169,0.82312,0.86095,...
    0.8946,0.9238,0.94879,0.96963,0.98582,0.99632,1];
Foil.S814.LS.x = [0,...
    0.00048,0.00607,0.01644,0.03097,0.04923,0.07077,0.09515,...
    0.12193,0.15072,0.18122,0.21322,0.24712,0.28389,0.32394,0.36753,...
    0.41483,0.46552,0.51909,0.57485,0.63189,0.68912,0.74529,0.79901,...
    0.84887,0.89348,0.93154,0.96197,0.98364,0.99606,1];
Foil.S814.US.z = [0,...
    0.00703,0.01892,0.0313,0.04378,0.05608,0.06791,0.07903,0.08921,...
    0.09821,0.1058,0.11175,0.11564,0.11696,0.11573,0.11251,0.10775,...
    0.10173,0.09473,0.08698,0.07873,0.07016,0.06146,0.05276,...
    0.04417,0.03567,0.02706,0.01848,0.01071,0.0047,0.00112,0];
Foil.S814.LS.z = [0,...
    -0.0047,-0.01746,-0.03159,-0.04646,-0.06162,-0.07662,-0.09096,...
    -0.10412,-0.11545,-0.12425,-0.12971,-0.13079,-0.12736,-0.1199,...
    -0.10887,-0.09511,-0.07962,-0.06328,-0.04703,-0.03173,-0.01818,...
    -0.00701,0.00134,0.00671,0.00917,0.0091,0.00701,0.00377,0.00102,0];


blade.RADIUS = [0.0200,    0.0238,    0.0277,    0.0315,    0.0354,    0.0392,    0.0430,    0.0469,    0.0507,    0.0545,    0.0584,    0.0622,    0.0661,    0.0699,    0.0737,    0.0776,    0.0814,    0.0853,    0.0891,    0.0929,    0.0968,    0.1006,    0.1044,    0.1083,    0.1121,    0.1160,    0.1198,    0.1236,    0.1275,    0.1313,    0.1352,    0.1390,    0.1428,    0.1467,    0.1505,    0.1543,    0.1582,    0.1620,    0.1659,    0.1697,    0.1735,    0.1774,    0.1812,    0.1851,    0.1889,    0.1927,    0.1966,    0.2004,    0.2042,    0.2081,    0.2119,    0.2158,    0.2196,    0.2234,    0.2273,    0.2311,    0.2349,    0.2388,    0.2426,    0.2465,    0.2503,    0.2541,    0.2580,    0.2618,    0.2657,    0.2695,    0.2733,    0.2772,    0.2810,    0.2848,    0.2887,    0.2925,    0.2964,    0.3002,    0.3040,    0.3079,    0.3117,    0.3156,    0.3194,    0.3232,    0.3271,    0.3309,    0.3347,    0.3386,    0.3424,    0.3463,    0.3501,    0.3539,    0.3578,    0.3616,    0.3655,    0.3693,    0.3731,    0.3770,    0.3808,    0.3846,    0.3885,    0.3923,    0.3962,    0.4000];
blade.THICKNESS = [24 24 20.7 18.7 17.6 16.6 15.6 14.6 13.6 12.6];
c = [0.0120,    0.0120,    0.0120,    0.0120,    0.0149,    0.0178,    0.0208,    0.0237,    0.0266,    0.0295,    0.0325,    0.0354,    0.0383,    0.0412,    0.0442,    0.0471,    0.0500,    0.0498,    0.0494,    0.0489,    0.0485,    0.0481,    0.0477,    0.0474,    0.0470,    0.0466,    0.0462,    0.0459,    0.0455,    0.0452,    0.0448,    0.0444,    0.0441,    0.0438,    0.0434,    0.0430,    0.0427,    0.0423,    0.0420,    0.0416,    0.0412,    0.0409,    0.0405,    0.0402,    0.0398,    0.0394,    0.0391,    0.0387,    0.0384,    0.0380,    0.0376,    0.0373,    0.0369,    0.0366,    0.0362,    0.0358,    0.0355,    0.0351,    0.0348,    0.0344,    0.0340,    0.0337,    0.0333,    0.0330,    0.0326,    0.0322,    0.0319,    0.0315,    0.0312,    0.0308,    0.0304,    0.0300,    0.0297,    0.0294,    0.0290,    0.0286,    0.0283,    0.0279,    0.0276,    0.0272,    0.0268,    0.0265,    0.0261,    0.0258,    0.0254,    0.0250,    0.0247,    0.0243,    0.0240,    0.0236,    0.0232,    0.0229,    0.0225,    0.0222,    0.0218,    0.0214,    0.0211,    0.0207,    0.0204,    0.0200];
c(4:17) = linspace(0.03,0.125,14);
blade.THETA = [15.0000,   15.0000,   15.0000,   15.0000,   15.0000,   15.0000,   15.0000,   15.0000,   15.0000,   15.0000,   15.0000,   15.0000,   15.0000,   15.0000,   15.0000,   15.0000,   14.7949,   14.2384,   13.6818,   13.1253,   12.5687,   12.0212,   11.5222,   11.0232,   10.5242,   10.0253,    9.5263,    9.1545,    8.7899,    8.4253,    8.0606,    7.6960,    7.3879,    7.1000,    6.8121,    6.5242,    6.2364,    5.9788,    5.7485,    5.5182,    5.2879,    5.0576,    4.8394,    4.6475,    4.4556,    4.2636,    4.0717,    3.8838,    3.7303,    3.5768,    3.4232,    3.2697,    3.1162,    2.9798,    2.8455,    2.7111,    2.5768,    2.4424,    2.3343,    2.2384,    2.1424,    2.0465,    1.9505,    1.8636,    1.7869,    1.7101,    1.6333,    1.5566,    1.4848,    1.4273,    1.3697,    1.3121,    1.2545,    1.1970,    1.1394,    1.0818,    1.0242,    0.9667,    0.9091,    0.8515,    0.7939,    0.7364,    0.6788,    0.6212,    0.5758,    0.5374,    0.4990,    0.4606,    0.4222,    0.3838,    0.3455,    0.3071,    0.2687,    0.2303,    0.1919,    0.1535,    0.1152,    0.0768,    0.0384,    0.0000];
blade.CHORD = c.*.4;
blade.THICKNESS = interp1([0.0 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0],blade.THICKNESS,blade.RADIUS/0.4);
FoilArea = trapz(fliplr(Foil.N63215.US.x),fliplr(Foil.N63215.US.z/0.15)) - trapz(Foil.N63215.LS.x,Foil.N63215.LS.z/0.15);
blade.AREA = (blade.CHORD.*0.01.*blade.THICKNESS).*FoilArea;

ChordData = [0.15	43.910256
    0.200913	42.948718
    0.248858	41.346154
    0.301370	40.064103
    0.351598	38.461538
    0.401826	36.858974
    0.447489	35.576923
    0.500000	34.134615
    0.598174	31.410256
    0.698630	28.525641
    0.801370	25.801282
    0.904110	22.916667
    1.000000	20.192308];

TwistData = [0.15	32.211538
    0.198630	27.083333
    0.248858	20.032051
    0.301370	15.224359
    0.349315	11.057692
    0.399543	7.532051
    0.449772	5.608974
    0.502283	4.166667
    0.600457	2.243590
    0.698630	0.641026
    0.799087	-0.320513
    0.899543	-1.282051
    1.000	-1.762821];

blade.RADIUS = linspace(0.15,1.0);
blade.THETA = interp1(TwistData(:,1),TwistData(:,2),blade.RADIUS);
blade.CHORD = interp1(ChordData(:,1),ChordData(:,2),blade.RADIUS)/1000;
blade.RADIUS = 0.175*linspace(0.15,1.0);
blade.THICKNESS = [];

        
 FoilArea = trapz((Foil.S814.US.x),(Foil.S814.US.z)) - trapz(Foil.S814.LS.x,Foil.S814.LS.z);
blade.AREA = blade.CHORD.*FoilArea;

rhoSteel = 2700; %kg/m^3
rhoH2O = 998;

MassBlade = blade.AREA * rhoSteel;
ratio = rhoSteel/rhoH2O;
apparentMass = 9.80665*(MassBlade - MassBlade/ratio);

% for i = 1:360
%     
%     
%     SelfMoment = trapz(blade.RADIUS,apparentMass .* blade.RADIUS.* blade.RADIUS .* sin(theta));
% end


% try BodySurface1 with +ve with 15 secs sim time 
for i = 1:size(CpHistory,1)
    inds = BodySurface0(:);
    Cp = CpHistory(i,inds)';
    q = 0.5.*rhoH2O.*sqrt(VCollocPts_x(inds).^2 + VCollocPts_y(inds).^2 + VCollocPts_z(inds).^2);
    
    
    Norms = [PanelNormalHistory0_x(inds,i) PanelNormalHistory0_y(inds,i) PanelNormalHistory0_z(inds,i)];
    
    CollocPts = [CollocPtHistory0_x(:,i) CollocPtHistory0_y(:,i) CollocPtHistory0_z(:,i)];
    F = -[q.*Area(inds).*Cp.*Norms(:,1) q.*Area(inds).*Cp.*Norms(:,2) q.*Area(inds).*Cp.*Norms(:,3)];
    th =0;
   
    
%     EulerAngles.x = EulerHist0_x(i); % roll
%     EulerAngles.y = EulerHist0_y(i); % pitch
%     EulerAngles.z = EulerHist0_z(i); % yaw
%     
%     cosphi = cos(EulerAngles.x); costhe = cos(EulerAngles.y); cospsi = cos(EulerAngles.z);
%     sinphi = sin(EulerAngles.x); sinthe = sin(EulerAngles.y); sinpsi = sin(EulerAngles.z);
%     a1 = costhe*cospsi;
%     a2 = costhe*sinpsi;
%     a3 = -sinthe;
%     b1 = sinphi * sinthe * cospsi - cosphi*sinpsi;
%     b2 = sinphi * sinthe * sinpsi + cosphi*cospsi;
%     b3 = sinphi*costhe;
%     c1 = cosphi * sinthe * cospsi + sinphi*sinpsi;
%     c2 = cosphi * sinthe * sinpsi - sinphi*cospsi;
%     c3 = cosphi*costhe;
%     TRANS = [ a1, b1, c1 ; a2, b2, c2 ; a3, b3, c3];
%     
%     Fxl = a1 * F(:,1) + a2 * F(:,2) + a3 * F(:,3);
%     Fyl = b1 * F(:,1) + b2 * F(:,2) + b3 * F(:,3);
%     Fzl = c1 * F(:,1) + c2 * F(:,2) + c3 * F(:,3);
    %scatter(CollocPts(inds,2),CollocPts(inds,3));
    %drawnow
    M = cross(F,[CollocPts(inds,1) CollocPts(inds,2) CollocPts(inds,3)]);
    SelfMoment(i) = FoilArea * (rhoSteel - rhoH2O) * 9.80665 * cos(deg2rad(th) + BodyRates0_x*Times(i)) * trapz(blade.RADIUS,blade.RADIUS.*blade.RADIUS.*blade.CHORD);
    
    
    Yl = CollocPts_y(inds);
    Zl = CollocPts_z(inds);
    
    YLr = (round(Yl*1000))/1000;
    
    YLru = unique(YLr);
    
    bnds = [0.026; YLru(1:end-1) + 0.5*diff(YLru); 0.175];
    dbnds = diff(bnds);
    for j = 1:length(YLru)
         Fop(j) = sum((YLr==YLru(j)).*F(:,1)/dbnds(j));
    end
    
    
    Mout_of_plane(i) = trapz(YLru,YLru.*Fop');
   % Min_plane(i) = trapz(YLru,YLru.*Fip');
    
    Mx(i) = sum(M(:,1));
    Fx(i) = sum(F(:,1));
    My(i) = sum(Yl.*F(:,1));
    Fy(i) = sum(F(:,2));
    Mz(i) = sum(M(:,3));
    Fz(i) = sum(F(:,3));
    CLift = sum(F(:,3))*cosd(8.5) - sum(F(:,1))*sind(8.5);
    CL(i) = sum(CLift)./12;
%    CDrag = Norms(inds,1).*Area(inds).*(Cp);
 %   CD(i) = sum(CDrag)./12;
    
    
%     scatter3(CollocPts(:,1), CollocPts(:,2), CollocPts(:,3));
% axis equal square
% drawnow
    
    
    %scatter(Yl,Zl)
    %axis equal
    %drawnow
    
    
%    Mucp0 = Mu(BodySurface0);
 %   CpCp0 = Cp;%(BodySurface0);
    
 %   r = 0.6;
    
    
    
   % CPress = zeros(size(Cp(1,:)));
    
    %for j = 1:size(R0,2)
   %     PressChord(j) = interp1(R0(:,j),C0(:,j),r,'cubic');
    %    CPress(j) = interp1(R0(:,j),-CpCp0(:,j),r,'cubic');
   % end
    %clf
    %plot(PressChord,CPress);
    
    %Cl(i) = trapz(PressChord,CPress);
    %drawnow;
    disp(i)
end
hold all


Times * (20.9500/2*pi)


%plot(PressChord,CPress);
return
figure


BodyPanPts = [C1;C2;C3;C4];
PtIDS = 1:4*length(C1);
PtIDS = reshape(PtIDS,length(C1),4);

p = patch('Vertices',BodyPanPts,...
    'Faces',PtIDS,'FaceVertexCData',Cp(:),...
    'FaceColor','flat','EdgeColor','none');
set(gcf,'Renderer','OpenGL')
hold all

WakePanPts = [WakePanC1_x WakePanC1_y WakePanC1_z;
    WakePanC2_x WakePanC2_y WakePanC2_z;
    WakePanC3_x WakePanC3_y WakePanC3_z;
    WakePanC4_x WakePanC4_y WakePanC4_z];

PtIDS = 1:4*length(WakePanC1_x);
PtIDS = reshape(PtIDS,length(WakePanC1_x),4);

p2 = patch('Vertices',WakePanPts,...
    'Faces',PtIDS,'FaceVertexCData',WakePanGamma(:),...
    'FaceColor','flat','EdgeColor','k');


% ProtoWakePanPts = [ProtoWakePointsX(:,1) ProtoWakePointsY(:,1) ProtoWakePointsZ(:,1);
%     ProtoWakePointsX(:,2) ProtoWakePointsY(:,2) ProtoWakePointsZ(:,2);
%     ProtoWakePointsX(:,3) ProtoWakePointsY(:,3) ProtoWakePointsZ(:,3);
%     ProtoWakePointsX(:,4) ProtoWakePointsY(:,4) ProtoWakePointsZ(:,4)];
% 
% 
% 
% 
% PtIDS = 1:4*length(ProtoWakePointsX(:,1));
% PtIDS = reshape(PtIDS,length(ProtoWakePointsX(:,1)),4);
% 
% p3 = patch('Vertices',ProtoWakePanPts,...
%    'Faces',PtIDS,'FaceVertexCData',ProtoWakeGamma(:),...
%    'FaceColor','flat','EdgeColor','k');

axis equal tight
view(3)



scatter3(VortonX_x(:),VortonX_y(:),VortonX_z(:),'+')
% 
figure
plot(Times,Force_x)
%plot(Times,CL,'-b')
hold all
%plot(Times,Cl,'-k')
%plot(Times,2*pi*AlphaHistory,'-r')