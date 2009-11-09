function [Aerofoil z] = NRELFoil(xin)



%%  S814 Geometry
Aerofoil.S814.US.x = [0,...
        0.00116,0.0083,0.02064,0.03771,0.05918,0.08475,0.11409,0.14685,...
        0.18266,0.22111,0.26177,0.30418,0.34829,0.39439,0.44237,0.49169,...
        0.54177,0.59199,0.64174,0.69037,0.73723,0.78169,0.82312,0.86095,...
        0.8946,0.9238,0.94879,0.96963,0.98582,0.99632,1];
Aerofoil.S814.LS.x = [0,...
        0.00048,0.00607,0.01644,0.03097,0.04923,0.07077,0.09515,...
        0.12193,0.15072,0.18122,0.21322,0.24712,0.28389,0.32394,0.36753,...
        0.41483,0.46552,0.51909,0.57485,0.63189,0.68912,0.74529,0.79901,...
        0.84887,0.89348,0.93154,0.96197,0.98364,0.99606,1];
Aerofoil.S814.US.z = [0,...
        0.00703,0.01892,0.0313,0.04378,0.05608,0.06791,0.07903,0.08921,...
        0.09821,0.1058,0.11175,0.11564,0.11696,0.11573,0.11251,0.10775,...
        0.10173,0.09473,0.08698,0.07873,0.07016,0.06146,0.05276,...
        0.04417,0.03567,0.02706,0.01848,0.01071,0.0047,0.00112,0];
Aerofoil.S814.LS.z = [0,...
        -0.0047,-0.01746,-0.03159,-0.04646,-0.06162,-0.07662,-0.09096,...
        -0.10412,-0.11545,-0.12425,-0.12971,-0.13079,-0.12736,-0.1199,...
        -0.10887,-0.09511,-0.07962,-0.06328,-0.04703,-0.03173,-0.01818,...
        -0.00701,0.00134,0.00671,0.00917,0.0091,0.00701,0.00377,0.00102,0];


%%  S809 Geometry


Aerofoil.S809.US.x = [0,0.0004,0.0058,0.0163,0.0316,0.0515,0.0757,0.1039,0.1358,0.171,0.2092,...
    0.2499,0.2926,0.3369,0.3822,0.4281,0.4738,0.52,0.568,0.6175,0.6672,...
    0.7161,0.7631,0.8076,0.8485,0.8854,0.9176,0.9452,0.968,0.9853,...
    0.9962,1,];
Aerofoil.S809.US.z = [0,0.0028,0.0117,0.0213,0.0314,0.0414,0.0513,0.0608,0.0697,0.0779,0.085,...
    0.0911,0.0959,0.0993,0.1011,0.101,0.0984,0.0924,0.0836,0.0738,0.064,...
    0.0546,0.0458,0.0376,0.0302,0.0234,0.0169,0.011,0.006,0.0024,0.0005,0,];
Aerofoil.S809.LS.x = [0,0.0014,0.0093,0.0232,0.0422,0.0658,0.0932,0.124,0.1575,0.1936,0.2318,...
    0.2713,0.3119,0.3533,0.3954,0.4383,0.4823,0.5284,0.5766,0.6265,0.6771,...
    0.7275,0.7767,0.8235,0.8668,0.9054,0.9385,0.9651,0.9845,0.9961,1,];
Aerofoil.S809.LS.z = [0,-0.005,-0.0127,-0.0216,-0.0314,-0.042,-0.053,-0.0641,-0.0747,-0.0845,...
    -0.0933,-0.1006,-0.1059,-0.1087,-0.1084,-0.1048,-0.0976,-0.087,-0.0744,...
    -0.0611,-0.0479,-0.0356,-0.0247,-0.0156,-0.0086,-0.0037,-0.0008,0.0005,...
    0.0006,0.0002,0,];


%%  NACA 00xx Geometry
x = fzero(@N00xx,[.5 1.5]);
Aerofoil.N00xx.Coeffts = 5*[.2969, -.1260, -.3516, .2843, -.1015];
Aerofoil.N0012.US.x = linspace(0,1);
Aerofoil.N0012.US.z = .12*N00xx(linspace(0,x)); 
Aerofoil.N0012.LS.x = linspace(0,1);
Aerofoil.N0012.LS.z = - Aerofoil.N0012.US.z;    

%%  output

z.N0012.US = interp1(Aerofoil.N0012.US.x,Aerofoil.N0012.US.z,xin,'cubic');
z.N0012.LS = interp1(Aerofoil.N0012.LS.x,Aerofoil.N0012.LS.z,xin,'cubic');
z.N0012.X = xin;
z.S814.US = interp1(Aerofoil.S814.US.x,Aerofoil.S814.US.z,xin,'cubic');
z.S814.LS = interp1(Aerofoil.S814.LS.x,Aerofoil.S814.LS.z,xin,'cubic');
z.S814.X = xin;
z.S809.US = interp1(Aerofoil.S809.US.x,Aerofoil.S809.US.z,xin,'cubic');
z.S809.LS = interp1(Aerofoil.S809.LS.x,Aerofoil.S809.LS.z,xin,'cubic');
z.S809.X = xin;
function y = N00xx(x)
    y = 5*[.2969, -.1260, -.3516, .2843, -.1015]*[sqrt(x);x;x.^2;x.^3;x.^4];