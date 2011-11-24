function blade = SetBladeData(blade)
% Set current data to the selected data set.
blade.REVERSE = false;
blade.isNREL = false;   %   for transisition piece
blade.isSOTON = false;  %   for transisition piece
blade.TransitionPiece = [];
blade.SWEEP = [];
switch blade.type;
    case 'NREL UAE' % User selects Peaks.
        %%  NRELBlade -- NREL data
        blade.RADIUS=[    0.5080
            0.5226
            0.5372
            0.5518
            0.5663
            0.5809
            0.5955
            0.6101
            0.6247
            0.6393
            0.6539
            0.6684
            0.6830
            0.6976
            0.7122
            0.7268
            0.7414
            0.7560
            0.7705
            0.7851
            0.7997
            0.8143
            0.8289
            0.8435
            0.8581
            0.8726
            0.8872
            0.9018
            0.9164
            0.9310
            0.9456
            0.9602
            0.9747
            0.9893
            1.0039
            1.0185
            1.0331
            1.0477
            1.0623
            1.0768
            1.0914
            1.1060
            1.1206
            1.1352
            1.1498
            1.1644
            1.1789
            1.1935
            1.2081
            1.2227
            1.2373
            1.2519
            1.2665
            1.2811
            1.2956
            1.3102
            1.3248
            1.3394
            1.3540
            1.3686
            1.3832
            1.3977
            1.4123
            1.4269
            1.4415
            1.4561
            1.4707
            1.4853
            1.4998
            1.5144
            1.5290
            1.5436
            1.5582
            1.5728
            1.5874
            1.6019
            1.6165
            1.6311
            1.6457
            1.6603
            1.6749
            1.6895
            1.7040
            1.7186
            1.7332
            1.7478
            1.7624
            1.7770
            1.7916
            1.8061
            1.8207
            1.8353
            1.8499
            1.8645
            1.8791
            1.8937
            1.9082
            1.9228
            1.9374
            1.9520;2.257;...
            2.343;2.562;2.867;3.172;3.185;3.476;3.781;4.023;4.086;4.391;4.696;4.78;5;5.305;5.532];
        blade.CHORD=[    0.2180
            0.2180
            0.2180
            0.2180
            0.2180
            0.2180
            0.2180
            0.2180
            0.2180
            0.2180
            0.2180
            0.2183
            0.2203
            0.2239
            0.2292
            0.2360
            0.2443
            0.2539
            0.2648
            0.2768
            0.2900
            0.3041
            0.3191
            0.3349
            0.3515
            0.3686
            0.3863
            0.4045
            0.4230
            0.4417
            0.4607
            0.4797
            0.4987
            0.5176
            0.5363
            0.5547
            0.5727
            0.5903
            0.6074
            0.6237
            0.6394
            0.6542
            0.6681
            0.6810
            0.6928
            0.7034
            0.7127
            0.7207
            0.7271
            0.7321
            0.7353
            0.7369
            0.7368
            0.7358
            0.7342
            0.7322
            0.7302
            0.7284
            0.7269
            0.7254
            0.7239
            0.7224
            0.7209
            0.7194
            0.7179
            0.7165
            0.7150
            0.7135
            0.7120
            0.7106
            0.7091
            0.7076
            0.7061
            0.7046
            0.7032
            0.7017
            0.7002
            0.6987
            0.6972
            0.6958
            0.6943
            0.6928
            0.6913
            0.6898
            0.6883
            0.6868
            0.6854
            0.6839
            0.6824
            0.6809
            0.6794
            0.6779
            0.6764
            0.6749
            0.6735
            0.6720
            0.6705
            0.6690
            0.6675
            0.6660;0.636;...
            0.627;0.605;0.574;0.543;0.542;0.512;0.482;0.457;0.451;0.42;0.389;0.381;0.358;0.328;0.305];
        
        RADIUS=[0.508;0.66;1.257;1.343;1.51;1.648;1.952;2.257;...
            2.343;2.562;2.867;3.172;3.185;3.476;3.781;4.023;4.086;4.391;4.696;4.78;5;5.305;5.532];
        blade.THETA=[0;0;20.04;18.074;14.292;11.909;7.979;5.308;4.715;...
            3.425;2.083;1.15;1.115;0.494;-0.015;-0.381;-0.475;-0.92;-1.352;-1.469;-1.775;-2.191;-2.5];
        blade.THETA = interp1(RADIUS,blade.THETA,blade.RADIUS,'cubic');
        blade.THETA = blade.THETA  - blade.THETA(end);
        blade.isNREL = true;
        blade.THICKNESS = [];
    case 'SOTON' % User selects Membrane.
        r = linspace(.2, 1, 17);
        blade.RADIUS = [    0.0200
    0.0238
    0.0277
    0.0315
    0.0354
    0.0392
    0.0430
    0.0469
    0.0507
    0.0545
    0.0584
    0.0622
    0.0661
    0.0699
    0.0737
    0.0776
    0.0814
    0.0853
    0.0891
    0.0929
    0.0968
    0.1006
    0.1044
    0.1083
    0.1121
    0.1160
    0.1198
    0.1236
    0.1275
    0.1313
    0.1352
    0.1390
    0.1428
    0.1467
    0.1505
    0.1543
    0.1582
    0.1620
    0.1659
    0.1697
    0.1735
    0.1774
    0.1812
    0.1851
    0.1889
    0.1927
    0.1966
    0.2004
    0.2042
    0.2081
    0.2119
    0.2158
    0.2196
    0.2234
    0.2273
    0.2311
    0.2349
    0.2388
    0.2426
    0.2465
    0.2503
    0.2541
    0.2580
    0.2618
    0.2657
    0.2695
    0.2733
    0.2772
    0.2810
    0.2848
    0.2887
    0.2925
    0.2964
    0.3002
    0.3040
    0.3079
    0.3117
    0.3156
    0.3194
    0.3232
    0.3271
    0.3309
    0.3347
    0.3386
    0.3424
    0.3463
    0.3501
    0.3539
    0.3578
    0.3616
    0.3655
    0.3693
    0.3731
    0.3770
    0.3808
    0.3846
    0.3885
    0.3923
    0.3962
    0.4000];
        blade.THICKNESS = [ 24 24 22.5 20.7 19.5 18.7 18.1 17.6 17.1 16.6 16.1 15.6 15.1 14.6 14.1 13.6 13.1 12.6];
        c = [     0.0375
    0.0375
    0.0375
    0.0377
    0.0403
    0.0453
    0.0522
    0.0607
    0.0701
    0.0801
    0.0901
    0.0997
    0.1084
    0.1157
    0.1212
    0.1244
    0.1250
    0.1244
    0.1235
    0.1223
    0.1212
    0.1202
    0.1193
    0.1184
    0.1175
    0.1165
    0.1156
    0.1147
    0.1138
    0.1129
    0.1120
    0.1111
    0.1102
    0.1094
    0.1085
    0.1076
    0.1067
    0.1058
    0.1049
    0.1040
    0.1031
    0.1022
    0.1013
    0.1004
    0.0995
    0.0986
    0.0977
    0.0968
    0.0959
    0.0950
    0.0941
    0.0932
    0.0923
    0.0914
    0.0905
    0.0896
    0.0887
    0.0878
    0.0869
    0.0860
    0.0851
    0.0842
    0.0833
    0.0824
    0.0815
    0.0806
    0.0797
    0.0788
    0.0779
    0.0770
    0.0761
    0.0751
    0.0742
    0.0734
    0.0725
    0.0716
    0.0707
    0.0698
    0.0689
    0.0680
    0.0671
    0.0662
    0.0653
    0.0644
    0.0635
    0.0626
    0.0617
    0.0608
    0.0599
    0.0590
    0.0581
    0.0572
    0.0563
    0.0554
    0.0545
    0.0536
    0.0527
    0.0518
    0.0509
    0.0500];
        blade.THETA = [   15.0000
   15.0000
   15.0000
   15.0000
   15.0000
   15.0000
   15.0000
   15.0000
   15.0000
   15.0000
   15.0000
   15.0000
   15.0000
   15.0000
   15.0000
   15.0000
   14.7949
   14.2384
   13.6818
   13.1253
   12.5687
   12.0212
   11.5222
   11.0232
   10.5242
   10.0253
    9.5263
    9.1545
    8.7899
    8.4253
    8.0606
    7.6960
    7.3879
    7.1000
    6.8121
    6.5242
    6.2364
    5.9788
    5.7485
    5.5182
    5.2879
    5.0576
    4.8394
    4.6475
    4.4556
    4.2636
    4.0717
    3.8838
    3.7303
    3.5768
    3.4232
    3.2697
    3.1162
    2.9798
    2.8455
    2.7111
    2.5768
    2.4424
    2.3343
    2.2384
    2.1424
    2.0465
    1.9505
    1.8636
    1.7869
    1.7101
    1.6333
    1.5566
    1.4848
    1.4273
    1.3697
    1.3121
    1.2545
    1.1970
    1.1394
    1.0818
    1.0242
    0.9667
    0.9091
    0.8515
    0.7939
    0.7364
    0.6788
    0.6212
    0.5758
    0.5374
    0.4990
    0.4606
    0.4222
    0.3838
    0.3455
    0.3071
    0.2687
    0.2303
    0.1919
    0.1535
    0.1152
    0.0768
    0.0384
         0];
        blade.isSOTON = true;
        blade.CHORD = c.*.4;
        blade.THICKNESS = [];
    case 'ESRU PoC 1' % User selects Sinc.
        blade.RADIUS=[0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1,1.050000,1.1,1.150000,1.2,1.250000;];
        blade.CHORD=[0.100,0.0965,0.0930,0.0895,0.0860,0.0825,0.0790,0.0755,0.0720,0.0685,0.0650,0.0615,0.0580,0.0545,0.0510,0.0475,0.0440,0.0405,0.0370,0.0335,0.0300;];
        blade.THETA=[35.48354,30.66771,26.60373,23.19068,20.33544,17.95276,15.96524,14.30328,12.90517,11.71700,10.69272,9.794129,8.990847,8.260345,7.587935,6.966768,6.397837,5.889976,5.459857,5.131995,4.938747;];
        blade.THICKNESS = [];
    case 'ESRU PoC 2' % User selects Sinc.
        blade.RADIUS=[0.2500;0.3000;0.3500;0.4000;0.4500;0.5000;0.5500;0.6000;0.6500;0.7000;0.7500;0.8000;0.8500;0.9000;0.9500;1;1.050;1.100;1.150;1.200;1.250;];
        blade.CHORD=0.01*[23.30;20.60;18.05;15.85;13.90;12.30;10.90;10;9.200;8.550;8;7.600;7.250;6.850;6.450;6.050;5.550;5;4.400;3.800;3.200;];
        blade.THETA= [26.40;23.25;20.50;18;15.80;13.75;12;10.55;9.550;8.700;8;7.400;6.850;6.350;5.800;5.250;4.600;4;3.400;2.800;2.200;];
        blade.THICKNESS = [];
        
    case 'Elliptic'
        blade.RADIUS = 0.1*linspace(-0.5,0.5);
        blade.THETA = 10*zeros(size(blade.RADIUS));
        blade.CHORD = sqrt(0.125*0.125 - linspace(-0.125,0.125));
        blade.THICKNESS = [];
    case 'Straight'
        blade.RADIUS = linspace(-6,6);
        blade.THETA = 10*zeros(size(blade.RADIUS));
        blade.CHORD = 1*ones(size(blade.RADIUS));
        blade.THICKNESS = [];
        
        
        th = linspace(0,pi/2);
        
        x = cos(th);
        y = sin(th);
        
        x = 1-x;
        
        x = [x linspace(1,4)];
        
        y = [y ones(size(y))];
        
        % x = [x fliplr(4-x)];
        % y = [y fliplr(y)];
        
        blade.RADIUS = linspace(-10,10);
        blade.THETA =   0*ones(size(blade.RADIUS));
        blade.CHORD =  1*ones(size(blade.RADIUS));
        blade.THICKNESS = [];
        
    case 'Wing'
        blade.RADIUS = linspace(-10,10)/5;
        blade.THETA =   [linspace(5,0,50) linspace(0,5,50)];
        
        theta = [4.4  4.3 3.8 3.5 3 2.9 2.5 2 1.9 1.5 0.6 0.1 -0.6 -2 -5.5];
        
        r = linspace(0.1, 0.975,numel(theta))*max(blade.RADIUS);
        theta = -[fliplr(theta) theta];
        r = [fliplr(-r) r];
        blade.THETA = interp1(r,theta,blade.RADIUS,'cubic','extrap');
        blade.THETA(20:40) = blade.THETA(20:40)-10;
        blade.THETA(60:80) = blade.THETA(60:80)-10;
        
        blade.CHORD =  [linspace(3.983,17.7,50)/17.7 linspace(17.7,3.983,50)/17.7];
        blade.SWEEP =  sind(33.5)*[linspace(10,0,50) linspace(0,10,50)]/5;
end


minrad = min(blade.RADIUS);
maxrad = max(blade.RADIUS);
blade.DistPanel.maxx = min(maxrad,min(maxrad, blade.Cutout.Tip));
blade.DistPanel.minx = max(minrad,max(blade.Cutout.Root, minrad));


if blade.isNREL && ~isempty(blade.y) && (blade.Cutout.Root < 1.257) && (blade.Cutout.Root >= 0.45)
    if (blade.Cutout.Root <= 0.660)
        blade.n1 = interp1(blade.y,1:length(blade.y),0.660,'nearest');
        blade.n2 = interp1(blade.y,1:length(blade.y),1.257,'nearest');
        
        blade.y(1:blade.n1) = linspace(blade.Cutout.Root,0.660,blade.n1);
        blade.y(blade.n1:blade.n2) = linspace(max(blade.Cutout.Root,0.660),1.257,1+blade.n2-blade.n1);
        
        
        blade.TransitionPiece = zeros(size(blade.y));
        blade.TransitionPiece(1:blade.n1) = 1;
        if (blade.n2-blade.n1) > 1
            blade.TransitionPiece(blade.n1:blade.n2) = [logspace(0,-1,blade.n2-blade.n1) 0];
        end
    end
end

if blade.isSOTON && ~isempty(blade.y) && (blade.Cutout.Root < 0.08) && (blade.Cutout.Root >= 0.02)
    if (blade.Cutout.Root < 0.08)
        
        blade.n2 = interp1(blade.y,1:length(blade.y),0.08,'nearest');
        
        %blade.y(1:blade.n1) = linspace(blade.Cutout.Root,0.05,blade.n1);
        blade.y(1:blade.n2+1) = linspace(max(blade.Cutout.Root,0.02),0.08,blade.n2+1);
        
        blade.y(blade.n2+1:end) = linspace(0.08,0.4,numel(blade.y(blade.n2+1:end)));
        blade.TransitionPiece = zeros(size(blade.y));
        blade.TransitionPiece(1:blade.n2) = 1;
        %     if (blade.n2-blade.n1) > 1
        %     blade.TransitionPiece(blade.n1:blade.n2) = [logspace(0,-1,blade.n2-blade.n1) 0];
        %     end
    end
end


if ~isempty(blade.y)
    blade.Radius = blade.y;
    
    blade.Chord = interp1(blade.RADIUS,blade.CHORD,blade.Radius,'linear');
    blade.Theta = interp1(blade.RADIUS,blade.THETA,blade.Radius,'linear');
    blade.Thickness = [];
    if ~isempty(blade.THICKNESS)
        blade.Thickness = interp1(blade.RADIUS,blade.THICKNESS,blade.Radius,'cubic');
    end
    if ~isempty(blade.SWEEP)
        blade.Sweep = interp1(blade.RADIUS,blade.SWEEP,blade.Radius,'cubic');
    else
            blade.Sweep = zeros(size(blade.Radius));
    end

end

PlotBlade(blade);

function PlotBlade(blade)
hold(blade.axes,'off')
if ~isempty(blade.y)
    scatter(blade.axes,blade.Radius,blade.Chord);
    hold(blade.axes,'on');
end
if ~isempty(blade.THICKNESS)
    plot(blade.axes,blade.RADIUS,blade.THICKNESS/100);
    hold(blade.axes,'on');
end
plotyy(blade.axes,blade.RADIUS,blade.CHORD,blade.RADIUS,blade.THETA);

