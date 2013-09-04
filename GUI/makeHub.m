function makeHub(Xhub, Yhub, Zhub)

NBlades = length(Xhub);

Angle = 2*pi/NBlades;

for B = 1:NBlades
    l = (1+(length(Xhub{B})-1)/2);
    
    [THETA,RHO,Z] = cart2pol(Zhub{B},Yhub{B},Xhub{B});
    Radius = mean(RHO);
    
    XhubLower = (Xhub{B}(1:l));
    YhubLower = (Yhub{B}(1:l));
    ZhubLower = (Zhub{B}(1:l));
    
    
    XhubUpper = flipud(Xhub{B}(l:end));
    YhubUpper = flipud(Yhub{B}(l:end));
    ZhubUpper = flipud(Zhub{B}(l:end));
    
    MeanX = 0.5*(XhubUpper + XhubLower);
    MeanY = 0.5*(YhubUpper + YhubLower);
    MeanZ = 0.5*(ZhubUpper + ZhubLower);
    
    
    % want a patch between upper surface and lower surface of next blade, ie at
    % Angle around hub. To make a smooth transition want to blend from US to
    % Meanline to LS during fitting. Assume X is correct
    
    n = 12;
    
    
    
    ThetaUS = atan2(ZhubUpper, YhubUpper);
    ThetaLS = atan2(ZhubLower, YhubLower);
    ThetaMean = 0.5*(ThetaUS+ThetaLS);
    
    TwistFit = FitData(MeanX,ThetaMean);
    
    dxfront = MeanX(end - 1) - MeanX(end);
    dxback = MeanX(1) - MeanX(2);
    
    Thtaback = TwistFit(MeanX(1) + dxback);
    Thtafront = TwistFit( MeanX(end) - dxfront);
    clear THETA RHO Z XH;
    
    
    for i = 1:l
        THETA(i,:) = linspace(ThetaUS(i), ThetaLS(i) + Angle,n)';
        XH(i,:) = linspace(XhubUpper(i), XhubLower(i),n)';
    end
    thtas = linspace(0,pi/2,n);

    R = ones(size(THETA))*Radius;
    
    XS = Radius * cos(thtas(1:end-1))';
    Xadd = MeanX(1) + repmat(XS,[1 n]);
    Xsub = MeanX(end) - repmat(flipud(XS),[1 n]);
    RFront = Radius * sin(thtas(1:end-1));
    RBack = fliplr(RFront);
    
    R = [repmat(RFront',[1 n]); R; repmat(RBack',[1 n])];
    
    ThetaAdd = zeros(size(Xadd));
    ThetaSub = zeros(size(Xsub));
    TwistFit = FitData(MeanX,ThetaMean);
    for i = 1:n
        TwistFit = FitData(XH(:,i),THETA(:,i));
        ThetaAdd(:,i) = TwistFit(Xadd(:,i));
        ThetaSub(:,i) = TwistFit(Xsub(:,i));% + THETA(end,i);
    end
    
    TwistFit = FitData(MeanX,ThetaMean);
    ThetaAdd(:,1) = TwistFit(Xadd(:,1));
    ThetaSub(:,1) = TwistFit(Xsub(:,1));% + THETA(end,i);
    ThetaAdd(:,end) = TwistFit(Xadd(:,end));
    ThetaSub(:,end) = TwistFit(Xsub(:,end));% + THETA(end,i);
    
    
    XH = [Xadd; XH; Xsub];
    %THETA = [ThetaAdd; THETA; ThetaSub];
    dThetadX = -(ThetaMean(end) - ThetaMean(1)) ./ (MeanX(end) - MeanX(1))

    ThetaFront = repmat(THETA(1,:),[n-1 1]) + (MeanX(1) - Xadd) * dThetadX;
    ThetaBack = repmat(THETA(1,:),[n-1 1]) + (MeanX(1) - Xsub) * dThetadX;
    %THETA = [repmat(THETA(1,:),[n-1, 1]); THETA; repmat(THETA(end,:),[n-1, 1])];
    THETA = [ThetaFront;THETA;ThetaBack];
    
    [Y Z X] = pol2cart(THETA,R,XH);
    
    
    
    Xout{B} = X;
    Yout{B} = Y;
    Zout{B} = Z;
    surf(X,Y,Z,'facecolor','w');
    set(gcf,'Renderer','OpenGL');
    clear X Y Z THETA RHO R XH
end

MakeNonWakeSheddingNeuFromSurface(Xout,Yout,Zout);


return;




[THETA,RHO,Z] = cart2pol(Zhub,Yhub,Xhub);

Radius = mean(RHO);
clear THETA RHO Z;

[x y z] = sphere(32);


X = z;
Y = x;
Z = y;

x = X;
y = Y;
z = Z;

midpoint = find(x(:,1)==0);


%MakeNonWakeSheddingNeuFromSurface((z-1),.25*x,.25*y)


theta_insert = repmat(atan2(MeanZ,MeanY),[1 size(x,2)]);
theta = [atan2(MeanZ(end),MeanY(end)) + zeros(size(x(1:midpoint,:)));flipud(theta_insert);atan2(MeanZ(1),MeanY(1)) + zeros(size(x(midpoint:end,:)))];

x_insert = repmat(MeanX,[1 size(x,2)]);
x = [min(MeanX) + Radius * x(1:midpoint,:);flipud(x_insert); max(MeanX) + Radius * x(midpoint:end,:)];

y_insert = repmat(y(midpoint,:),[size(x_insert,1) 1]);
y = Radius * [y(1:midpoint,:);y_insert; y(midpoint:end,:)];
z_insert = repmat(z(midpoint,:),[size(x_insert,1) 1]);
z = Radius * [z(1:midpoint,:);z_insert; z(midpoint:end,:)];


clear X Y Z;
X{1} = x;
Y{1} = y.*cos(theta) + z.*sin(theta);
Z{1} = y.*sin(theta) - z.*cos(theta);

MakeNonWakeSheddingNeuFromSurface(X,Y,Z);


    function Fit = FitData(xdata,ydata)
% Set up chord fit and options.
[xData, yData] = prepareCurveData(xdata, ydata);
ft = fittype( 'pchipinterp' );
opts = fitoptions( ft );
opts.Normalize = 'on';
Fit = fit( xData, yData, ft, opts );