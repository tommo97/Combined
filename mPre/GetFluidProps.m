function handles = GetFluidProps(handles)

%%  Use temperature and fluid type to get viscocity
tempr = str2num(get(handles.fluid_temp,'String'));

switch handles.fluid.type
    case 'Air' % User selects Peaks.
        %%  NRELBlade -- NREL data
        handles.fluid.rho = 1.225;
        c = 120;
        t0 = 291.15;
        mu0 = 18.27e-6;
        mu = Sutherland(mu0,c,tempr,t0);
    case 'Fresh Water' 
        T = [100 80 60 40 30 25 22 20 15 10  4 0];
        R = [958.4  971.8  983.2  992.2  995.6502  997.0479  997.7735  998.2071  999.1026 999.7026  999.9720  999.8395];
        handles.fluid.rho = interp1(T,R,tempr-273.15,'cubic');
        a = 0.585;
        b = -12;
        c = 0.03361;
        d = 1.235;
        t = tempr - 273.15;
        handles.fluid.nu = ITTC(a,b,c,d,t);
    case 'Sea Water'
        handles.fluid.rho = 1027;
        a = 0.659;
        b = 1;
        c = -0.05076;
        d = 1.7688;
        t = tempr - 273.15;
        handles.fluid.nu = ITTC(a,b,c,d,t);
    otherwise
        handles.fluid.rho = 1;
end

set(handles.rho,'String',num2str(handles.fluid.rho));
set(handles.nu,'String',num2str(handles.fluid.nu));

function mu = Sutherland(mu0,c,t,t0)
mu = mu0 * ((t0 + c)/(t + c)) * (t/t0)^1.5;

function nu = ITTC(a,b,c,d,t)
nu = 1e-6 * (((a*1e-3) * (t - b) - c) * (t - b) + d);