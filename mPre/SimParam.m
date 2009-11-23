function [handles CaseInfo] = SimParam(handles)

handles.CaseData.name = get(handles.name,'String');
handles.CaseData.dtinit = str2double(get(handles.dtinit,'String'));
handles.CaseData.dtout = str2double(get(handles.dtout,'String'));
handles.CaseData.tmax = str2double(get(handles.tmax,'String'));
handles.CaseData.nss = str2num(get(handles.nss,'String'));
handles.CaseData.maxp = str2num(get(handles.maxP,'String'));
handles.CaseData.cflmax = str2double(get(handles.cflmax,'String'));
handles.CaseData.scale = str2double(get(handles.scale,'String'));
handles.CaseData.vinf = [str2double(get(handles.uinf,'String')) str2double(get(handles.vinf,'String')) str2double(get(handles.winf,'String'))];
handles.CaseData.nu = str2double(get(handles.nu,'String'));
handles.CaseData.rho = str2double(get(handles.rho,'String'));

t = get(handles.scheme,'String');
handles.CaseData.fvmscheme  = t(get(handles.scheme,'Value'));

if strcmp(handles.CaseData.fvmscheme,'O(2) MUSCL')
    handles.CaseData.limiter = get(handles.minmodbeta,'String');
else
    t = get(handles.limiter,'String');
    handles.CaseData.limiter  = t(get(handles.limiter,'Value'));
end

t = get(handles.integrator,'String');
handles.CaseData.integrator  = t(get(handles.integrator,'Value'));

CaseInfo = struct([]);
CaseInfo = AddVal2Case(CaseInfo,'Input File',handles.CaseData.name,'string');
CaseInfo = AddVal2Case(CaseInfo,'dtInit',handles.CaseData.dtinit,'double');
CaseInfo = AddVal2Case(CaseInfo,'dtOut',handles.CaseData.dtout,'double');
CaseInfo = AddVal2Case(CaseInfo,'tMax',handles.CaseData.tmax,'double');
CaseInfo = AddVal2Case(CaseInfo,'NumSubSteps',handles.CaseData.nss,'int');
CaseInfo = AddVal2Case(CaseInfo,'CFLmax',handles.CaseData.cflmax,'double');
CaseInfo = AddVal2Case(CaseInfo,'MaxP',handles.CaseData.maxp,'int');
CaseInfo = AddVal2Case(CaseInfo,'Scheme',handles.CaseData.fvmscheme,'String');
CaseInfo = AddVal2Case(CaseInfo,'Limiter',handles.CaseData.limiter,'String');
CaseInfo = AddVal2Case(CaseInfo,'Integrator',handles.CaseData.integrator,'String');