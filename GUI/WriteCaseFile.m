function handles = WriteCaseFile(handles)

name = handles.fullname;
MaxP = int16(str2num(get(handles.maxP,'String')));
Scale = str2num(get(handles.scale,'String'));
numbodies = 1;
if get(handles.blades_as_bodies,'Value')
    numbodies = int16(str2num(get(handles.nblades,'String')));
end


Vinf{1} = [str2num(get(handles.uinf,'String')) str2num(get(handles.vinf,'String')) str2num(get(handles.winf,'String'))];


for i = 1:numbodies
    Attitudes{i} = [0 0 0];%[str2num(get(handles.psi,'String')) str2num(get(handles.theta,'String')) str2num(get(handles.phi,'String'))];
    Origin{i} = [str2num(get(handles.Ox,'String')) str2num(get(handles.Oy,'String')) str2num(get(handles.Oz,'String'))];
    Rates{i} = [str2num(get(handles.p,'String')) str2num(get(handles.q,'String')) str2num(get(handles.r,'String'))];
    Vels{i} = [str2num(get(handles.u,'String')) str2num(get(handles.v,'String')) str2num(get(handles.w,'String'))];
    axis = [str2num(get(handles.rotor_axis_x,'String')) str2num(get(handles.rotor_axis_y,'String')) str2num(get(handles.rotor_axis_z,'String'))];
    axis = axis./(sqrt(dot(axis,axis)));
    rotorspeed = 2*pi*str2num(get(handles.rpm,'String'))/60;
    Rates{i} = Rates{i} + rotorspeed*axis;
    names{i} = [get(handles.name,'String') '_' num2str(i)];
end
    
%%  Write case file 
                        
fid_cas = fopen(['../case_files/' name '.cas'],'wt');
fname = [name '.neu'];
CaseData = struct([]);
CaseData = AddVal2Case(CaseData,'Input Neu File',fname,'string');
CaseData = AddVal2Case(CaseData,'max P',MaxP,'int');
CaseData = AddVal2Case(CaseData,'Scale',Scale,'double');
CaseData = AddVal2Case(CaseData,'Num Bodies',numbodies,'int');
CaseData = AddVal2Case(CaseData,'Body Attitudes', Attitudes,'double');
CaseData = AddVal2Case(CaseData,'Body CGs',Origin, 'double');
CaseData = AddVal2Case(CaseData,'Body Rates',Rates,'double');
CaseData = AddVal2Case(CaseData,'Body Kin-Vels',Vels,'double');
CaseData = AddVal2Case(CaseData,'Case Name',name,'string');
CaseData = AddVal2Case(CaseData,'dt Init',get(handles.dtinit,'String'),'string');
CaseData = AddVal2Case(CaseData,'Num Sub-Steps',get(handles.nss,'String'),'string');
CaseData = AddVal2Case(CaseData,'Vinf',Vinf,'double');
CaseData = AddVal2Case(CaseData,'ds',get(handles.DS,'String'),'string');





% BodyData = struct([]);
% BodyData = AddVal2Case(BodyData,'# Bodies',handles.BodyData.numbodies,'int');
% BodyData = AddVal2Case(BodyData,'Attitudes', handles.BodyData.attitudes,'double');
% BodyData = AddVal2Case(BodyData,'CGs',handles.BodyData.origin, 'double');
% BodyData = AddVal2Case(BodyData,'Body Rates',handles.BodyData.rates,'double');
% BodyData = AddVal2Case(BodyData,'Body Vels',handles.BodyData.vels,'double');
% BodyData = AddVal2Case(BodyData,'Name',handles.BodyData.names,'string');


% CaseInfo = struct([]);
% CaseInfo = AddVal2Case(CaseInfo,'Input File',fname,'string');
% CaseInfo = AddVal2Case(CaseInfo,'dtInit',handles.CaseData.dtinit,'double');
% CaseInfo = AddVal2Case(CaseInfo,'dtOut',handles.CaseData.dtout,'double');
% CaseInfo = AddVal2Case(CaseInfo,'tMax',handles.CaseData.tmax,'double');
% CaseInfo = AddVal2Case(CaseInfo,'NumSubSteps',handles.CaseData.nss,'int');
% CaseInfo = AddVal2Case(CaseInfo,'CFLmax',handles.CaseData.cflmax,'double');
% CaseInfo = AddVal2Case(CaseInfo,'MaxP',handles.CaseData.maxp,'int');
% CaseInfo = AddVal2Case(CaseInfo,'Scheme',handles.CaseData.fvmscheme,'String');
% CaseInfo = AddVal2Case(CaseInfo,'Limiter',handles.CaseData.limiter,'String');
% CaseInfo = AddVal2Case(CaseInfo,'Integrator',handles.CaseData.integrator,'String');







fprintf(fid_cas,'#          Input file generated by MATLAB Version:  %s\n', version);
fprintf(fid_cas,'#          %s\n', datestr(now));
fprintf(fid_cas,'#\t%s\t%s\t%s\t%s\t%s\n', 'dtInit','dtOut','nss','CFLmax','tMax');
fprintf(fid_cas,'#\t%s\t%s\t%s\t%s\t%s\n', 'dtInit','dtOut','nss','CFLmax','tMax');
for i = 1:CaseData.NumVals
    l = numel(CaseData.Data{i}.Name);
    d = 24 - l;
    spc = repmat(' ',1,d);
    fprintf(fid_cas,'%s:%s%s;\n',CaseData.Data{i}.Name, spc, CaseData.Data{i}.Val);
end

fclose(fid_cas);

fid = fopen(['../case_files/' name '.cas']);

tline = fgetl(fid);

handles.casemsg{1} = 'The following is a transcript of the case file:';
linen = 2;
handles.casemsg{linen} = ['>   '  tline];
while ischar(tline)
    linen = linen + 1;
    tline = fgetl(fid);
    handles.casemsg{linen} = ['>   '  tline];
end

fclose(fid);
