function handles = run_main(handles)
here = pwd;
cd ..
[status, result] = system('./main  > dump')


cd(here)