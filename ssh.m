isdone = 0;
while isdone < 2
[a b] = system('ssh lap05140@masternode.mecheng.strath.ac.uk ls output/5062/dump_00*.*');
if (str2num(b(end-5:end-3)) == 300) 
    isdone = isdone + 1;
    system('ssh lap05140@masternode.mecheng.strath.ac.uk scp output/5062/dump_00*.* tom@130.159.102.221:~/Desktop/PanelProject/output/HPC1')
end
[a b] = system('ssh lap05140@masternode.mecheng.strath.ac.uk ls output/10309/dump_00*.*');
if (str2num(b(end-5:end-3)) == 300) 
    isdone = isdone + 1;
    system('ssh lap05140@masternode.mecheng.strath.ac.uk scp output/10309/dump_00*.* tom@130.159.102.221:~/Desktop/PanelProject/output/HPC2')
end
pause(1)
disp(isdone)
end
