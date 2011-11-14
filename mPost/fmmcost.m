clear all
close all
clc

mmax = 12;
%   Calculating leaf moments
count = 0;
for pmax = 1:10
    for k1 = 1:pmax
        
        for k2 = 1:pmax-k1
            
            %for k3 = 1:(pmax-(k1+k2))
                
                count = count + 1;
            %end
        end
    end
    cst(pmax) = count;
end



%   Calculating parent moments

for pmax = 1:10
    for mlevs = 1:mmax
        count = 0;
        for m = 1:mlevs
            
            for ix = 1:2
                
                for iy = 1:2
                    
                    for iz = 1:2
                        
                        for n1 = 1:pmax
                            
                            for n2 = 1:(pmax - n1)
                                
                                for n3 = 1:(pmax - (n1+n2))
                                    
                                    for k1 = 1:(n1 + 1)
                                        
                                        for k2 = 1:(n2 + 1)
                                            
                                            for k3 = 1:(n3 + 1)
                                                
                                                count = count + 1;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        disp(count);
        cstb(pmax,mlevs) = count;
    end
end



%   Calculating parent moments directly
count = 0;

for pmax = 1:10
    for mlevs = 1:mmax
        count = 0;
        nels = 8^(mlevs);
        for k1 = 1:pmax
            
            for k2 = 1:pmax-k1
                
                for k3 = 1:(pmax-(k1+k2))
                    
                    count = count + 1;
                end
            end
        end
        cstd(pmax,mlevs) = nels*count;
    end
end

%   For all branches at all layers do ISB - 189 in ISB
cstbr = zeros(10,mmax-1);
for maxp = 1:10
    for mlevs = 1:mmax-1
        count = 0;
        for ix = 1:2
            for iy = 1:2
                for iz = 1:2
                    for k1 = 1:maxp
                        
                        for k2 = 1:(maxp - k1)
                            
                            for k3 = 1:(maxp - (k1 + k2))
                                
                                for n1 = k1:maxp
                                    
                                    for n2 = k2:(maxp - n1)
                                       
                                        for n3 = k3:(maxp - (n1 + n2 ))
                                            
                                            count = count + 1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                    
                end
            end
        end
        disp([maxp mlevs 189*8^mlevs*count])
        cstbr(maxp,mlevs) = 189*8^mlevs*count;
    end
end
            
            
            
