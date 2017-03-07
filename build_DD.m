function [DD] = build_DD(n,BLKSZ)
    d = BLKSZ*n+1;
    DD = spalloc((BLKSZ*n)^2,d^2,(BLKSZ*n)^2);
    l = 1;
    for iii=1:BLKSZ*n
        for jjj=1:BLKSZ*n
            DD(l,(jjj-1)*d+iii) = 1;
            l = l+1;
        end
    end
end