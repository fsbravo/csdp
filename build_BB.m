function [BB_col] = build_BB(n,BLKSZ)
    BB_col = cell(1,BLKSZ);
    d = BLKSZ*n+1;
    for j=1:BLKSZ
        BB = spalloc(n,d^2,2*n);
        for iii=1:n
            temp = zeros(d); temp(iii+(j-1)*n,d) = 1/2; temp(d,iii+(j-1)*n) = 1/2;
            BB(iii,:) = temp(:)';
        end
        BB_col{j} = BB;
    end
end