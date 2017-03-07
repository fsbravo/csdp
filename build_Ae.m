function [Ae,invAeAeT,be] = build_Ae(n,BLKSZ)
    n_cons = BLKSZ*(BLKSZ+1)+1+2*BLKSZ*n;
    d = n*BLKSZ+1;
    Ae = zeros(n_cons,d^2);
    be = zeros(n_cons,1);
    l = 1;
    % bottom right corner is 1
    Ae(l,:) = [zeros(1,d^2-1) 1];
    be(l) = 1;
    l = l+1;
    % diagonal blocks: off-diagonal sums to 0
    for iii=1:BLKSZ
        idx = (iii-1)*n+1:iii*n;
        temp = zeros(d); temp(idx,idx) = ones(n)-eye(n);
        Ae(l,:) = temp(:)'/norm(temp,'fro');
        be(l) = 0;
        l = l+1;
    end
    % diagonal blocks: diagonal sums to 1
    for iii=1:BLKSZ
        idx = (iii-1)*n+1:iii*n;
        temp = zeros(d); temp(idx,d) = ones(n,1);
        Ae(l,:) = temp(:)'/norm(temp,'fro');
        be(l) = 1/norm(temp,'fro');
        l = l+1;
    end
    % off-diagonal blocks: diagonal sums to 0
    for iii=1:BLKSZ
        idx_i = (iii-1)*n+1:iii*n;
        for jjj=iii+1:BLKSZ
            idx_j = (jjj-1)*n+1:jjj*n;
            temp = zeros(d); temp(idx_i,idx_j) = eye(n); temp = temp+temp';
            Ae(l,:) = temp(:)'/norm(temp,'fro');
            be(l) = 0;
            l = l+1;
        end
    end
    % off-diagonal blocks: off-diagonal sums to 1
    for iii=1:BLKSZ
        idx_i = (iii-1)*n+1:iii*n;
        for jjj=iii+1:BLKSZ
            idx_j = (jjj-1)*n+1:jjj*n;
            temp = zeros(d); temp(idx_i,idx_j) = ones(n)-eye(n); temp = temp+temp';
            Ae(l,:) = temp(:)'/norm(temp,'fro');
            be(l) = 2/norm(temp,'fro');
            l = l+1;
        end
    end
    % diagonal elements are the same as u and v
    for iii=1:BLKSZ*n
        temp = zeros(d); temp(iii,iii) = 1; temp(iii,d) = -0.5; temp(d,iii) = -0.5;
        Ae(l,:) = temp(:)'/norm(temp,'fro');
        be(l) = 0;
        l = l+1;
    end 
    for iii=1:BLKSZ*n
        temp = zeros(d); temp(iii,d) = 1; temp(d,iii) = -1;
        Ae(l,:) = temp(:)'/norm(temp,'fro');
        be(l) = 0;
        l = l+1;
    end 
    Ae = sparse(Ae);
    invAeAeT = inv(Ae*Ae');
end