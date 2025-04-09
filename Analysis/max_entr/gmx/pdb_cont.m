function [P,C]=pdb_cont(i,file,d0,hipair,beta)
    now=datetime('now');
    fprintf('Begin to read .pdb files: %s\n',datestr(now));

    pdb=fopen(sprintf('%s%d.pdb',file,i),'r');
    lspdb = textscan(pdb, '%s', 'Delimiter', '\n');
    lspdb = lspdb{1};
    fclose(pdb);
    
    xx = [];
    yy = [];
    zz = [];
    natom = [];
    t = [];
    
    % Parse the lines of the PDB file
    for i = 1:length(lspdb)
        l_pdb = strsplit(strtrim(lspdb{i}));
        if strcmp(l_pdb{1}, 'ATOM')
            xx(end+1) = str2double(l_pdb{6});
            yy(end+1) = str2double(l_pdb{7});
            zz(end+1) = str2double(l_pdb{8});
            natom(end+1) = str2double(l_pdb{2});
        elseif strcmp(l_pdb{1}, 'TITLE')
            t(end+1) = str2double(l_pdb{end});
        end
    end
    
    natom = max(natom);
    nsample = length(t);
    
    x = zeros(nsample, natom);
    y = zeros(nsample, natom);
    z = zeros(nsample, natom);
    
    % Reshape the coordinates for each sample
    for i = 1:nsample
        for j = 1:natom
            x(i, j) = xx((i - 1) * natom + j);
            y(i, j) = yy((i - 1) * natom + j);
            z(i, j) = zz((i - 1) * natom + j);
        end
    end

    [nfr, N] = size(x);
    nhipair=size(hipair,1);
    r = zeros(nfr,N,3);
    r(:,:,1) = x;
    r(:,:,2) = y;
    r(:,:,3) = z;
    if ~isempty(hipair)
        P = zeros(1,nhipair);
        C = zeros(nhipair,nhipair);
        for i = 1:nfr
            D = squareform(pdist(squeeze(r(i,:,:)),'euclidean'));
            PD = 0.5*(1+tanh(beta*(D-d0)));
            row = hipair(:,1);
            col = hipair(:,2);
            P_temp = PD(sub2ind(size(PD),row,col));
            P = P + P_temp;
            C = C + P_temp' * P_temp;
        end
    else
        Pv = zeros(1,N^2);
        C = zeros(N^2,N^2);
        for i = 1:nfr
            D = squareform(pdist(squeeze(r(i,:,:)),'euclidean'));
            PD = 0.5*(1+tanh(beta*(D-d0)));
            P_temp = zeros(1,N^2);
            for j = 1:N^2
                [n, m] = ind2sub([N,N],j);
                P_temp(j) = PD(n,m);
            end
            Pv = Pv + P_temp;
            C = C + P_temp' * P_temp;
        end
        P=reshape(Pv,N,N);
    end
    P = P / nfr;
    C = C / nfr;
    disp('Complete P and C');

end