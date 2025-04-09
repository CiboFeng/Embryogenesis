function max_entr_cont(file,r0,beta,tf1,tf2)
    r0=str2double(r0);
    beta=str2double(beta);
    data=load(tf1);
    idx=data.idx;
    hipair=data.hipair+1;
    nhipair=size(hipair,1);
    
    parpool(length(idx));
    P = zeros(1,nhipair);
    C = zeros(nhipair,nhipair);
    parfor i = 1:length(idx)
        [P_temp,C_temp] = pdb_cont(idx(i),file,r0,hipair,beta);
        P=P+P_temp
        C=C+C_temp
    end
    P=P/length(idx);
    C=C/length(idx);
    
    save(tf2, 'P', 'C','-v7.3');
    now = datetime('now');
    fprintf('Complete P and C %s\n', datestr(now));
end