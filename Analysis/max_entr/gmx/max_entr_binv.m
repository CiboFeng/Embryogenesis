function max_entr_binv(bfile,binvfile)
    data=load(bfile);
    B=data.B;
    B_=pinv(B);
    save(binvfile,'B_');
end