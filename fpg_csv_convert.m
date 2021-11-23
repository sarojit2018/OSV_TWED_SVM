addpath(genpath('C:\Users\Sarojit Auddya\mcyt\mcyt\bt'));

data_path = '';
MyPath = userpath;
MyDir = MyPath(1:strfind(MyPath,';')-1);
MyWorkDir = genpath(MyDir);
addpath(MyWorkDir, '-end');

traindata = cell(100,50);

forg_tr = cell(100,25);
gen_tr = forg_tr;

FAR = 0;
FRR = 0;


eer_vals = zeros(20,1);
far_vals = zeros(100);
frr_vals = zeros(100);

for i=0:9
        for j=0:9
            fname = sprintf('00%d%d',i,j);
            fname1 = [fname];
            dir1 = dir(fname1);
            for k=3:52
                   filename = dir1(k).name;
                   filename1 = sprintf('00%d%d/%s',i,j,filename);
                   [x y z az in pps] = FPG_Signature_Read(filename1,0,0);
                   
                   si = size(x); %size of each vector
                   tstamp = zeros(si); %timestamp
                   for p=1:si
                        tstamp(p) = p/pps;
                   end
                   traindata{i*10+j+1,k-2} = [x y z az in tstamp];
                   if(k<28)
                       forg_tr{i*10+j+1,k-2}=[x y z az in tstamp];
                       [si_x,si_y] = size(forg_tr{i*10+j+1,k-2});
                       mat = forg_tr{i*10+j+1,k-2};
                       str = sprintf('00%d%d_mod',i,j);
                       %mkdir str;
                       filename_mod = sprintf('00%d%df%d.csv',i,j,k-3);
                       str_file = sprintf('00%d%d_mod/%s',i,j,filename_mod);
                       %writematrix(mat,str_file);
                       %fid = fopen(str_file,'wt');
                       csvwrite(filename_mod,mat);
                       %fclose(fid);
                   else
                       gen_tr{i*10+j+1,k-27}=[x y z az in tstamp];
                       %[si_x,si_y] = size(gen_tr{i*10+j+1,k-2});
                       mat = gen_tr{i*10+j+1,k-27};
                       str = sprintf('00%d%d_mod',i,j);
                       %mkdir str;
                       filename_mod = sprintf('00%d%dg%d.csv',i,j,k-3);
                       str_file = sprintf('00%d%d_mod/%s',i,j,filename_mod);
                       %writematrix(mat,str_file);
                       csvwrite(filename_mod,mat);
                   end
            end
        end
end



