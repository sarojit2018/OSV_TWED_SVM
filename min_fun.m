function f = min_fun(x)
traindata = cell(1,50);
trainX = cell(1,10);
forg_tr = cell(1,25);
gen_tr = forg_tr;

for i=0:0
        for j=0:0
            fname = sprintf('00%d%d',i,j);
            fname1 = [fname];
            dir1=dir(fname1);
            for k=3:52
                   filename = dir1(k).name;
                   filename1 = sprintf('00%d%d/%s',i,j,filename);
                   [x y z az in pps] = FPG_Signature_Read(filename1,0,0);
                   
                   traindata{i*10+j+1,k-2} = [x y z az in];
                   if(k<28)
                       forg_tr{i*10+j+1,k-2}=[x y z az in];
                   else
                       gen_tr{i*10+j+1,k-27}=[x y z az in];
                   end
            end
        end
end


%disp(gen_tr{1,1});

for i=1:25
    [a b] = size(gen_tr{1,i});
    mat = gen_tr{1,i};
    for j=2:a
       mat(j-1,:) = mat(j,:) - mat(j-1,:);
    end
    gen_tr{1,i} = mat;
end



for i=1:25
    [a b] = size(forg_tr{1,i});
    mat = forg_tr{1,i};
    for j=2:a
       mat(j-1,:) = mat(j,:) - mat(j-1,:);
    end
    forg_tr{1,i} = mat;
end

%disp(gen_tr{1,1});

for i=1:5
    trainX{1,i} = gen_tr{1,i};
end

for i=1:5
    trainX{1,5+i} = forg_tr{1,i};
end

trainY = [1 1 1 1 1 -1 -1 -1 -1 -1];

    val = 0;
    val1 = 0;
    for lv1=1:10
        val1 = val1 + x(lv1);
        for lv2=1:10
            mat1 = trainX{1,lv1};
            mat2 = trainX{1,lv2};
            val = val + 0.5*x(lv1)*x(lv2)*trainY(lv1)*trainY(lv2)*edrkernel(mat1,mat2);     
        end
    end
    f = val - val1;
end