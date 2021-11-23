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
epsilon = 0.0001;

eer_vals = zeros(20,1);
far_vals = zeros(100);
frr_vals = zeros(100);


dtw_values_gen_gen = zeros(100,25);
twed_values_gen_gen = zeros(100,25);


dtw_values_for_gen = zeros(100,25);
twed_values_for_gen = zeros(100,25);

diff_dtw_avg = zeros(100);
diff_twed_avg = zeros(100);


for i=0:9
        for j=0:9
            fname = sprintf('00%d%d',i,j);
            fname1 = [fname];
            dir1=dir(fname1);
            for k=3:52
                   filename = dir1(k).name;
                   filename1 = sprintf('00%d%d/%s',i,j,filename);
                   [x y z az in pps] = FPG_Signature_Read(filename1,0,0);
                   
                   si = size(x); %size of each vector
                   tstamp = zeros(si); %timestamp
                   for p=1:si
                        tstamp(p) = p/pps;
                   end
                   traindata{i*10+j+1,k-2} = [x y z az in];
                   if(k<28)
                       forg_tr{i*10+j+1,k-2}=[x y z az in];
                   else
                       gen_tr{i*10+j+1,k-27}=[x y z az in];
                   end
            end
        end
end



for u_num = 1:100
    max_gen_data_points = 0; 
    min_gen_data_points = 99999;   
for i=1:25
    [a b] = size(gen_tr{u_num,i});
    
    if(a > max_gen_data_points) 
        max_gen_data_points = a;
    end
    if(a < min_gen_data_points)
        min_gen_data_points = a;
    end
    
    %min_gen_data_points = int8(min_gen_data_points);

    
    % matrix storing values for a ith genuine signature
    gen_mat = gen_tr{u_num,i};  
    
    % calculating delta
    for j=2:a
       gen_mat(j-1,1:b-1) = gen_mat(j,1:b-1) - gen_mat(j-1,1:b-1);%calculate deltas of each feature x,y,p,in,az
       gen_mat(j-1,b) = gen_mat(j-1,b) - gen_mat(1,b);%calculate the time elapsed from timestamp zero
    end
    
    % last row is redundant so removed
    gen_mat(a,:) = [];
    
    [a b] = size(gen_mat);
    
    delta2_x = zeros(a,1); %secondary delat_x values
    delta2_y = zeros(a,1); %secondary delta_y values
    
    for j=2:a
       delta2_x(j-1,1) = gen_mat(j,1) - gen_mat(j-1,1);
       delta2_y(j-1,1) = gen_mat(j,2) - gen_mat(j-1,2);
    end
    
    new_gen_mat = [gen_mat delta2_x delta2_y];%contains secondary deltas along with the primary deltas
    
    new_gen_mat(a,:) = []; %deleting the redundant last row
    
     
    [a b] = size(new_gen_mat);
    
    len_delta1 = zeros(a,1); %primary delta length
    len_delta2 = zeros(a,1); %secondary delta length
   
   
    for j = 1:a
        len_delta1(j,1) = sqrt(new_gen_mat(j,1)*new_gen_mat(j,1) + new_gen_mat(j,2)*new_gen_mat(j,2) + epsilon);
        len_delta2(j,1) = sqrt(new_gen_mat(j,6)*new_gen_mat(j,6) + new_gen_mat(j,7)*new_gen_mat(j,7) + epsilon);
    end
    
    sin_alpha = zeros(a,1);
    cos_alpha = zeros(a,1);
    
    for j =1:a 
       sin_alpha(j,1) = new_gen_mat(j,2)/len_delta1(j,1);
       cos_alpha(j,1) = new_gen_mat(j,1)/len_delta1(j,1);
    end
    
    gen_mat_delta1_delta2 = new_gen_mat;
    gen_mat_delta1_only = gen_mat;
    gen_mat_all = [new_gen_mat sin_alpha cos_alpha len_delta1 len_delta2];
    % updation of global matrix 
    gen_tr{u_num,i} = gen_mat_all(:,1:5);
end

 mdpoints = min_gen_data_points;
 %min data points are interpolated/resampled at these many points
 min_gen_data_points = floor(mdpoints);
 
 setglobalsize(min_gen_data_points);


for i=1:25
    [a b] = size(forg_tr{u_num,i});
    
    % matrix storing values for ith forged signature
    forg_mat = forg_tr{u_num,i};
    
    % calculating delta for the forged signature
    for j=2:a
       forg_mat(j-1,1:b-1) = forg_mat(j,1:b-1) - forg_mat(j-1,1:b-1);
       forg_mat(j-1,b) = forg_mat(j-1,b) - forg_mat(1,b);
    end
    % last row is redundant so removed
    forg_mat(a,:) = [];
    
    
    [a b] = size(forg_mat);
    
    delta2_x = zeros(a,1); %secondary delat_x values
    delta2_y = zeros(a,1); %secondary delta_y values
    
    for j=2:a
       delta2_x(j-1,1) = forg_mat(j,1) - forg_mat(j-1,1);
       delta2_y(j-1,1) = forg_mat(j,2) - forg_mat(j-1,2);
    end
    
    new_forg_mat = [forg_mat delta2_x delta2_y];%contains secondary deltas along with the primary deltas
    
    new_forg_mat(a,:) = []; %deleting the redundant last row
    
     
    [a b] = size(new_forg_mat);
    
    len_delta1 = zeros(a,1); %primary delta length
    len_delta2 = zeros(a,1); %secondary delta length
   
   
    for j = 1:a
        len_delta1(j,1) = sqrt(new_forg_mat(j,1)*new_forg_mat(j,1) + new_forg_mat(j,2)*new_forg_mat(j,2) + epsilon);
        len_delta2(j,1) = sqrt(new_forg_mat(j,6)*new_forg_mat(j,6) + new_forg_mat(j,7)*new_forg_mat(j,7) + epsilon);
    end
    
    sin_alpha = zeros(a,1);
    cos_alpha = zeros(a,1);
    
    for j =1:a 
       sin_alpha(j,1) = new_forg_mat(j,2)/len_delta1(j,1);
       cos_alpha(j,1) = new_forg_mat(j,1)/len_delta1(j,1);
    end
    
    forg_mat_delta1_delta2 = new_forg_mat;
    forg_mat_delta1_only = forg_mat;
    forg_mat_all = [new_forg_mat sin_alpha cos_alpha len_delta1 len_delta2];
    
    % updating the global matrix for ith forged signature of user
    forg_tr{u_num,i} = forg_mat_all(:,1:5);
end

% For Normalization
% calculating mean and standard deviation for genuine signatures
for i=1:25
  gen_mat = gen_tr{u_num,i};
  [a,num_feat] =  size(gen_mat);
  
  % mean_val contains mean for all features using all samples
  mean_val = zeros(num_feat);
  for j=1:num_feat
      mean_val(j) = mean(gen_mat(:,j));
  end
  
  % std_val contains std deviation for all features using all samples
  std_val = zeros(num_feat);
  for j=1:num_feat
      std_val(j) = std(gen_mat(:,j));
  end
  
  % undergoing normalization of ith genuine signature
  for j=1:num_feat%columns(features)
     for k=1:a%rows
         gen_mat(k,j) = (gen_mat(k,j) - mean_val(j))/std_val(j);
     end
  end
  
  % updating the global matrix with normalized values
  gen_tr{u_num,i} = gen_mat;
end

% For Normalization
% calculating mean and standard deviation for forged signatures
for i=1:25
  forg_mat = forg_tr{u_num,i};
  [a,num_feat] =  size(forg_mat);
  
  % mean_val contains mean for all features using all samples
  mean_val = zeros(num_feat);
  for j=1:num_feat
      mean_val(j) = mean(forg_mat(:,j));
  end
  
  % std_val contains std deviation for all features using all samples
  std_val = zeros(num_feat);
  for j=1:num_feat
      std_val(j) = std(forg_mat(:,j));
  end
  
  % undergoing n    ormalization of ith forged signature
  for j=1:num_feat%columns(features)
     for k=1:a %rows
         forg_mat(k,j) = (forg_mat(k,j) - mean_val(j))/std_val(j);
     end
  end
  
  % updating the global matrix with normalized values
  forg_tr{u_num,i} = forg_mat;
end

%fprintf('Sizes of genuine signatures before interpolation');
%for i=1:25
%    disp(size(gen_tr{1,i}));
%end


% For Interpolation/resampling of genuine signatures
% interpolation of genuine signatures
for i=1:25
    gen_mat = gen_tr{u_num,i};
    
    [gen_size,num_feat] = size(gen_mat);
    
    new_gen_mat = zeros(min_gen_data_points,num_feat);
    
    for lv=1:num_feat
        x = gen_mat(:,lv);
        y = resample(x,min_gen_data_points,gen_size);
        new_gen_mat(:,lv) = y;
    end
    
    gen_tr{u_num,i} = new_gen_mat;    
end

%fprintf('Sizes of genuine signatures after interpolation');
%for i=1:25
%    disp(size(gen_tr{1,i}));
%end

% For Interpolation
% interpolation/resampling of forged signatures
for i=1:25
    forged_mat = forg_tr{u_num,i};
    
    %cur_size -> no of sample points presently in ith sign & k -> features
    [forg_size,num_feat] = size(forged_mat);
   
    new_forg_mat = zeros(min_gen_data_points,num_feat);
    
    for lv = 1:num_feat
        x = gen_mat(:,lv);
        y = resample(x,min_gen_data_points,gen_size);
        new_forg_mat(:,lv) = y;
    end  
    
    forg_tr{u_num,i} = new_forg_mat;
end

cnt = 1;
for i=1:5
    for j=1:5
        dtw_values_gen_gen(u_num,cnt) = dtw(gen_tr{u_num,i},gen_tr{u_num,j});
        dtw_values_for_gen(u_num,cnt) = dtw(forg_tr{u_num,i},gen_tr{u_num,j});
        twed_values_gen_gen(u_num,cnt) = time_warp_edit_dist(gen_tr{u_num,i},gen_tr{u_num,j});
        twed_values_for_gen(u_num,cnt) = time_warp_edit_dist(forg_tr{u_num,i},gen_tr{u_num,j});
         if(i==j) 
           dtw_values_gen_gen(u_num,cnt) = 0;
           twed_values_gen_gen(u_num,cnt) = 0;
         end
        cnt = cnt+1;
    end
end

gen_gen_avg_dist_dtw = 0;
for_gen_avg_dist_dtw = 0;
gen_gen_avg_dist_twed = 0;
for_gen_avg_dist_twed = 0;

gen_gen_avg_dist_dtw = sum(dtw_values_gen_gen(u_num,:))/20;
for_gen_avg_dist_dtw = sum(dtw_values_for_gen(u_num,:))/25;

gen_gen_avg_dist_twed = sum(twed_values_gen_gen(u_num,:))/20;
for_gen_avg_dist_twed = sum(twed_values_for_gen(u_num,:))/25;



disp("Difference between gen-gen and for-gen average for DTW");
disp(for_gen_avg_dist_dtw - gen_gen_avg_dist_dtw);


disp("Difference between gen-gen and for-gen average for TWED");
disp(for_gen_avg_dist_twed - gen_gen_avg_dist_twed);


diff_dtw_avg(u_num) = abs(for_gen_avg_dist_dtw - gen_gen_avg_dist_dtw);
diff_twed_avg(u_num) = abs(for_gen_avg_dist_twed - gen_gen_avg_dist_twed);
end

users = [1:100];
plot(users,diff_dtw_avg);
hold on;
plot(users,diff_twed_avg);
hold off;

