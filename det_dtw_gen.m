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

eer_vals = zeros(100,1);
far_vals = zeros(100);
frr_vals = zeros(100);

 
far = [];
frr = [];

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
                   traindata{i*10+j+1,k-2} = [x y z az in tstamp];
                   if(k<28)
                       forg_tr{i*10+j+1,k-2}=[x y z az in tstamp];
                   else
                       gen_tr{i*10+j+1,k-27}=[x y z az in tstamp];
                   end
            end
        end 
end


%disp(gen_tr{1,1});
% max no. of points needed for each signature after interpolation

% min no of points in genuine signatures of the given user.


lp_var = 1;


targets = zeros(20*40,1);
outputs = zeros(20*40,1);


for u_num = 1:1
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
    
    [a b] = size(gen_mat);
    
    % calculating delta
    for j=2:a
       gen_mat(j-1,1:b-1) = gen_mat(j,1:b-1) - gen_mat(j-1,1:b-1);%calculate deltas of each feature x,y,p,in,az
       gen_mat(j-1,b) = gen_mat(j-1,b) - gen_mat(1,b);%calculate the time elapsed from timestamp zero
    end
    
    % last row is redundant so removed
    gen_mat(a,:) = [];
    tstamp = gen_mat(:,b);
    [st_size,dummy] = size(tstamp);
    tstamp(st_size,:) = [];
    gen_mat(:,b) = [];
    
    [a b] = size(gen_mat);
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
    
    %AR parameter estimation using yule-walker equations
    [gen_si_x,gen_si_y] = size(gen_mat_all);
    
    gen_AR_params = zeros(gen_si_x,gen_si_y);
   
    
    for loop_var1=1:gen_si_y
       feat_samples = gen_mat_all(:,loop_var1);
       feat_samples_row = reshape(gen_mat_all,[],1);
       AR_params = aryule(feat_samples_row,gen_si_x-1); %calculating AR parameters
       AR_params_col = reshape(AR_params,1,[]);
       gen_AR_params(:,loop_var1) = AR_params_col;
    end
    
     %z-normalization
  [a,num_feat] =  size(gen_mat_all);
  
  % mean_val contains mean for all features using all samples
  mean_val = zeros(num_feat);
  for j=1:num_feat
      mean_val(j) = mean(gen_mat_all(:,j));
  end
  
  % std_val contains std deviation for all features using all samples
  std_val = zeros(num_feat);
  for j=1:num_feat
      std_val(j) = std(gen_mat_all(:,j));
  end
  
  % undergoing normalization of ith genuine signature
  for j=1:num_feat%columns(features)
     for k=1:a%rows
         gen_mat_all(k,j) = (gen_mat_all(k,j) - mean_val(j))/std_val(j);
     end
  end
  
   
     %disp(gen_AR_params);
     
    %updation of global matrix    
    
    %disp(size(tstamp));
    %disp(size(gen_mat_all));
    gen_tr{u_num,i} = [gen_mat_all(:,1:5) gen_AR_params(:,1:5) tstamp];
    [dummy_x,dummy_y] = size(gen_tr{u_num,i});
    %setglobal_latent(dummy_y);1
    
end

 mdpoints = min_gen_data_points;
 %min data points are interpolated/resampled at these many points
 min_gen_data_points = floor(mdpoints);
 
 setglobalsize(min_gen_data_points);
 

disp(max_gen_data_points);
disp(min_gen_data_points);


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
    tstamp = forg_mat(:,b);
    [st_size,dummy] = size(tstamp);
    tstamp(st_size,:) = [];
    forg_mat(:,b) = [];
    
    [a b] = size(forg_mat);
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
     
    %AR parameter estimation using yule-walker equations
    [forg_si_x,forg_si_y] = size(forg_mat_all);
    
    forg_AR_params = zeros(forg_si_x,forg_si_y);
    
    for loop_var1=1:forg_si_y
       feat_samples = forg_mat_all(:,loop_var1);
       feat_samples_row = reshape(forg_mat_all,[],1);
       AR_params = aryule(feat_samples_row,forg_si_x-1);  %calculation of AR parameters
       AR_params_col = reshape(AR_params,1,[]);
       forg_AR_params(:,loop_var1) = AR_params_col;
    end
    
    %disp(gen_AR_params);
    
  %z-normalization
  [a,num_feat] =  size(forg_mat_all);
  
  %mean_val contains mean for all features using all samples
  mean_val = zeros(num_feat);
  for j=1:num_feat
      mean_val(j) = mean(forg_mat_all(:,j));
  end
  
  % std_val contains std deviation for all features using all samples
  std_val = zeros(num_feat);
  for j=1:num_feat
      std_val(j) = std(forg_mat_all(:,j));
  end
  
  % undergoing normalization of ith genuine signature
  for j=1:num_feat%columns(features)
     for k=1:a%rows
         forg_mat_all(k,j) = (forg_mat_all(k,j) - mean_val(j))/std_val(j);
     end
  end
    
    
    % updating the global matrix for ith forged signature of user
    %disp(size(forg_mat_all));
    %disp(size(tstamp));
    forg_tr{u_num,i} = [forg_mat_all(:,1:5) forg_AR_params(:,1:5) tstamp];
    [dummy_x,dummy_y] = size(forg_tr{u_num,i});
    %setglobal_latent(dummy_y);
end


% For Normalization
% end of Normalization

%fprintf('Sizes of genuine signatures before interpolation');
%for i=1:25
%    disp(size(gen_tr{1,i}));
%end

% For Interpolation/resampling of genuine signatures
% interpolation of genuine signatures
%fprintf('Sizes of genuine signatures after interpolation');
%for i=1:25
%    disp(size(gen_tr{1,i}));
%end

% For Interpolation
% interpolation/resampling of forged signatures
X = zeros(50,num_feat*min_gen_data_points);
Y = zeros(50,1);

for i=1:25
    Y(i) = 1;
    Y(25+i) = -1;
end
%the first few signatures are to be used as train

real_Y = [Y(6:50,:)];

average_threshold = 0;
for i=1:5
    for j = 1:5
       dtw_dist = dtw(gen_tr{u_num,i},gen_tr{u_num,j});
       average_threshold = max(average_threshold, dtw_dist);
    end
end


for th = average_threshold:0.05:(average_threshold+10)

testY = zeros(45,1);


for i=6:25
    avg_dtw_dist = 0.0;
    for j=1:5
        avg_dtw_dist = avg_dtw_dist + dtw(gen_tr{u_num,i},gen_tr{u_num,j});
    end
    avg_dtw_dist = avg_dtw_dist/5;
    
    if(avg_dtw_dist <= average_threshold)
        testY(i-5,1) = 1;
        
    else
        testY(i-5,1) = -1;
        FRR = FRR + 1;
    end
end


for i=1:25
    avg_dtw_dist = 0.0;
    for j=1:5
        avg_dtw_dist = avg_dtw_dist + dtw(forg_tr{u_num,i},gen_tr{u_num,j});
    end
    avg_dtw_dist = avg_dtw_dist/5;
    
    if(avg_dtw_dist <= average_threshold)
        testY(i+20,1) = 1;
        FAR = FAR + 1;
        
    else
        testY(i+20,1) = -1;
    end
end


disp("False Acceptance Rate");
disp(FAR*100/500);
disp("False Rejection Rate")
disp(FRR*100/400);
end
end