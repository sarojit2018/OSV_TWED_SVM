function twed = time_warp_edit_dist(X,Y)
    twed = 0;
    
    prop = 0.1;%deletion cost(lambda) proportionality cost
    %prop = getgloballambda;
    
    disp(prop);
    
    [si_x,k] = size(X);
    [si_y,k] = size(Y);
    
    %disp("#########****************");
    %disp(k);
    
    lambda = prop*k;
    gamma = 0.00000001; %stiffness parameter(=0 means distance only on features and no temporal element)
    p = 2;
    twedmat = zeros(si_x,si_y);
    
    twedmat(1,1) = 0;
    
    tstamp_x = X(:,k);
    tstamp_y = Y(:,k);
    
    X(:,k) = [];
    Y(:,k) = [];
    
    %disp(size(X));
    %disp(size(Y));
    
    for i=2:si_x
        twedmat(i,1) = realmax; 
    end
    
    for i=2:si_y
        twedmat(1,i) = realmax; 
    end
    
    for i=2:si_x
       for j=2:si_y
           x_prev = X(i-1,:);
           y_prev = Y(j-1,:);
           x_cur = X(i,:);
           y_cur = Y(j,:);
           
           tx_prev = tstamp_x(i-1);
           ty_prev = tstamp_y(j-1);
           
           tx_cur = tstamp_x(i);
           ty_cur = tstamp_y(j);
           
           
           
           val1 = twedmat(i-1,j) + lpnorm(x_prev,x_cur,p) + gamma*lpnorm(tx_cur,tx_prev,p) +  lambda;
           val2 = twedmat(i-1,j-1) + (lpnorm(x_prev,y_prev,p) + gamma*lpnorm(tx_prev,ty_prev,p) + gamma*lpnorm(tx_cur,ty_cur,p) + lpnorm(x_cur,y_cur,p));
           val3 = twedmat(i,j-1) + lpnorm(y_prev,y_cur,p) + gamma*lpnorm(ty_cur,ty_prev,p) +lambda;
           
           v = [val1,val2,val3];
           
           twedmat(i,j) = min(v);
       end
    end
    
    twed_len = 1;
    i = si_x;
    j = si_y;
   while(i>1 || j>1)
      if(i==1)
         twed_len = twed_len + 1;
         j = j - 1;
         
      elseif(j==1)
         twed_len = twed_len + 1;
         i = i - 1;
         
      else
          %minim = realmax;
           v = [twedmat(i-1,j),twedmat(i,j-1),twedmat(i-1,j-1)];
           minim = min(v)
           twed_len = twed_len + 1;
           if(minim==twedmat(i-1,j))
              i = i - 1;
           elseif(minim==twedmat(i,j-1))
              j = j - 1;
           else
               i = i - 1;
               j = j - 1;
      end
      end
   end
   
   disp(minim);
    
   %disp("xxxxxxxxxx");
   %disp(twed_len);
    
    twed = twedmat(si_x,si_y)/twed_len;
end