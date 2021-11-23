function d = dtw(X,Y)
    d = 0;
    p = 2;
    [si_x,k] = size(X);
    [si_y,k] = size(Y);
    %disp(size(X));
    %disp(size(Y));
    
    dtw_mat = zeros(si_x,si_y);
    dtw_mat(1,1) = 0;
    for i = 2:si_x
        dtw_mat(i,1) = realmax;
    end
    
    for i = 2:si_y
       dtw_mat(1,i) = realmax;
    end
    
    for i = 2:si_x
       for j = 2:si_y
          x = X(i,:);
          y = Y(j,:);
          
          dist = lpnorm(x,y,p);
          v = [dtw_mat(i-1,j-1),dtw_mat(i-1,j),dtw_mat(i,j-1)];
          min_val = min(v); 
          dtw_mat(i,j) = min_val + dist;
       end
    end
    
    dtw_len = 1;
    
    [si_x,k] = size(X);
    [si_y,k] = size(Y);
   
    
    i = si_x;
    j = si_y;
    
    %disp(i);
    %disp(j);    
        
    %disp("########");
    %disp(i);
    %disp(j);
    
   while(i~=1 || j~=1)
      if(i==1)
         dtw_len = dtw_len + 1;
         j = j - 1;
         
      elseif(j==1)
         dtw_len = dtw_len + 1;
         i = i - 1;
         
      else
           v = [dtw_mat(i-1,j),dtw_mat(i,j-1),dtw_mat(i-1,j-1)];
           mini = min(v);
           dtw_len = dtw_len + 1;
           if(mini==dtw_mat(i-1,j))
              i = i - 1;
           elseif(mini==dtw_mat(i,j-1))
              j = j - 1;
           else
               i = i - 1;
               j = j - 1;
      end
      end
   end
    %disp("Warping path Length");
    %disp(dtw_len);
    d = dtw_mat(si_x,si_y)/dtw_len;
end