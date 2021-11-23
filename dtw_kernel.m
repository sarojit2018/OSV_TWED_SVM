function G = dtw_kernel(U,V)
    [si_u k] = size(U);
   [si_v k]= size(V);
   
   alpha = 0.5;
   beta = 2;
 
   nf = getglobal_latent;
   G = zeros(si_u,si_v);
    
   for i = 1:si_u
      for j = 1:si_v
          sig1 = U(i,:);
          sig2 = V(j,:);
          
          [z,size_sig] = size(sig1);
          %disp(size_sig);
          gen_sig1 = reshape(sig1,[size_sig/nf,nf]);
          gen_sig2 = reshape(sig2,[size_sig/nf,nf]);
          
          
          gen_sig1(:,nf) = [];
          gen_sig2(:,nf) = [];
          
           %disp(size(gen_sig1));
          %disp(size(gen_sig2));
          
          dtw_dist = dtw(gen_sig1,gen_sig2);
          
          val = -power(dtw_dist,beta)*alpha;
          
          G(i,j) = exp(val);
          
      end
   end

end