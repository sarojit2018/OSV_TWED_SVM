function G = twedkernel(U,V)
    [si_u k] = size(U);
   [si_v k]= size(V);
   
   alpha = 0.5;
   
   %alpha = getglobalsize;
   
   beta = 2;
 
   G = zeros(si_u,si_v);
    
   nf = getglobal_latent;
   disp("#############");
   disp(nf);
   
   for i = 1:si_u
      for j = 1:si_v
          sig1 = U(i,:);
          sig2 = V(j,:);
          
          [z,size_sig] = size(sig1);
          %disp(size_sig);
          gen_sig1 = reshape(sig1,[size_sig/nf,nf]);
          gen_sig2 = reshape(sig2,[size_sig/nf,nf]);
          
          twedr = time_warp_edit_dist(gen_sig1,gen_sig2);
          
          val = -power(twedr,beta)*alpha;
          
          G(i,j) = exp(val);
          
      end
   end

end