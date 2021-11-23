function twed = TWED(U,V)
    [si_u k] = size(U);
   [si_v k]= size(V);
   
   e = getglobal;
   alpha = 1;
   beta = 1;
 
    si = getglobalsize;
   G = zeros(si_u,si_v);
    
   for i = 1:si_u
      for j = 1:si_v
          sig1 = U(i,:);
          sig2 = V(j,:);
          
          [z,size_sig] = size(sig1);
          %disp(size_sig);
          gen_sig1 = reshape(sig1,[size_sig/5,5]);
          gen_sig2 = reshape(sig2,[size_sig/5,5]);
          
          edr = editdist(gen_sig1,gen_sig2,e);
          
          val = 2*si-edr;
          
          G(i,j) = beta*power(val,alpha);
          
      end
   end

end