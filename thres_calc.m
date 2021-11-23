function thres = thres_calc(gen_X)
    thres = 0;
    local_maxima = 0;

    si = size(gen_X);
    for i=1:si
        for j = i+1:si
            sig1 = gen_X(i,:);
            sig2 = gen_X(j,:);
            [d,s] = size(sig1);
            %disp(s);
            a = s/5;
            gen_sig1 = reshape(sig1,[a,5]);
            gen_sig2 = reshape(sig2,[a,5]);
            for k=1:a
               local_minima = realmax;
               for l=1:a
                    sample_k = gen_sig1(k,:);
                    sample_l = gen_sig2(l,:);
                    k_l_dist = dist(sample_k,sample_l);
                    if(local_minima > k_l_dist)
                        local_minima = k_l_dist;
                    end
               end
               
               if(local_maxima < local_minima) 
                  local_maxima = local_minima; 
               end
            end
            
            if(thres < local_maxima)
                thres = local_maxima;
            end
        end
    end

end