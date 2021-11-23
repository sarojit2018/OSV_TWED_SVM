path_str = 'C:\Users\Sarojit Auddya\Documents\Audacity\160102077_';

vowels= ['a','e','i','o','u'];

for i=1:5
   file_path_prefix = vowels(i)+'\'+'160102077_'+vowels(i)+'_';
   template_file = file_path_prefix + '1'+'.wav';
   [y_tp,fs] = audioread(template_file);
   [tp_coeff,tp_delta1,tp_delta2,tp_loc] = mfcc(y_tp,fs);
   tp_f1 = tp_coeff(:,2:end);
   tp_f2 = tp_delta1(:,2:end);
   tp_f3 = tp_delta2(:,2:end);
   feat = [tp_f1,tp_f2,tp_f3];
   tp_feat = transpose(feat);
   dtw_scores = zeros(24);
   for j=2:25
       file_path = file_path_prefix + num2str(j)+'.wav';
       [y,fs] = audioread(file_path);
       [coeff,delta1,delta2,loc] = mfcc(y,fs);
       f1 = coeff(:,2:end);
       f2 = delta1(:,2:end);
       f3 = delta2(:,2:end);
       f = [f1,f2,f3];
       feat = transpose(f);
       [dtw_scores(j-1),x,y] = dtw(feat,tp_feat);
       plot(x,y);
       hold on;
   end
   hold off;
end