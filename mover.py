import os
import shutil
path = 'C:/Users/Sarojit Auddya/mcyt/mcyt/'
#os.mkdir('tempdir')

for i in range(10):
	for j in range(10):
		foldername = "00" + str(i) + str(j)
		filename = ""
		dest = foldername + "_mod"
		os.mkdir(dest)
		for k in range(50):
			if(k<25):
				filename = foldername+ "f" + str(k) + ".csv"
			else:
				filename = foldername+ "g" + str(k) + ".csv"
			shutil.move(filename,dest)





