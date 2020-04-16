from astropy.table import Table
import os
dir='full-data'
flist=os.listdir(dir)

for f in flist:
    if f.find('.csv')>0:
        fileIn=os.path.join(dir,f)
        fileOut=os.path.join(dir,f.replace('.csv','.fits'))
        print(fileIn,fileOut)
        data=Table.read(fileIn)
        data.write(fileOut)

