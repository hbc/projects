import os.path
import urllib
import gzip
import csv
import requests

dataDir = '/Volumes/ody/consults/lk_FOY/data/'
dataFH = '/Volumes/ody/consults/lk_FOY/meta/unified.metadata.refined.PBMC.tab'

with open(dataFH, 'rb') as tabfile:
    # read tab delimited file in with csv library
    inf = csv.reader(tabfile, delimiter="\t")
    # skip ehader
    next(inf)
    # load remaining lines into a list
    lol=list(inf)

# iterate through rows of list 
for row in lol:
    FTPloc= row[5]
    writeDir = row[8]
    #grab the name of the gzipped CEL file from the ftp url
    gzcelFH = os.path.basename(FTPloc)
    print(gzcelFH)
    #split off the gz file extension to get the uncompressed CEL file name
    celFH = os.path.splitext(gzcelFH)[0]
    
    # if no FTP url, goto next line
    if FTPloc == 'NA':
        continue
    # if an uncompressed CEL file already exists at the correct location , goto next line
    if os.path.isfile(os.path.join(dataDir, writeDir, celFH)):
        continue

    # download the gzipped CEL file from the FTP url and output to correct location and filename
    urllib.urlretrieve(FTPloc, os.path.join(dataDir, writeDir, gzcelFH))
    
    # load in the gzipped CEL file and uncompress
    gzippedCELdata = gzip.open(os.path.join(dataDir, writeDir, gzcelFH), 'rb')
    uncompressedCELdata = gzippedCELdata.read()
        
    # write out unzipped CEL file
    with open(os.path.join(dataDir, writeDir, celFH), 'w') as file:
        file.write(uncompressedCELdata)

    # remove gzipped CEL file
    os.remove(os.path.join(dataDir, writeDir, gzcelFH))
    print(celFH)

print "done"