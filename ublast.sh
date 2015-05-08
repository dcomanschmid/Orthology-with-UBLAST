####################################################################
# genome-wide ortholog search with "UBLAST"                        #
# UBLAST is implemented in USEARCH                                 #
#	- faster and more sensitive than BLAST                     #
# http://www.drive5.com/usearch/manual/install.html                #
#	1. Install usearch                                         #
# 	2. Run -ublast                                             #
# 	3. Find RHB (reciprocal best hits)                         #
#                                                                  #
# Diana Coman Schmid                                               #
# Eawag 2015                                                       #
# diana.comanschmid@eawag.ch                                       #
####################################################################


### 1. Install usearch    
# download the binary file from http://www.drive5.com/usearch/download.html 

cd /home/dianacs/Software/usearch/

# provide read/write access
chmod +x /home/dianacs/Software/usearch/usearch8.0.1623_i86linux32

# rename the binary to "usearch8" 
mv usearch8.0.1623_i86linux32 usearch8

# add (temporary) binary to PATH

export PATH="$PATH:/home/dianacs/Software/usearch/"

### Run UBLAST

cd /media/dianacs/data1/ortho_usearch/


# create UDB databases from the FASTA files 

usearch8 -makeudb_ublast /media/dianacs/data1/zebrafish/proteome/Danio_rerio.Zv9.pep.all.fa -output Danio_rerioZv9pepall.udb

usearch8 -makeudb_ublast /media/dianacs/data1/rainbowtrout/proteome/Oncorhynchus_mykiss_pep.fa -output Oncorhynchus_mykiss_pep.udb

# run UBLAST for zebrafish=fasta vs. rainbowtrout=udb 

usearch8 -ublast /media/dianacs/data1/zebrafish/proteome/Danio_rerio.Zv9.pep.all.fa -db Oncorhynchus_mykiss_pep.udb -evalue 1e-6 -blast6out ZT.u8

# run UBLAST for rainbowtrout=fasta vs. zebrafish=udb 

usearch8 -ublast /media/dianacs/data1/rainbowtrout/proteome/Oncorhynchus_mykiss_pep.fa -db Danio_rerioZv9pepall.udb -evalue 1e-6 -blast6out TZ.u8
