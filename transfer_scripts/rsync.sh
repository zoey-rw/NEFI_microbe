# script for rsync between pecan2 and scc data directories.

# sending ITS data from SCC to pecan
rsync -av -O -e ssh /projectnb/talbot-lab-data/NEFI_data/ITS/scc_gen/ zrwerbin@pecan2.bu.edu:/fs/data3/caverill/NEFI_data/ITS/scc_gen/

# sending ITS data from pecan to SCC
rsync -av -O zrwerbin@pecan2.bu.edu:/fs/data3/caverill/NEFI_data/ITS/pecan_gen/ /projectnb/talbot-lab-data/NEFI_data/ITS/pecan_gen/

# sending 16S data from SCC to pecan
rsync -av -O -e ssh /projectnb/talbot-lab-data/NEFI_data/16S/scc_gen/ zrwerbin@pecan2.bu.edu:/fs/data3/caverill/NEFI_data/16S/scc_gen/

# sending 16S data from pecan to SCC
rsync -av -O zrwerbin@pecan2.bu.edu:/fs/data3/caverill/NEFI_data/16S/pecan_gen/ /projectnb/talbot-lab-data/NEFI_data/16S/pecan_gen/

# set permissions to read/write/execute for everyone.
chmod 777 -R /projectnb/talbot-lab-data/NEFI_data/