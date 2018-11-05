rsync -av -O -e ssh /projectnb/talbot-lab-data/NEFI_data/ITS/scc_gen/ zrwerbin@pecan2.bu.edu:/fs/data3/caverill/NEFI_data/ITS/scc_gen/

rsync -av -O  --no-p zrwerbin@pecan2.bu.edu:/fs/data3/caverill/NEFI_data/ITS/pecan_gen/ /projectnb/talbot-lab-data/NEFI_data/ITS/pecan_gen/

rsync -av -O -e ssh /projectnb/talbot-lab-data/NEFI_data/16S/scc_gen/ zrwerbin@pecan2.bu.edu:/fs/data3/caverill/NEFI_data/16S/scc_gen/

rsync -av -O --no-p zrwerbin@pecan2.bu.edu:/fs/data3/caverill/NEFI_data/16S/pecan_gen/ /projectnb/talbot-lab-data/NEFI_data/16S/pecan_gen/