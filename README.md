# SCNS
Software and computing for nuclear and subnuclear physics

#How to run the Code:
The machines of the CNAF have the possibility of installing CAFana on their files following the 
First time you set up a new version of sbnana:
1. source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh newmrb
2. cd /storage/gpfs_data/icarus/local/users/USERNAME/
3. ups list -aK+ sbnana <– List all the existing versions of SBNAna, choose one with “e20:prof” tag
4. mkdir sbnana_v09_78_06
5. cd sbnana_v09_78_06
6. setup sbnana v09_78_06 -q e20:prof
7. export MRB_PROJECT=sbnana
8. mrb newDev
9. source /storage/gpfs_data/icarus/local/users/USERNAME/sbnana_v09_78_06/localProducts_sbnana_v09_78_06_e20_prof/setup
