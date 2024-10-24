
# LFQ_timsTOFPro_diaPASEF_Ecoli_03.d.zip

docker run --platform linux/amd64 -v ${PWD}:/data \
  ghcr.io/pride-archive/aspera \
  ascp -i /home/aspera/.aspera/cli/etc/asperaweb_id_dsa.openssh \
  -TQ -P33001 \
  prd_ascp@fasp.ebi.ac.uk:/pride/data/archive/2022/02/PXD028735/LFQ_timsTOFPro_diaPASEF_Condition_A_Sample_Alpha_02.d.zip  \
  /data
