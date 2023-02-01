#emmax_in.tped/emmax_in.tfam/emmax_in.nosex
#make kinship matrix
/public/home/leimengyu/emmax-beta-07Mar2010/emmax-kin /public/home/leimengyu/emmax-beta-07Mar2010/emmax_in -v -d 10
#GWAS
/public/home/leimengyu/emmax-beta-07Mar2010/emmax -t /public/home/leimengyu/emmax-beta-07Mar2010/emmax_in -o /public/home/leimengyu/emmax-beta-07Mar2010/LL_emmax -p /public/home/leimengyu/emmax-beta-07Mar2010/LL.txt -k /public/home/leimengyu/emmax-beta-07Mar2010/emmax_in.BN.kinf
