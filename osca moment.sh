##osca
./osca_Mac --efile ln.txt --gene-expression --make-bod --no-fid --out lnprofile
#qc
osca --befile myprofile --sd-min 0.02 --make-bod --out newprofile
osca --befile myprofile --missing-ratio-probe 0.01 --make-bod --out newprofile
./osca_Mac --befile lnprofile --pheno agri1.phen --linear --out my
#mlm
./osca_Mac --moa-exact --befile ln2profile --pheno agri1.phen --out mymlmexactagri1
./osca_Mac --moa-exact --befile myprofile --pheno my.phen --task-num 1000 --task-id 1 --thread-num 10 --out my
for line in `cat ./phen`; do
    ./osca_Mac --moa-exact --befile lnprofile --pheno $line --out mymlmexact$line 
done
