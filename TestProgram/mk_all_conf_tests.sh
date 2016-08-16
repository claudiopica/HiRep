#!/bin/bash
shopt -s nullglob

conf_dir="TestConfigurations-test"
make_dir="../Make"

flg_L1="${conf_dir}/MkFlags_L1_*"
flg_L2="${conf_dir}/MkFlags_L2_*"
flg_L3="${conf_dir}/MkFlags_L3_*"
flg_L4="${conf_dir}/MkFlags_L4_*"

confn=0;
for i in $flg_L1
do
  for j in $flg_L2
  do
    for k in $flg_L3
    do
      for l in $flg_L4
      do
        ((confn++))
        echo "Testing Configuration #$confn"
	cat $i $j $k $l > ${make_dir}/MkFlags
        make cleanall >/dev/null
        make -j4 tests >mk_log 2>/dev/null
        echo "Running tests"
      done
    done
  done
done

shopt -u nullglob



