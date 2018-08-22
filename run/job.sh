#!/bin/bash

source ../env.sh

echo "$@"
echo "pwd -P: $(pwd -P)"
echo "TMPDIR = $TMPDIR"
echo "LD_LIBRARY_PATH = $LD_LIBRARY_PATH"

# files=$(python -c "
# import os, json
# with open('$1','r+') as f:
#     j = json.load(f)
#     input = j['input'][0]
#     file = ['files'][0]
#     print file
#     input['original_files'] = [ file ]
#     file = '$TMPDIR/'+os.path.basename(file)
#     print file
#     json.dump(j,f)
# ")
# # echo "ifname = $ifname"
#
# # echo ${files[0]}
# # echo ${files[1]}
# echo ${files}
#
# exit
#
# # flag_file=~/.condor_tmp_copy
# # while [ -a $flag_file ]; do sleep 1; done
# # touch $flag_file
# # tmp_file=$(cp -v $ifname $TMPDIR/)
# # echo $tmp_file
# # tmp_file=$(sed "s/[\`\']//g;s/.*-> *//" <<< $tmp_file)
# # rm $flag_file
#
# flag_file=~/.condor_tmp_copy
# while [ -a $flag_file ]; do sleep 1; done
# touch $flag_file
# cp -v ${files[0]} $TMPDIR/
# rm $flag_file
#
# # path='/home/ivanp/work/ntuple_analysis'
# # ${path}/bin/analyses/test ${path}/runcards/test.json $1
#
# # rm $tmp_file
# rm ${files[1]}

path='/home/ivanp/work/ntuple_analysis'
${path}/bin/analyses/test ${path}/runcards/test.json $1 \
  --tmp-dir=$TMPDIR

