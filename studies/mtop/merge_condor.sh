#!/bin/bash

dir=out

for x in $(
  find $dir/ -name '*.root' |
  sed 's:.*/\(.*\)_[0-9]\+\.root$:\1:' |
  sort -u
); do
  echo $x

  cat > .merge_${x}.sh << SCRIPT
#!/bin/bash
export PATH=${PATH}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}

cd $(pwd -P)
hadd .merge_${x}.root $dir/${x}_*.root
/home/ivanp/work/ntuple_analysis/bin/merge_root ${x}.root .merge_${x}.root

if [ ! -s ".merge_${x}.err" ]; then # if empty
  rm -f .merge_${x}.err .merge_${x}.out
fi
rm -f .merge_${x}.root
rm -- \$0
SCRIPT
  chmod +x .merge_${x}.sh

  condor_submit - > /dev/null << JOB
Universe   = vanilla
Executable = $(pwd -P)/.merge_${x}.sh
Output     = $(pwd -P)/.merge_${x}.out
Error      = $(pwd -P)/.merge_${x}.err
Queue
JOB

done


