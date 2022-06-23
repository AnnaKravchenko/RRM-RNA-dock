#!/bin/bash -i

# copy of the regular docking script 
# here i'm testing DOCKING WITH RETRAINTS 
# so here I dock just one given by user motif 

trap "kill -- -$BASHPID; $ATTRACTDIR/shm-clean" ERR EXIT
$ATTRACTDIR/shm-clean
source $MY_CONDA/etc/profile.d/conda.sh

conda activate attract

cd $1
datadir=$1
motif=$2
nb_cpu=$3
tmp=$4
proteinr=$5 
restraints=$6 

Nstart=100

if [ ! -s $RANDSEARCH//randsearch-$Nstart.dat ];then
    python $ATTRACTTOOLS/randsearch.py 2 $Nstart --fix-receptor > $RANDSEARCH//randsearch-$Nstart.dat
fi

dock(){
  motif=$1
  np=$2
  list="$LIBRARY/$motif-clust1.0r.list"
  listaa="$LIBRARY/$motif-clust1.0-aa.list"
  ligandr=`head -n 1 $list`
  ligandaa=`head -n 1 $listaa`
  Nconf=`cat $list|wc -l`

  gridparams=" --grid 1 receptorgrid.gridheader"
  parals="--np $np --chunks $np"
  
  paramsprep="$ATTRACTDIR/../attract.par $proteinr $ligandr --ens 2 $list --rest $restraints --ghost --only-rot $parals"

  # --ghost option ignores protein and all the attract.par, so it should place frag EXACTLY on the ghost, e.g. :
  #params="$ATTRACTDIR/../attract.par $proteinr $ligandr --fix-receptor --ghost --ens 2 $list --gravity 2 $gridparams --rest $restraints"

  params="$ATTRACTDIR/../attract.par $proteinr $ligandr --fix-receptor --ens 2 $list --gravity 2 $gridparams --rest $restraints"
  scoreparams="$ATTRACTDIR/../attract.par $proteinr $ligandr --score --fix-receptor --ens 2 $list --rcut 50 $parals"
 
  echo '**************************************************************'
  echo 'calculate receptorgrid grid'
  echo '**************************************************************'
  echo $ligandr
  awk '{print substr($0,58,2)}' $ligandr | sort -nu > receptorgrid.alphabet

  $ATTRACTDIR/make-grid-omp $proteinr $ATTRACTDIR/../attract.par 5.0 7.0 receptorgrid.gridheader  --shm --alphabet receptorgrid.alphabet

  echo '**************************************************************'
  echo 'Generate starting structures...'
  echo '**************************************************************'

  if [ ! -s $RANDSEARCH/randsearch-$Nstart-ens-$Nconf.dat ];then
      # option 'all' (instead of 'random') puys all awailable conformers in each starting point. so docking takes LOOOOOONGEEEEEEER like a lot 
      python $ATTRACTTOOLS/ensemblize.py $RANDSEARCH/randsearch-$Nstart.dat $Nconf 2 all > $RANDSEARCH/randsearch-$Nstart-ens-$Nconf.dat
  fi

  start=$RANDSEARCH/randsearch-$Nstart-ens-$Nconf.dat
  echo $start
  echo '**'
  echo 'orientation'
  echo '**'

  python $ATTRACTDIR/../protocols/attract.py $start $paramsprep --vmax 50 --output $motif-rot.dat

  echo '**'
  echo 'minimisation'
  echo '**'
  python $ATTRACTDIR/../protocols/attract.py $motif-rot.dat $params --vmax 100 --output $motif.dat $parals

  echo "**************************************************************"
  echo "scoring $motif"
  echo "**************************************************************"

  python $ATTRACTDIR/../protocols/attract.py $motif.dat $scoreparams --output $motif.score $parals
  
  $ATTRACTDIR/shm-clean

  conda activate mypython3
  python3 $SCRD/select-dat-perrank2.py $motif.dat --score $motif.score --percent 50 --outpscore $tmp/$motif-pc.score > $tmp/$motif-pc.dat

  rm -f $motif.dat
  rm -f $motif.score

  conda activate attract
  python2 $ATTRACTTOOLS/fill-energies.py $tmp/$motif-pc.dat $tmp/$motif-pc.score  > $tmp/$motif-scored.dat

  rm -f $tmp/$motif-pc.score
  rm -f $tmp/$motif-pc.dat

  #sys.exit()

  echo "**************************************************************"
  echo "sort & deredundant $motif"
  echo "**************************************************************"

  python $ATTRACTTOOLS/sort.py $tmp/$motif-scored.dat > $tmp/$motif-sorted.dat

  $ATTRACTDIR/fix_receptor $tmp/$motif-sorted.dat 2 --ens 0 $Nconf | python $ATTRACTTOOLS/fill.py /dev/stdin $tmp/$motif-sorted.dat > $tmp/$motif-sorted.dat-fixre
  
  rm -f $tmp/$motif-scored.dat
  rm -f $tmp/$motif-sorted.dat

  $ATTRACTDIR/deredundant $tmp/$motif-sorted.dat-fixre 2 --ens 0 $Nconf | python $ATTRACTTOOLS/fill-deredundant.py /dev/stdin $tmp/$motif-sorted.dat-fixre > $motif-sorted-dr.dat
  
  rm -f $motif.dat $motif.score
  rm -f $tmp/$motif-sorted.dat-fixre
  

  $ATTRACTTOOLS/top $motif-sorted-dr.dat 10000000 > $motif-e7.dat
  rm -f $motif-sorted-dr.dat

  set -u +e
}

dock $motif $nb_cpu

#for j in `cat motif.list`;do  #motif.list is created by extractfrag.sh
#  dock $j $nb_cpu
#done

# to make restr.txt one need numbers of lines with atoms that have to fall in the same place 
# U 6 beads
# C 6 beads 

# G 7 beads 
# A 7 beads 

# only base - 4.5A for each 
# base + backbone - tighter for the base, more free for the tail 

# ignoring 1st atom of the base if its G or A 