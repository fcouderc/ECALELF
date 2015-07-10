#!/bin/bash
queue=$1
nEventsPerToy=factorizedSherpaFixed_DataSeedFixed_smooth_autobin_smearingEt_${queue}
#nEventsPerToy=factorized
nSmearToy=30 #it quantifies the over-sampling
commonCut=Et_25-trigger-noPF

if [ -z "${queue}" ];then local=y; fi

for regionsFile in scaleStep2smearing_1
  do
  dir=test/dato/fitres/toys/${regionsFile}

for scale in 1.00 1.01 0.99 1.02 0.98 1.05 0.95
do
for constAlpha in  0.02-0.00 0.04-0.00 0.01-0.00 
  do
  const=`echo $constAlpha | cut -d '-' -f 1`
  alpha=`echo $constAlpha | cut -d '-' -f 2`

  echo "[[[[[[[[[[[ const:alpha ]]]]]]]]]]] ${const}:${alpha}"
  baseDir=${dir}/${nEventsPerToy}/${const}-${alpha}-${scale}/
  mkdir -p $baseDir 

  cat data/validation/toyMC-smearingMethod.dat > ${baseDir}/invertedToyMC.dat

#mcToy.txt is the initFile passed to ZFitter
  cat > $baseDir/${alphaConst}/mcToy.txt <<EOF
scale_absEta_1_1.4442-gold-${commonCut} 	=  ${scale} +/- 0.0050000 L(0.9 - 1.1) // [GeV]
EOF


  cat data/validation/toyMC-smearingMethod.dat > ${baseDir}/toyMC.dat

  for nToys in `seq 1 100`; 
  do 
      newDir=${baseDir}/${alphaConst}/${nSmearToy}/${nToys}/
      mkdir -p $newDir 
      ls $newDir
      
      if [ -z "${local}" ];then
	  bsub -q ${queue} -oo ${newDir}/stdout.log -eo ${newDir}/stderr.log -J "$regionsFile-$const-$alpha" \
	      "cd /afs/cern.ch/user/g/gfasanel/new_version_ECALELF/CMSSW_7_5_0_pre4/src/Calibration/ZFitter; eval \`scramv1 runtime -sh\`; uname -a;  echo \$CMSSW_VERSION; ./bin/ZFitter.exe -f ${baseDir}/toyMC.dat \
      --regionsFile=data/regions/${regionsFile}.dat \
	--commonCut=${commonCut} --smearerFit --outDirFitResData=$newDir \
        --constTermFix --alphaGoldFix --smearerType=profile --noPU --smearingEt \
	--initFile=${baseDir}/${alphaConst}/mcToy.txt --profileOnly --plotOnly --runToy --nSmearToy=${nSmearToy} \
        --autoBin        > ${newDir}/log2.log"
      else
	  echo "#============================================================ Toy = $nToys"
	  ./bin/ZFitter.exe -f ${baseDir}/toyMC.dat \
	      --regionsFile=data/regions/${regionsFile}.dat \
	      --commonCut=${commonCut} --smearerFit --outDirFitResData=$newDir \
	      --constTermFix --smearerType=profile --noPU  --smearingEt \
	      --initFile=${baseDir}/${alphaConst}/mcToy.txt \
	      --plotOnly --profileOnly --runToy --eventsPerToy=0 --nSmearToy=${nSmearToy} --autoBin |tee ${newDir}/nSmearToy_${nSmearToy}.log
	  exit 0
      fi
    done
    wait 
  done
  wait
done

done

exit 0

#Qui sotto va capito cosa vuoi fare tu

for file in  test/dato/fitres/Hgg_Et-toys/scaleStep2smearing_9/factorizedSherpaFixed_DataSeedFixed_smooth_cmscaf1nd/0.01-0.00/15/*/log2.log; do grep DUMP $file > `dirname $file`/dumpNLL.dat; done

for index in `seq 2 50`; do sed -i "s|\[DUMP NLL\]|$index\t |" test/dato/fitres/Hgg_Et-toys/scaleStep2smearing_9/factorizedSherpaFixed_DataSeedFixed_smooth_cmscaf1nd/0.01-0.00/15/$index/dumpNLL.dat; done

cat test/dato/fitres/Hgg_Et-toys/scaleStep2smearing_9/factorizedSherpaFixed_DataSeedFixed_smooth_cmscaf1nd/0.01-0.00/15/*/dumpNLL.dat > dumpNLL.dat

cat  dumpNLL.dat  | awk -F '\t' '{cat[$2]+=$3;cat2[$2]+=$3*$3; n[$2]++;};END{for(i in cat){print i, n[i], cat[i]/n[i], cat2[i]/n[i]-cat[i]/n[i]*cat[i]/n[i]}}' | awk '($5 > 5){print $0}' 
