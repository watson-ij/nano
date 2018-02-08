# nanoAOD setup
quick setup for analysis
```
scram p -n nanoAOD CMSSW CMSSW_9_4_2
cd nanoAOD/src
cmsenv
git clone git@github.com:watson-ij/nano.git 
scram b -j 20
```


# nanoAOD customiser for prod
run with  --customise nano/nanoAOD/nano_cff.customise
 - added cmeson producer
 - added additional muon variables

```
scram p -n nanoAOD CMSSW CMSSW_9_4_2
cd nanoAOD/src
cmsenv
git-cms-init -q
git cms-merge-topic cms-nanoAOD:master
git checkout -b nanoAOD cms-nanoAOD/master
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools

git clone git@github.com:watson-ij/nano.git 

scram b -j 20
```
