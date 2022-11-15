#!/bin/bash
# 'RunTrackCutStudyG4.sh'
# Derek Anderson
# 11.11.2022
#
# Short to script to run the macro
# 'Fun4All_G4_sPHENIX_ForTrackCutStudy.C'
# over a (small) set of files.

# declare i/o lists
declare -a inFiles
declare -a embFiles
declare -a outFiles

# input files
inFiles[0]="\"https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/G4Hits_sHijing_0_20fm-0000000040-00020.root\""
inFiles[1]="\"https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/G4Hits_sHijing_0_20fm-0000000040-00021.root\""
inFiles[2]="\"https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/G4Hits_sHijing_0_20fm-0000000040-00022.root\""
inFiles[3]="\"https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/G4Hits_sHijing_0_20fm-0000000040-00023.root\""
inFiles[4]="\"https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/G4Hits_sHijing_0_20fm-0000000040-00024.root\""
inFiles[5]="\"https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/G4Hits_sHijing_0_20fm-0000000040-00025.root\""
inFiles[6]="\"https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/G4Hits_sHijing_0_20fm-0000000040-00026.root\""
inFiles[7]="\"https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/G4Hits_sHijing_0_20fm-0000000040-00027.root\""
inFiles[8]="\"https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/G4Hits_sHijing_0_20fm-0000000040-00028.root\""
inFiles[9]="\"https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/G4Hits_sHijing_0_20fm-0000000040-00029.root\""

# input embed files
embFiles[0]="\"https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/G4Hits_sHijing_0_20fm-0000000040-00020.root\""
embFiles[1]="\"https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/G4Hits_sHijing_0_20fm-0000000040-00021.root\""
embFiles[2]="\"https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/G4Hits_sHijing_0_20fm-0000000040-00022.root\""
embFiles[3]="\"https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/G4Hits_sHijing_0_20fm-0000000040-00023.root\""
embFiles[4]="\"https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/G4Hits_sHijing_0_20fm-0000000040-00024.root\""
embFiles[5]="\"https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/G4Hits_sHijing_0_20fm-0000000040-00025.root\""
embFiles[6]="\"https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/G4Hits_sHijing_0_20fm-0000000040-00026.root\""
embFiles[7]="\"https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/G4Hits_sHijing_0_20fm-0000000040-00027.root\""
embFiles[8]="\"https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/G4Hits_sHijing_0_20fm-0000000040-00028.root\""
embFiles[9]="\"https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/G4Hits_sHijing_0_20fm-0000000040-00029.root\""

# output files
outFiles[0]="\"sPhenixG4_test00020.root\""
outFiles[1]="\"sPhenixG4_test00021.root\""
outFiles[2]="\"sPhenixG4_test00022.root\""
outFiles[3]="\"sPhenixG4_test00023.root\""
outFiles[4]="\"sPhenixG4_test00024.root\""
outFiles[5]="\"sPhenixG4_test00025.root\""
outFiles[6]="\"sPhenixG4_test00026.root\""
outFiles[7]="\"sPhenixG4_test00027.root\""
outFiles[8]="\"sPhenixG4_test00028.root\""
outFiles[9]="\"sPhenixG4_test00029.root\""

# other parameters
nEvt=1
nSkip=0

# loop over files
(( nFile=0 ))
for input in ${inFiles[@]}; do
  root -b -q "Fun4All_G4_sPHENIX_ForTrackCutStudy.C($nEvt, $input, ${outFiles[$nFile]}, ${embFiles[$nFile]}, $nSkip)"
  (( nFile++ ))
done

# delete arrays
unset inFiles
unset embFiles
unset outFiles

# end -------------------------------------------------------------------------
