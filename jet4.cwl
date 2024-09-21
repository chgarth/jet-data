#! /usr/bin/env cwl-runner 
cwlVersion: v1.2
class: CommandLineTool

baseCommand: [ 'sh', 'script.sh' ]

requirements:
  DockerRequirement:
    dockerImageId: gerris-custom
  ResourceRequirement:
    coresMin: 16
    coresMax: 16
  InitialWorkDirRequirement:
    listing:
      - $(inputs.postproc)
      - entryname: script.sh
        entry: |
          gerris3D -m -s 1 $1 > tmp.gfs
          mpirun -n 16 --oversubscribe --allow-run-as-root gerris3D tmp.gfs


inputs:
  gfsfile:
    type: File
    inputBinding:
      position: 1
  postproc:
    type: File
    inputBinding:
      position: 2
    default: 
      class: File
      path: gerris-merge-nrrd.py

outputs: 
  cdb:
    type: File
    outputBinding: 
      glob: data.csv
  data:
    type: File[]
    outputBinding:
      glob: [ 'v_????.vti' ] 