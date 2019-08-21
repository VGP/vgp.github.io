#!/bin/sh

if [ -z $NCBI_API_KEY ] ; then
  echo NCBI_API_KEY not set.
  exit
fi

module load edirect

if [ ! -e genbank.xml ] ; then
  echo Fetching genbank.xml.

  esearch -db bioproject -q 'PRJNA489243' \
  | \
  elink -db bioproject -target assembly -name bioproject_assembly_all \
  | \
  esummary \
  > genbank.xml
fi

if [ ! -e genbank.xml.map ] ; then
  echo Parsing genbank.xml.

  #  Extract elements 'Genbank', 'AssemblyName' and 'AssemblyType' from all
  #  'DocumentSummary' elements IF it contains a `Synonym` element (which is
  #  where the 'Genbank' element is).

  xtract \
    -input genbank.xml \
    -pattern DocumentSummary -if Synonym -element Genbank AssemblyName AssemblyType \
  | \
  sort -k2,2 \
  > genbank.map
fi

exit 0
