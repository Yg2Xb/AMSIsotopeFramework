#!/bin/bash
#nohup root -b -l -q ../rootlogon.C 'buildCDFLookupTables.cpp("Beryllium")' > be.log 2>&1 &
#nohup root -b -l -q ../rootlogon.C 'buildCDFLookupTables.cpp("Boron")' > b.log 2>&1 &
#nohup root -b -l -q ../rootlogon.C 'buildCDFLookupTables.cpp("Carbon")' > c.log 2>&1 &
#nohup root -b -l -q ../rootlogon.C 'buildCDFLookupTables.cpp("Nitrogen")' > n.log 2>&1 &
#nohup root -b -l -q ../rootlogon.C 'buildCDFLookupTables.cpp("Oxygen")' > o.log 2>&1 &
root -b -l -q ../rootlogon.C 'buildCDFLookupTables.cpp("Lithium")'