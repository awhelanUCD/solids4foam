#!/bin/bash

FOVALUE=$1
SNVALUE=$2
#YOUNGSMOD=$3
FNVALUE=$3
FCVALUE=$4
EPSILONN=$4

#sed -i "s/YOUNGSMOD/${YOUNGSMOD}/" mechanicalProperties
#sed -i "s/PRATIO/${PRATIO}/" mechanicalProperties

sed -i "s/FOVALUE/${FOVALUE}/" mechanicalProperties
sed -i "s/SNVALUE/${SNVALUE}/" mechanicalProperties
sed -i "s/FNVALUE/${FNVALUE}/" mechanicalProperties
sed -i "s/FCVALUE/${FCVALUE}/" mechanicalProperties
sed -i "s/EPSILONN/${EPSILONN}/" mechanicalProperties





