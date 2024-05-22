#!/bin/bash

function main() {
  echo $1
  ssh root@$1 'apt-get update && apt-get -y upgrade'
  ssh root@$1 'yum -y update'
  exit 0
}

main