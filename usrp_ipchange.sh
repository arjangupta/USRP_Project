#!/bin/bash
sudo stop network-manager
sudo ifconfig eth0 192.168.40.1 netmask 255.255.255.0
