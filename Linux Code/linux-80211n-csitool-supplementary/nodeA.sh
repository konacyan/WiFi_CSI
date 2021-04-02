#!/usr/bin/sudo /bin/bash
echo "=================================================="
echo -e "\033[47;30m Run nodeA \033[0m"
echo -e "\033[31m Node A Sending Data\033[0m"
sudo injection/setup_inject.sh
sudo injection/random_packets


echo -e "\033[31m Node A Receiving Data\033[0m"
sudo injection/setup_monitor_csi.sh
data_name=$1
#touch Data/$1
sudo netlink/log_to_file Data/B2A$1.dat

echo -e "\033[31m Node A Sending Data\033[0m"
sudo injection/setup_inject.sh
sudo injection/random_packets

echo -e "\033[47;30m Run nodeA end \033[0m"
