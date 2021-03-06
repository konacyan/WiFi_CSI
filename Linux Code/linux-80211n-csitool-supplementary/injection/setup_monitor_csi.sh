#!/usr/bin/sudo /bin/bash

#sudo service network-manager stop
sudo modprobe -r iwldvm iwlwifi mac80211
while [ $? -ne 0 ]
do
    sudo modprobe -r iwldvm iwlwifi mac80211
done

modprobe -r iwlwifi mac80211 cfg80211
while [ $? -ne 0 ]
do
    modprobe -r iwlwifi mac80211 cfg80211
done

sudo modprobe iwlwifi connector_log=0x1
while [ $? -ne 0 ]
do
    sudo modprobe iwlwifi connector_log=0x1
done


:<<EOF
if [ "$#" -ne 2 ]; then
    echo "Going to use default settings!"
    chn=64
    bw=HT20
else
    chn=$1
    bw=$2
fi
EOF
chn=60
bw=HT40+
sudo service network-manager restart
sudo service network-manager stop
sudo iwconfig wlan0 mode monitor
#iwconfig wlan0 mode monitor 2>/dev/null 1>/dev/null

for a in {1..5} 
do
	if [ $? -ne 0 ];then
    	sudo service network-manager restart
    	sudo service network-manager stop
    	sudo iwconfig wlan0 mode monitor
		sleep 1	
		iwconfig	
	fi
done

#while [ $? -ne 0 ]
#do
#    sudo service network-manager restart
#    sudo service network-manager stop
#    sudo iwconfig wlan0 mode monitor 

#done

ifconfig wlan0 up 2>/dev/null 1>/dev/null
while [ $? -ne 0 ]
do
  ifconfig wlan0 up 2>/dev/null 1>/dev/null
done

iw wlan0 set channel $chn $bw
echo Channel Num:$chn, Channel Bandwidth:$bw

