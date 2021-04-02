#!/usr/bin/sudo /bin/bash
:<<EOF
modprobe -r iwlwifi mac80211 cfg80211
modprobe iwlwifi debug=0x40000
ifconfig wlan0 2>/dev/null 1>/dev/null
while [ $? -ne 0 ]
do
	        ifconfig wlan0 2>/dev/null 1>/dev/null
done
iw dev wlan0 interface add mon0 type monitor
iw mon0 set channel $1 $2
ifconfig mon0 up
EOF


:<<EOF
#!/usr/bin/sudo /bin/bash
service network-manager stop
WLAN_INTERFACE=$1
SLEEP_TIME=2
modprobe iwlwifi debug=0x40000
if [ "$#" -ne 3 ]; then
    echo "Going to use default settings!"
    chn=64
    bw=HT20
else
    chn=$2
    bw=$3
fi
sleep $SLEEP_TIME
ifconfig $WLAN_INTERFACE 2>/dev/null 1>/dev/null
while [ $? -ne 0 ]
do
    ifconfig $WLAN_INTERFACE 2>/dev/null 1>/dev/null
done
sleep $SLEEP_TIME
echo "Add monitor mon0....."
iw dev $WLAN_INTERFACE interface add mon0 type monitor
sleep $SLEEP_TIME
echo "Bringing $WLAN_INTERFACE down....."
ifconfig $WLAN_INTERFACE down
while [ $? -ne 0 ]
do
    ifconfig $WLAN_INTERFACE down
done
sleep $SLEEP_TIME
echo "Bringing mon0 up....."
ifconfig mon0 up
while [ $? -ne 0 ]
do
    ifconfig mon0 up
done
sleep $SLEEP_TIME
echo "Set channel $chn $bw....."
iw mon0 set channel $chn $bw
EOF


:<<EOF
# service network-manager stop
sudo modprobe -r iwldvm iwlwifi mac80211
modprobe -r iwlwifi mac80211 cfg80211
modprobe iwlwifi debug=0x40000
if [ "$#" -ne 2 ]; then
    echo "Going to use default settings!"
    chn=64
    bw=HT20
else
    chn=$1
    bw=$2
fi
ifconfig wlan0 2>/dev/null 1>/dev/null
while [ $? -ne 0 ]
do
            ifconfig wlan0 2>/dev/null 1>/dev/null
done
iw dev wlan0 interface add mon0 type monitor

ifconfig wlan0 down
while [ $? -ne 0 ]
do
    ifconfig wlan0 down
done
ifconfig mon0 up
while [ $? -ne 0 ]
do
           ifconfig mon0 up
done

iw mon0 set channel $chn $bw
EOF




sudo service network-manager stop
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

modprobe iwlwifi debug=0x40000
while [ $? -ne 0 ]
do
           modprobe iwlwifi debug=0x40000
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


ifconfig wlan0 2>/dev/null 1>/dev/null
while [ $? -ne 0 ]
do
            ifconfig wlan0 2>/dev/null 1>/dev/null
done
iw dev wlan0 interface add mon0 type monitor

ifconfig wlan0 down
while [ $? -ne 0 ]
do
    ifconfig wlan0 down
done
ifconfig mon0 up
while [ $? -ne 0 ]
do
           ifconfig mon0 up
done

iw mon0 set channel $chn $bw
echo Channel Number: $chn, Channel Bandwidth: $bw

sudo echo 0x1c911 | sudo tee /sys/kernel/debug/ieee80211/phy0/iwlwifi/iwldvm/debug/monitor_tx_rate

