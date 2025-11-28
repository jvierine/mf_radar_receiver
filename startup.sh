sudo sysctl -w net.core.wmem_max=2500000
sudo sysctl -w net.core.rmem_max=100000000


nohup thor.py -m 192.168.20.2 -d "A:A A:B" -c ch3,ch4 -f 2.78e6 -r 1e6 /data1/mfraw/ &> log.20.2  &
nohup thor.py -m 192.168.10.2 -d "A:A A:B" -c ch1,ch2 -f 2.78e6 -r 1e6 /data1/mfraw/ &> log.10.2 &
