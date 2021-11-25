# gr-lora
This is the first open-source implementation of real-time [LoRa](https://en.wikipedia.org/wiki/LoRa) PHY collision decoding.
The collision decoding algorithm **Pyramid** used here comes from the following paper:

> Xu, Zhenqiang and Xie, Pengjin and Wang, Jiliang. Pyramid: Real-Time LoRa Collision Decoding with Peak Tracking. In Proceedings of IEEE INFOCOM. 2021.

This work is developed based on Matt Knight's [gr-lora](https://github.com/jkadbear/gr-lora) with many improvements on reverse engineering.
A MATLAB script for prototype verification is also provided ([LoRaPHY](https://github.com/jkadbear/LoRaPHY)).

## Installation
### Docker Installation
The easiest way to install `gr-lora` is using docker. Simply run
```
git clone https://github.com/jkadbear/gr-lora.git .
./docker_run.sh # or sudo ./docker_run.sh
```

The script will download the docker container and show you a shell with `gr-lora` installed:
```
root@6f8d0966cba4:/src/examples#
```

## Usage
Run `gnuradio-companion` in the above shell.
Open the grc file `/src/examples/rx_file_collision.grc` and run it.
You will see the following output:
```
* MESSAGE DEBUG PRINT PDU VERBOSE *
()
pdu_length = 12
contents = 
0000: 06 30 f0 01 02 03 04 05 06 05 08 01 
***********************************
* MESSAGE DEBUG PRINT PDU VERBOSE *
()
pdu_length = 11
contents = 
0000: 05 30 00 07 07 07 07 07 e7 6b 01 
***********************************
```
From the collision signal we successfully decode two packets with data `01 02 03 04 05 06` and `07 07 07 07 07`!
(Other bytes are header, CRC checksum and integrity indicator.)

## TODO
- Decoding multiple channels simultaneously
- Implement upper layers (LoRaWAN)

## References

1. Zhenqiang Xu, Pengjin Xie, Jiliang Wang. Pyramid: Real-Time LoRa Collision Decoding with Peak Tracking. In Proceedings of IEEE INFOCOM. 2021: 1-9.
2. Zhenqiang Xu, Shuai Tong, Pengjin Xie, Jiliang Wang. FlipLoRa: Resolving Collisions with Up-Down Quasi-Orthogonality. In Proceedings of IEEE SECON. 2020: 1-9.
3. Shuai Tong, Zhenqiang Xu, Jiliang Wang. CoLoRa: Enabling Multi-Packet Reception in LoRa. In Proceedings of IEEE INFOCOM. 2020: 2303-2311.
4. Shuai Tong, Jiliang Wang, Yunhao Liu. Combating packet collisions using non-stationary signal scaling in LPWANs. In Proceedings of Proceedings of ACM MobiSys. 2020: 234-246.
5. Yinghui Li, Jing Yang, Jiliang Wang. DyLoRa: Towards Energy Efficient Dynamic LoRa Transmission Control. In Proceedings of IEEE INFOCOM. 2020: 2312-2320.
6. Qian Chen, Jiliang Wang. AlignTrack: Push the Limit of LoRa Collision Decoding. In Proceedings of IEEE ICNP. 2021.
7. Jinyan Jiang, Zhenqiang Xu, Jiliang Wang. Long-Range Ambient LoRa Backscatter with Parallel Decoding. In Proceedings of ACM MobiCom. 2021.
8. Shuai Tong, Zilin Shen, Yunhao Liu, Jiliang Wang. Combating Link Dynamics for Reliable LoRa Connection in Urban Settings. In Proceedings of ACM MobiCom. 2021.
9. Chenning Li, Hanqing Guo, Shuai Tong, Xiao Zeng, Zhichao Cao, Mi Zhang, Qiben Yan, Li Xiao, Jiliang Wang, Yunhao Liu. NELoRa: Towards Ultra-low SNR LoRa Communication with Neural-enhanced Demodulation. In Proceedings of ACM SenSys. 2021.
