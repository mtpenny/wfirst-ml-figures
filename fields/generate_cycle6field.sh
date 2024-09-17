#!/bin/bash

codedir=${0%/*}
#Argument processing
newrot=9999
roll=0
centers=0
bore=0
for i in "$@"; do
    ms="^-rot=([+-]{0,1}[0-9\.]+)$"
    if [[ $i =~ $ms ]]; then
	echo
	newrot=${BASH_REMATCH[1]}
    fi
    ms="^-roll=([+-]{0,1}[0-9\.]+)$"
    if [[ $i =~ $ms ]]; then
	echo
	roll=${BASH_REMATCH[1]}
    fi
    if [[ $i =~ ^-{1,2}h$ ]] || [[ $i =~ ^-{1,2}help$ ]]; then
	echo "Usage: $0 {-rot=<angle(deg)>} {-roll=<angle(deg)>} {-bore //coordinates are for boresight center} {-h //help} {-c //only centers}"
	exit
    fi
    if [[ $i =~ ^-{1,2}c$ ]]; then
	centers=1
    fi
    if [[ $i =~ ^-{1,2}bore$ ]]; then
	bore=1
    fi
done
#echo $newrot > /dev/stderr

awk -v newrot=$newrot -v cent=$centers -v roll=$roll -v bore=$bore -v CONVFMT="%.9g" -v OFMT="%.9g" 'BEGIN{
pi=2*atan2(1,0)
BORE_ROT=(300+roll)*pi/180.0
BORE_OFF=0.479
CHIPPIX_X=4088
CHIPPIX_Y=4088
PIXSCALE_X=0.11
PIXSCALE_Y=0.11
ox[0]=-0.5; ox[1]=-0.5; ox[2]=0.5; ox[3]=0.5; ox[4]=-0.5; ox[5]=0;
oy[0]=-0.5; oy[1]=0.5; oy[2]=0.5; oy[3]=-0.5; oy[4]=-0.5; oy[5]=0;
}{
  if($1~"BORE_ROT"){BORE_ROT=($2+roll)*pi/180.0; if(newrot!=9999) BORE_ROT=newrot*pi/180.0}
  if($1~"BORE_OFF") BORE_OFF=$2
  if($1~"CHIPPIX"){CHIPPIX_X=$2; CHIPPIX_Y=$(NF)}
  if($1~"PIXSCALE"){PIXSCALE_X=$2; PIXSCALE_Y=$(NF)}

  if(NF==3&&!($0~"#")){
    
    for(i=(cent==0?0:5);i<(cent==0?5:6);i++){
       x = $2 + ox[i]*CHIPPIX_X*PIXSCALE_X/3600.0; y = $3 + oy[i]*CHIPPIX_Y*PIXSCALE_Y/3600.0;
       nx = (bore==1?BORE_OFF*cos(BORE_ROT):0) + x*cos(0.5*pi+BORE_ROT) - y*sin(0.5*pi+BORE_ROT)
       ny = (bore==1?BORE_OFF*sin(BORE_ROT):0) + x*sin(0.5*pi+BORE_ROT) + y*cos(0.5*pi+BORE_ROT)
       print $1,nx,ny
    }
    print ""
  }
}' $codedir/sca_centers.txt
