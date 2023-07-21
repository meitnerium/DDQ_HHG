#! /bin/bash
echo "set terminal gif animate delay 10" > animchi.gp
echo "set output \"animate.gif\"" >> animchi.gp
i=2001
max=42000
while [ $i -le $max ] ; do
  echo "plot \"./fort.$i\" u 1:2 w l, \"./fort.$i\" u 1:3 w l" >> animchi.gp
  let i=$i+100
done

#i=85001
#max=125000

#while [ $i -lt $max ] ; do
#  echo "plot \"./fort.$i\" u 1:2 w l" > anim.gp

#done

