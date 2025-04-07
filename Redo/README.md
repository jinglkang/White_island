# Rename the sample ID for data uploading 
```bash
# kangjingliang@KangdeMacBook-Pro-2 一  4 07 11:05:18 ~/Desktop
less 1.txt | perl -alne '$hash{"Cs"}="Control";$hash{"Vn"}="CO2 vent1";$hash{"Vs"}="CO2 vent2";@a=split /\_/;$hash1{"Cs"}="C";$hash1{"Vn"}="V1";$hash1{"Vs"}="V2";$id=$a[0]."_".$hash1{$a[1]}."_".$a[2];print "$a[0]\t$hash{$a[1]}\t$id"' >2.txt
```
