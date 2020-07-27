#!/bin/bash
rm *_small.png || true

for i in *.png; 
do
   convert -resize  25% $i $(echo $i| sed 's/.png//')_small.png;
done
