# Human reads removal

This image is defined to remove human reads with a 
Human genome (masked).

 - [Download the reference](https://drive.google.com/file/d/0B3llHR93L14wd0pSSnFULUlhcUk/edit), decompress it and place it in `nohost/`
 - Build the index with **bbmap** (for example running build_reference.sh)
 - Build the image: `sudo singularity build nohost.simg nohost.def`
