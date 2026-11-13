# Linux build (author machine paths): prepares the source tree, then builds everything in an Ubuntu 22.04 chroot.  Releases are
# built on Ubuntu 20.04 (glibc 2.31) instead; a 22.04 chroot's GPU programs need glibc 2.34, so run this with MAX_GLIBC=2.34.
J=/root/jammy
mountpoint -q $J/proc || mount -t proc proc $J/proc
mountpoint -q $J/dev || mount --bind /dev $J/dev
OUT_FATBIN=/root/nf/nf_cuda_fatbin.h bash /mnt/c/Users/Alp/Documents/passiflora/build/prepare-src.sh
mkdir -p $J/pk/src && rm -rf $J/pk/src/* $J/root/pkb && cp /root/nf/build_src/GetDecicsSrc/* $J/pk/src/
cp "$(dirname "$0")/build-linux-chroot-inner.sh" $J/pk/inner.sh
chroot $J /bin/bash /pk/inner.sh > /root/nf/pk_inner.log 2>&1; tail -12 /root/nf/pk_inner.log
mkdir -p /root/nf/jammy_out && cp $J/root/pkb/GetDecics_cpu $J/root/pkb/GetDecics_gpu $J/root/pkb/GetDecics_cuda /root/nf/jammy_out/ && ls -la /root/nf/jammy_out
