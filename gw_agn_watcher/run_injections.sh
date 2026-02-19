#!/usr/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=16GB
#SBATCH --time=48:00:00

module load miniforge3
source activate cosmo
export LAL_DATA_PATH=/home/imcmah/data/lalsuite-waveform-data

lalapps_inspinj \
	--output inj.xml \
	--f-lower 25 \
	--waveform TaylorF2threePointFivePN \
	--t-distr uniform --time-step 1 --gps-start-time 1000000000 --gps-end-time 1000100000 \
	--m-distr log --min-mass1 1.2 --max-mass1 75 --min-mass2 1.2 --max-mass2 75 \
	--d-distr volume --min-distance 1 --max-distance 6000e3 \
	--l-distr random \
	--i-distr uniform \
	--enable-spin --min-spin1 0.0 --max-spin1 1.0 --min-spin2 0.0 --max-spin2 1.0 --aligned \
	--disable-milkyway

bayestar-sample-model-psd \
	-o psd.xml \
	--H1=aLIGO175MpcT1800545 --L1=aLIGO175MpcT1800545 --V1=aLIGOAdVO4T1800545

bayestar-realize-coincs \
	-o coinc.xml \
	inj.xml --reference-psd psd.xml \
	--detector H1 L1 V1 \
	--measurement-error gaussian-noise \
	--snr-threshold 2.0 \
	--net-snr-threshold 6.0 \
	--min-triggers 2 \

bayestar-localize-coincs coinc.xml

ligolw_sqlite --preserve-ids --replace --database coinc.sqlite coinc.xml

ligo-skymap-stats \
	-o bayestar.tsv \
	--database coinc.sqlite \
	*.fits \
	--contour 50 90 \
	--area 10 100 \
	-j 8

for file in ./*.fits; do
    base="${file%.fits}"
    output="${base}_downsample.fits.gz"
    ligo-skymap-flatten --nside 8 "$file" "$output"
done

zip downsamples.zip ./*_downsample.fits.gz
zip multiorder_skymaps.zip ./*.fits
rm ./*_downsample.fits.gz
rm ./*.fits
