#!/bin/sh
#set -e

select_method=$1 # 0 -> ROF; 1 -> RDPOF; 2 -> LK
alpharof=$2
gammarof=$3
alphardpof=$4
gammardpof=$5
vx=$6
vy=$7
epipolar=$8


rdp_method="rdpof1D.exe"
rof_method="rof1D.exe"
lk_method="lucas_kanade.exe"

I1="input_0.png"
I2="input_1.png"

# check input files
gt_disp_tif="input_2.tiff"

if [ -f $gt_disp_tif ]; then
	$bin/check_input_images.py $I1 $I2 $gt_disp_tif;
else
	$bin/check_input_images.py $I1 $I2;
fi
if [ $? -eq 5 ]; then
	
	cat algo_info.txt | cut -d '=' -f2 >> demo_failure.txt;
	exit 0;
fi

# La commande suivante permet d'arrêter le script dès l'apparition d'une erreur.
set -e;

if [ $epipolar -eq 1 ];then
	$bin/epipolar.sh
fi

out_disp_tif="disp.tiff"
gt_disp_tif="input_2.tiff"
gt_disp_png="ground_truth.png"
out_disp_png="disp_map.png"
ncores="32"

case $select_method in
	
    0)
	$rof_method $I1 $I2 $out_disp_tif 8 $alpharof $gammarof 100 0.75 0.0001 1 10 0 0 0.25 mask_nan.png
    ;;
    1)
	$rdp_method $I1 $I2 $out_disp_tif 8 3 $alphardpof $gammardpof 0.2 100 0.75 0.0001 1 10 0 0 0.25 mask_nan.png
    ;;
    2)
	$lk_method $I1 $I2 $out_disp_tif $ncores $vx $vy 6 0.4 15 -1 0
    ;;

esac

if  [ ! -e $out_disp_tif ];then
	cp stderr.txt demo_failure.txt
	exit 0
fi

##Generate Ground Truth (if it is given) and returns message for the results in the DDL
if [ -f $gt_disp_tif ];then
    identify_image=$(identify $gt_disp_tif | head -n1);
    type_of_image=$(echo $identify_image | awk -F" " '{print $bin}');
	if [ "$type_of_image" == "PFM" ]; then
		$bin/save_pfm_as_image.py $gt_disp_tif
	fi
	
	$bin/generate_disp_maps.py $out_disp_tif $out_disp_png $gt_disp_tif
	$bin/generate_disp_maps.py $gt_disp_tif $gt_disp_png
	$bin/calculate_error.py $gt_disp_tif $out_disp_tif
	if [ $epipolar -eq 1 ];then
		echo "epipolar_gt='yes'" >> algo_info.txt
	else	
		echo "no_epipolar_gt='yes'" >> algo_info.txt
	fi
else
	$bin/generate_disp_maps.py $out_disp_tif $out_disp_png
    if [ $epipolar -eq 1 ];then
		echo "epipolar_no_gt='yes'" >> algo_info.txt
	else	
		echo "no_epipolar_no_gt='yes'" >> algo_info.txt
	fi
fi

