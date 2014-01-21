#!/bin/bash -l
#SBATCH -J 3PTI_warm
#SBATCH -p devel
#SBATCH -n 4
#SBATCH -t 00:15:00 
#SBATCH -A g2013239

module load intel openmpi

mpirun ~/glob/private/Q/Q5T/Qdyn5p 3pti_min.inp > 3pti_min.log
mpirun ~/glob/private/Q/Q5T/Qdyn5p 3pti_warm1K.inp>3pti_warm1K.log
mpirun ~/glob/private/Q/Q5T/Qdyn5p 3pti_warm20K.inp>3pti_warm20K.log
mpirun ~/glob/private/Q/Q5T/Qdyn5p 3pti_warm100K.inp>3pti_warm100K.log
mpirun ~/glob/private/Q/Q5T/Qdyn5p 3pti_warm200K.inp>3pti_warm200K.log
mpirun ~/glob/private/Q/Q5T/Qdyn5p 3pti_warm300K.inp>3pti_warm300K.log

mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_100.inp >fep_100.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_098.inp >fep_098.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_096.inp >fep_096.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_094.inp >fep_094.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_092.inp >fep_092.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_090.inp >fep_090.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_088.inp >fep_088.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_086.inp >fep_086.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_084.inp >fep_084.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_082.inp >fep_082.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_080.inp >fep_080.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_078.inp >fep_078.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_076.inp >fep_076.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_074.inp >fep_074.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_072.inp >fep_072.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_070.inp >fep_070.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_068.inp >fep_068.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_066.inp >fep_066.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_064.inp >fep_064.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_062.inp >fep_062.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_060.inp >fep_060.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_058.inp >fep_058.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_056.inp >fep_056.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_054.inp >fep_054.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_052.inp >fep_052.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_050.inp >fep_050.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_048.inp >fep_048.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_046.inp >fep_046.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_044.inp >fep_044.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_042.inp >fep_042.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_040.inp >fep_040.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_038.inp >fep_038.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_036.inp >fep_036.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_034.inp >fep_034.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_032.inp >fep_032.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_030.inp >fep_030.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_028.inp >fep_028.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_026.inp >fep_026.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_024.inp >fep_024.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_022.inp >fep_022.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_020.inp >fep_020.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_018.inp >fep_018.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_016.inp >fep_016.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_014.inp >fep_014.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_012.inp >fep_012.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_010.inp >fep_010.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_008.inp >fep_008.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_006.inp >fep_006.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_004.inp >fep_004.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_002.inp >fep_002.log
mpirun /home/alexb/glob/private/Q/Q5T/Qdyn5p fep_000.inp >fep_000.log
