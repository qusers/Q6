# Q6: A comprehensive simulation package for molecular dynamics simulations and 
# free energy calculations, including empirical valence bond simulations, 
# linear interaction energy calculations, and free energy perturbation.
# 
# Copyright © 2017 Johan Åqvist, John Marelius, Shina Caroline Lynn Kamerlin and Paul Bauer
# 
# This program is free software; you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free 
# Software Foundation; either version 2 of the License, or any later version.
# 
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to the Free Software Foundation, Inc., 51 Franklin 
# Street, Fifth Floor, Boston, MA  02110-1301, USA. Also add information on 
# how to contact you by electronic and paper mail.
################################################################################
#  Q Build Auto Incrementing Buildnumber                                       #
################################################################################

MAYOR_VER=6
MINOR_VER=0
OBJECTS    = $(shell ls *.o)
BUILD_NUMBER_LDFLAGS  = -DBUILD_USERNAME=\"$(USER)\"
BUILD_NUMBER_LDFLAGS += -DBUILD_HOSTNAME=\"$(shell hostname)\"
BUILD_NUMBER_LDFLAGS += -DBUILD_DATE=\"$$(date +'%Y%m%d')\"
BUILD_NUMBER_LDFLAGS += -DBUILD_SOURCE=\"$(shell git rev-parse HEAD)\"
BUILD_NUMBER_LDFLAGS += -DBUILD_NUMBER=\"$(MAYOR_VER).$(MINOR_VER).$$(cat $(BUILD_NUMBER_FILE))\"
VERSIO2='${shell $(FC) --version|head -n 1| sed "s/\n//g" | sed "s/ /\ /g" |sed "s/)//g" |sed "s/(//g"}'
BUILD_NUMBER_LDFLAGS += -DBUILD_COMPILER=\"$(VERSIO2)\"


$(BUILD_NUMBER_FILE): $(OBJECTS)
	@if ! test -f $(BUILD_NUMBER_FILE); then echo 0 > $(BUILD_NUMBER_FILE); fi
	@echo $$(($$(cat $(BUILD_NUMBER_FILE)) + 1)) > $(BUILD_NUMBER_FILE)


