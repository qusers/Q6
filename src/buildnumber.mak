################################################################################
#  Q Build Auto Incrementing Buildnumber                                       #
################################################################################

MAYOR_VER=5
MINOR_VER=10
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


