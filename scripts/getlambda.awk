#!/usr/bin/awk -f

{
        if ( $2=="Part" && $3=="3:" ) {
                startnow=2
        }


	if ( startnow==1 ) {
		print " ",$1,"      ",$3,"  ",$7
	}

	if ( $2=="Lambda(1)" && $4=="Energy" && $5=="gap" ) {
		startnow=1
		print "#        Lambda	energygap	pts"
	}

}
