#!/usr/bin/awk -f

{
	if ( startnow==1 ) {
                print "	",$2,"      ",$4,"  ",$8,"	",$3
        }

	if ( $2=="bin" && $3=="energy" && $4=="gap" ) {
		startnow++
#		print startnow
		print "#	energygap	dG(norm)	r_xy	dG(nonnorm)"
	}
#	if ( startnow==1 && $3=="energy" ) {
#		 startnow++
#		print startnow
#	}
}
