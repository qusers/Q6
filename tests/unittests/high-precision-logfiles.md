#create high precisoin logfiles, enabled with a compilation flag.
to shorten tests, and to see temperature drifts early, the precison of output, in 
temperature, engries etc. should be increased to single, doulbe or quad precision, depending on the precision 
the code was  compiled (new feature that is curretly implemented by irek).
to not break compatibilyty with logfile reading  utilities, this should be made an optional flag.

(it should be a switch/keyword at runtime, so that only when we run a test or when the user wants huge logfiles, the precision is increased.
