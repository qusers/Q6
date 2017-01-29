#!/bin/bash

python convert_amber2q.py
diff qamber14_gen.lib ../Q/ff-qamber14/qamber14.lib > tmp.lib.diff

diff qamber14_gen.prm ../Q/ff-qamber14/qamber14.prm > tmp.prm.diff

echo
echo "## Comparing generated library vs pre-made"
diff tmp.lib.diff qamber14.lib.diff
diff tmp.prm.diff qamber14.prm.diff
echo "## End library diffs"

