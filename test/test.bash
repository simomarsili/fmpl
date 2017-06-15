rm test.txt.scores;
../src/fmpl -i test.txt 2>/dev/null ;
if ! cmp test.txt.scores test.txt.SCORES >/dev/null 2>&1; then
    echo "test failed"
else
    echo "test ok"
fi
