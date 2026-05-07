cat commands.txt | xargs -d'\n' --max-procs=3 -I CMD bash -c CMD
cat commands2.txt | xargs -d'\n' --max-procs=3 -I CMD bash -c CMD
