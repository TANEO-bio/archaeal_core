#!usr/bin/env python
import subprocess


subprocess.Popen("mkdir faa2", shell=True)
for i in range(0, 32):
    cmd = "python faa1.py {}".format(i)
    subprocess.Popen(cmd, shell=True)
