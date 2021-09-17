import subprocess
import sys

textfile = sys.argv[1]
inputfile = sys.argv[2]
variable = sys.argv[3]
#Probably have to add inputfile and outputfile names here as well as sys.argv[2], sys.argv[3]
#outputfile = sys.argv[3]

with open(file=textfile, mode="r") as f:
    lines = f.readlines()[1:]
    for i, line in enumerate(lines):
        c1, c2, c3 = line.split(sep="\t")

        print("processing line {}".format(i))

        subprocess.run(
            [
                "ncks",
                "-v",
                "{}".format(variable),
                "-d",
                "longitude,{}".format(float(c3)),
                "-d",
                "latitude,{}".format(float(c2)),
                inputfile,
                "{}".format(variable) + "_1980_2014_{}".format(c1) + ".nc"
                
            ]
        )

